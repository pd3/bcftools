/*  vcfcall.c -- SNP/indel variant calling from VCF/BCF.

    Copyright (C) 2013-2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <time.h>
#include <zlib.h>
#include <stdarg.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <ctype.h>
#include "bcftools.h"
#include "call.h"
#include "prob1.h"
#include "ploidy.h"

void error(const char *format, ...);

#ifdef _WIN32
#define srand48(x) srand(x)
#define lrand48() rand()
#endif

#define CF_NO_GENO      1
#define CF_INS_MISSED   (1<<1)
#define CF_CCALL        (1<<2)
#define CF_GVCF         (1<<3)
//                      (1<<4)
//                      (1<<5)
#define CF_ACGT_ONLY    (1<<6)
#define CF_QCALL        (1<<7)
#define CF_ADJLD        (1<<8)
#define CF_NO_INDEL     (1<<9)
#define CF_ANNO_MAX     (1<<10)
#define CF_MCALL        (1<<11)
#define CF_PAIRCALL     (1<<12)
#define CF_QCNT         (1<<13)
#define CF_INDEL_ONLY   (1<<14)

typedef struct
{
    int flag;   // combination of CF_* flags above
    int output_type;
    htsFile *bcf_in, *out_fh;
    char *bcf_fname, *output_fname;
    char **samples;             // for subsampling and ploidy
    int nsamples, *samples_map; // mapping from output sample names to original VCF
    char *regions, *targets;    // regions to process
    int regions_is_file, targets_is_file;

    char *samples_fname;
    int samples_is_file;
    int *sample2sex;    // mapping for ploidy. If negative, interpreted as -1*ploidy
    int *sex2ploidy, *sex2ploidy_prev, nsex;
    ploidy_t *ploidy;

    bcf1_t *missed_line;
    call_t aux;     // parameters and temporary data
    gvcf_t gvcf;

    int argc;
    char **argv;

    //  int flag, prior_type, n1, n_sub, *sublist, n_perm;
    //  uint32_t *trio_aux;
    //  char *prior_file, **subsam;
    //  uint8_t *ploidy;
    //  double theta, pref, indel_frac, min_smpl_frac, min_lrt;
    // Permutation tests
    //  int n_perm, *seeds;
    //  double min_perm_p;
    //  void *bed;
}
args_t;

static char **add_sample(void *name2idx, char **lines, int *nlines, int *mlines, char *name, char sex, int *ith)
{
    int ret = khash_str2int_get(name2idx, name, ith);
    if ( ret==0 ) return lines;

    hts_expand(char*,(*nlines+1),*mlines,lines);
    int len = strlen(name);
    lines[*nlines] = (char*) malloc(len+3);
    memcpy(lines[*nlines],name,len);
    lines[*nlines][len]   = ' ';
    lines[*nlines][len+1] = sex;
    lines[*nlines][len+2] = 0;
    *ith = *nlines;
    (*nlines)++;
    khash_str2int_set(name2idx, strdup(name), *ith);
    return lines;
}

static char **parse_ped_samples(call_t *call, char **vals, int nvals)
{
    int i, j, mlines = 0, nlines = 0;
    kstring_t str = {0,0,0};
    void *name2fam = khash_str2int_init(), *name2idx = khash_str2int_init();
    char **lines = NULL;
    for (i=0; i<nvals; i++)
    {
        str.l = 0;
        kputs(vals[i], &str);
        char *col_ends[5], *tmp = str.s;
        j = 0;
        while ( *tmp && j<5 )
        {
            if ( isspace(*tmp) )
            {
                *tmp = 0;
                ++tmp;
                while ( isspace(*tmp) ) tmp++;  // allow multiple spaces
                col_ends[j] = tmp-1;
                j++;
                continue;
            }
            tmp++;
        }
        if ( j!=5 ) break;

        int ifam, ret = khash_str2int_get(name2fam, str.s, &ifam);
        family_t *fam = ret==0 ? &call->fams[ifam] : NULL;
        if ( !fam )
        {
            call->nfams++;
            hts_expand(family_t, call->nfams, call->mfams, call->fams);
            fam = &call->fams[call->nfams-1];
            fam->name = strdup(str.s);
            for (j=0; j<3; j++) fam->sample[j] = -1;
            khash_str2int_set(name2fam, fam->name, call->nfams-1);
        }

        char sex = col_ends[3][1]=='1' ? 'M' : 'F';
        lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[0]+1, sex, &j);
        if ( strcmp(col_ends[1]+1,"0") )    // father
        {
            if ( fam->sample[CHILD]>=0 ) error("Multiple childs in %s [%s,%s]\n", str.s, lines[j],lines[fam->sample[CHILD]]);
            fam->sample[CHILD] = j;
            if ( fam->sample[FATHER]>=0 ) error("Two fathers in %s?\n", str.s);
            lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[1]+1, 'M', &fam->sample[FATHER]);
        }
        if ( strcmp(col_ends[2]+1,"0") )    // mother
        {
            if ( fam->sample[MOTHER]>=0 ) error("Two mothers in %s?\n", str.s);
            lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[2]+1, 'F', &fam->sample[MOTHER]);
        }
    }
    free(str.s);
    khash_str2int_destroy(name2fam);
    khash_str2int_destroy_free(name2idx);

    if ( i!=nvals ) // not a ped file
    {
        if ( i>0 ) error("Could not parse samples, not a PED format.\n");
        return NULL;
    }
    assert( nlines==nvals );
    for (i=0; i<call->nfams; i++)
        assert( call->fams[i].sample[0]>=0 && call->fams[i].sample[1]>=0 && call->fams[i].sample[2]>=0 ); // multiple childs, not a trio

    // for (i=0; i<call->nfams; i++)
    //     fprintf(stderr,"mother=%s, father=%s, child=%s\n", lines[call->fams[i].sample[MOTHER]],lines[call->fams[i].sample[FATHER]],lines[call->fams[i].sample[CHILD]]);

    return lines;
}


/*
 *  Reads sample names and their ploidy (optional) from a file.
 *  Alternatively, if no such file exists, the file name is interpreted
 *  as a comma-separated list of samples. When ploidy is not present,
 *  the default ploidy 2 is assumed.
 */
static void set_samples(args_t *args, const char *fn, int is_file)
{
    int i, nlines;
    char **lines = hts_readlist(fn, is_file, &nlines);
    if ( !lines ) error("Could not read the file: %s\n", fn);

    char **smpls = parse_ped_samples(&args->aux, lines, nlines);
    if ( smpls )
    {
        for (i=0; i<nlines; i++) free(lines[i]);
        free(lines);
        lines = smpls;
    }

    args->samples_map = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr)); // for subsetting
    args->sample2sex  = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr));
    int dflt_sex_id = ploidy_add_sex(args->ploidy, "F");
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++) args->sample2sex[i] = dflt_sex_id;

    int *old2new = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr));
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++) old2new[i] = -1;

    int nsmpl = 0, map_needed = 0;
    for (i=0; i<nlines; i++)
    {
        char *ss = lines[i];
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", lines[i]);
        if ( *ss=='#' ) continue;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        char x = *se, *xptr = se; *se = 0;

        int ismpl = bcf_hdr_id2int(args->aux.hdr, BCF_DT_SAMPLE, ss);
        if ( ismpl < 0 ) { fprintf(stderr,"Warning: No such sample in the VCF: %s\n",ss); continue; }

        ss = se+1;
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) ss = "2";   // default ploidy
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        if ( se==ss ) { *xptr = x; error("Could not parse: \"%s\"\n", lines[i]); }

        if ( ss[1]==0 && (ss[0]=='0' || ss[0]=='1' || ss[0]=='2') )
            args->sample2sex[nsmpl] = -1*(ss[0]-'0');
        else
            args->sample2sex[nsmpl] = ploidy_add_sex(args->ploidy, ss);

        if ( ismpl!=nsmpl ) map_needed = 1;
        args->samples_map[nsmpl] = ismpl;
        old2new[ismpl] = nsmpl;
        nsmpl++;
    }

    for (i=0; i<args->aux.nfams; i++)
    {
        int j, nmiss = 0;
        family_t *fam = &args->aux.fams[i];
        for (j=0; j<3; j++)
        {
            fam->sample[i] = old2new[fam->sample[i]];
            if ( fam->sample[i]<0 ) nmiss++;
        }
        assert( nmiss==0 || nmiss==3 );
    }
    free(old2new);

    if ( !map_needed ) { free(args->samples_map); args->samples_map = NULL; }

    args->nsamples = nsmpl;
    args->samples = lines;
}

static void init_missed_line(args_t *args)
{
    int i;
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++)
    {
        args->aux.gts[i*2]   = bcf_gt_missing;
        args->aux.gts[i*2+1] = bcf_int32_vector_end;
    }
    args->missed_line = bcf_init1();
    bcf_update_genotypes(args->aux.hdr, args->missed_line, args->aux.gts, 2*bcf_hdr_nsamples(args->aux.hdr));
    bcf_float_set_missing(args->missed_line->qual);
}

static void print_missed_line(bcf_sr_regions_t *regs, void *data)
{
    args_t *args = (args_t*) data;
    call_t *call = &args->aux;
    bcf1_t *missed = args->missed_line;

    if ( args->flag & CF_GVCF ) error("todo: Combine --gvcf and --insert-missed\n");

    char *ss = regs->line.s;
    int i = 0;
    while ( i<args->aux.srs->targets_als-1 && *ss )
    {
        if ( *ss=='\t' ) i++;
        ss++;
    }
    if ( !*ss ) error("Could not parse: [%s] (%d)\n", regs->line.s,args->aux.srs->targets_als);

    missed->rid  = bcf_hdr_name2id(call->hdr,regs->seq_names[regs->prev_seq]);
    missed->pos  = regs->start;
    bcf_update_alleles_str(call->hdr, missed,ss);

    bcf_write1(args->out_fh, call->hdr, missed);
}

static void init_data(args_t *args)
{
    args->aux.srs = bcf_sr_init();

    // Open files for input and output, initialize structures
    if ( args->targets )
    {
        if ( bcf_sr_set_targets(args->aux.srs, args->targets, args->targets_is_file, args->aux.flag&CALL_CONSTR_ALLELES ? 3 : 0)<0 )
            error("Failed to read the targets: %s\n", args->targets);

        if ( args->aux.flag&CALL_CONSTR_ALLELES && args->flag&CF_INS_MISSED )
        {
            args->aux.srs->targets->missed_reg_handler = print_missed_line;
            args->aux.srs->targets->missed_reg_data = args;
        }
    }
    if ( args->regions )
    {
        if ( bcf_sr_set_regions(args->aux.srs, args->regions, args->regions_is_file)<0 )
            error("Failed to read the targets: %s\n", args->regions);
    }

    if ( !bcf_sr_add_reader(args->aux.srs, args->bcf_fname) ) error("Failed to open %s: %s\n", args->bcf_fname,bcf_sr_strerror(args->aux.srs->errnum));
    args->aux.hdr = bcf_sr_get_header(args->aux.srs,0);

    int i;
    if ( args->samples_fname )
    {
        set_samples(args, args->samples_fname, args->samples_is_file);
        args->nsex = ploidy_nsex(args->ploidy);
        args->sex2ploidy = (int*) calloc(args->nsex,sizeof(int));
        args->sex2ploidy_prev = (int*) calloc(args->nsex,sizeof(int));
        args->aux.ploidy = (uint8_t*) malloc(args->nsamples);
        for (i=0; i<args->nsamples; i++) args->aux.ploidy[i] = 2;
        for (i=0; i<args->nsex; i++) args->sex2ploidy_prev[i] = 2;
    }

    if ( args->samples_map )
    {
        args->aux.hdr = bcf_hdr_subset(bcf_sr_get_header(args->aux.srs,0), args->nsamples, args->samples, args->samples_map);
        if ( !args->aux.hdr ) error("Error occurred while subsetting samples\n");
        for (i=0; i<args->nsamples; i++)
            if ( args->samples_map[i]<0 ) error("No such sample: %s\n", args->samples[i]);
        if ( !bcf_hdr_nsamples(args->aux.hdr) ) error("No matching sample found\n");
    }
    else
    {
        args->aux.hdr = bcf_hdr_dup(bcf_sr_get_header(args->aux.srs,0));
        for (i=0; i<args->nsamples; i++)
            if ( bcf_hdr_id2int(args->aux.hdr,BCF_DT_SAMPLE,args->samples[i])<0 )
                error("No such sample: %s\n", args->samples[i]);
    }

    args->out_fh = hts_open(args->output_fname, hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));

    if ( args->flag & CF_QCALL )
        return;

    if ( args->flag & CF_MCALL )
        mcall_init(&args->aux);

    if ( args->flag & CF_CCALL )
        ccall_init(&args->aux);

    if ( args->flag&CF_GVCF )
    {
        bcf_hdr_append(args->aux.hdr,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
        args->gvcf.rid  = -1;
        args->gvcf.line = bcf_init1();
        args->gvcf.gt   = (int32_t*) malloc(2*sizeof(int32_t)*bcf_hdr_nsamples(args->aux.hdr));
        for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++)
        {
            args->gvcf.gt[2*i+0] = bcf_gt_unphased(0);
            args->gvcf.gt[2*i+1] = bcf_gt_unphased(0);
        }
    }

    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "QS");
    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "I16");

    bcf_hdr_append_version(args->aux.hdr, args->argc, args->argv, "bcftools_call");
    bcf_hdr_write(args->out_fh, args->aux.hdr);

    if ( args->flag&CF_INS_MISSED ) init_missed_line(args);
}

static void destroy_data(args_t *args)
{
    if ( args->flag & CF_CCALL ) ccall_destroy(&args->aux);
    else if ( args->flag & CF_MCALL ) mcall_destroy(&args->aux);
    else if ( args->flag & CF_QCALL ) qcall_destroy(&args->aux);
    int i;
    for (i=0; i<args->nsamples; i++) free(args->samples[i]);
    if ( args->aux.fams )
    {
        for (i=0; i<args->aux.nfams; i++) free(args->aux.fams[i].name);
        free(args->aux.fams);
    }
    if ( args->missed_line ) bcf_destroy(args->missed_line);
    if ( args->gvcf.line ) bcf_destroy(args->gvcf.line);
    ploidy_destroy(args->ploidy);
    free(args->sex2ploidy);
    free(args->sex2ploidy_prev);
    free(args->gvcf.gt);
    free(args->gvcf.dp);
    free(args->samples);
    free(args->samples_map);
    free(args->sample2sex);
    free(args->aux.ploidy);
    bcf_hdr_destroy(args->aux.hdr);
    hts_close(args->out_fh);
    bcf_sr_destroy(args->aux.srs);
}

void parse_novel_rate(args_t *args, const char *str)
{
    if ( sscanf(str,"%le,%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del,&args->aux.trio_Pm_ins)==3 )  // explicit for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = 1 - args->aux.trio_Pm_del;
        args->aux.trio_Pm_ins  = 1 - args->aux.trio_Pm_ins;
    }
    else if ( sscanf(str,"%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del)==2 )   // dynamic for indels
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_ins  = -1;    // negative value for dynamic calculation
    }
    else if ( sscanf(str,"%le",&args->aux.trio_Pm_SNPs)==1 )  // same for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = -1;
        args->aux.trio_Pm_ins  = -1;
    }
    else error("Could not parse --novel-rate %s\n", str);
}

static int parse_format_flag(const char *str)
{
    int flag = 0;
    const char *ss = str;
    while ( *ss )
    {
        const char *se = ss;
        while ( *se && *se!=',' ) se++;
        if ( !strncasecmp(ss,"GQ",se-ss) ) flag |= CALL_FMT_GQ;
        else if ( !strncasecmp(ss,"GP",se-ss) ) flag |= CALL_FMT_GP;
        else
        {
            fprintf(stderr,"Could not parse \"%s\"\n", str);
            exit(1);
        }
        if ( !*se ) break;
        ss = se + 1;
    }
    return flag;
}

static void set_ploidy(args_t *args, bcf1_t *rec)
{
    ploidy_query(args->ploidy,(char*)bcf_seqname(args->aux.hdr,rec),rec->pos,args->sex2ploidy,NULL,NULL);

    int i;
    for (i=0; i<args->nsex; i++)
        if ( args->sex2ploidy[i]!=args->sex2ploidy_prev[i] ) break;

    if ( i==args->nsex ) return;    // ploidy same as previously

    for (i=0; i<args->nsamples; i++)
    {
        if ( args->sample2sex[i]<0 )
            args->aux.ploidy[i] = -1*args->sample2sex[i];
        else
            args->aux.ploidy[i] = args->sex2ploidy[args->sample2sex[i]];
    }

    int *tmp = args->sex2ploidy; args->sex2ploidy = args->sex2ploidy_prev; args->sex2ploidy_prev = tmp;
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with samtools mpileup.\n");
    fprintf(stderr, "         This command replaces the former \"bcftools view\" caller. Some of the original\n");
    fprintf(stderr, "         functionality has been temporarily lost in the process of transition to htslib,\n");
    fprintf(stderr, "         but will be added back on popular demand. The original calling model can be\n");
    fprintf(stderr, "         invoked with the -c option.\n");
    fprintf(stderr, "Usage:   bcftools call [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "File format options:\n");
    fprintf(stderr, "   -o, --output <file>             write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type <b|u|z|v>     output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "       --ploidy <file>             space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n");
    fprintf(stderr, "   -r, --regions <region>          restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file <file>       restrict to regions listed in a file\n");
    fprintf(stderr, "   -s, --samples <list>            list of samples to include [all samples]\n");
    fprintf(stderr, "   -S, --samples-file <file>       PED file or a file with an optional column with sex (see man page for details) [all samples]\n");
    fprintf(stderr, "   -t, --targets <region>          similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file <file>       similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "   -A, --keep-alts                 keep all possible alternate alleles at variant sites\n");
    fprintf(stderr, "   -f, --format-fields <list>      output format fields: GQ,GP (lowercase allowed) []\n");
    fprintf(stderr, "   -g, --gvcf <minDP>              output gVCF blocks of homozygous REF calls\n");
    fprintf(stderr, "   -i, --insert-missed             output also sites missed by mpileup but present in -T\n");
    fprintf(stderr, "   -M, --keep-masked-ref           keep sites with masked reference allele (REF=N)\n");
    fprintf(stderr, "   -V, --skip-variants <type>      skip indels/snps\n");
    fprintf(stderr, "   -v, --variants-only             output variant sites only\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus/variant calling options:\n");
    fprintf(stderr, "   -c, --consensus-caller          the original calling method (conflicts with -m)\n");
    fprintf(stderr, "   -C, --constrain <str>           one of: alleles, trio (see manual)\n");
    fprintf(stderr, "   -m, --multiallelic-caller       alternative model for multiallelic and rare-variant calling (conflicts with -c)\n");
    fprintf(stderr, "   -n, --novel-rate <float>,[...]  likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]\n");
    fprintf(stderr, "   -p, --pval-threshold <float>    variant if P(ref|D)<FLOAT with -c [0.5]\n");
    fprintf(stderr, "   -P, --prior <float>             mutation rate (use bigger for greater sensitivity) [1.1e-3]\n");

    // todo (and more)
    // fprintf(stderr, "\nContrast calling and association test options:\n");
    // fprintf(stderr, "       -1 INT    number of group-1 samples [0]\n");
    // fprintf(stderr, "       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [%g]\n", args->aux.min_lrt);
    // fprintf(stderr, "       -U INT    number of permutations for association testing (effective with -1) [0]\n");
    // fprintf(stderr, "       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [%g]\n", args->aux.min_perm_p);
    fprintf(stderr, "\n");
    exit(-1);
}

int main_vcfcall(int argc, char *argv[])
{
    char *ploidy_fname = NULL;
    args_t args;
    memset(&args, 0, sizeof(args_t));
    args.argc = argc; args.argv = argv;
    args.aux.prior_type = -1;
    args.aux.indel_frac = -1;
    args.aux.theta      = 1.1e-3;
    args.aux.pref       = 0.5;
    args.aux.min_perm_p = 0.01;
    args.aux.min_lrt    = 1;
    args.flag           = CF_ACGT_ONLY;
    args.output_fname   = "-";
    args.output_type    = FT_VCF;
    args.aux.trio_Pm_SNPs = 1 - 1e-8;
    args.aux.trio_Pm_ins  = args.aux.trio_Pm_del  = 1 - 1e-9;

    int c;
    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"gvcf",1,0,'g'},
        {"format-fields",1,0,'f'},
        {"output",1,0,'o'},
        {"output-type",1,0,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"keep-alts",0,0,'A'},
        {"insert-missed",0,0,'i'},
        {"skip-Ns",0,0,'N'},            // now the new default
        {"keep-masked-refs",0,0,'M'},
        {"skip-variants",1,0,'V'},
        {"variants-only",0,0,'v'},
        {"consensus-caller",0,0,'c'},
        {"constrain",1,0,'C'},
        {"multiallelic-caller",0,0,'m'},
        {"pval-threshold",1,0,'p'},
        {"prior",1,0,'P'},
        {"novel-rate",1,0,'n'},
        {"ploidy",1,0,1},
        {0,0,0,0}
    };

    char *tmp = NULL;
    while ((c = getopt_long(argc, argv, "h?o:O:r:R:s:S:t:T:ANMV:vcmp:C:n:P:f:ig:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : ploidy_fname = optarg; break;
            case 'g':
                args.flag |= CF_GVCF;
                args.gvcf.min_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse, expected integer argument: -g %s\n", optarg);
                break;
            case 'f': args.aux.output_tags |= parse_format_flag(optarg); break;
            case 'M': args.flag &= ~CF_ACGT_ONLY; break;     // keep sites where REF is N
            case 'N': args.flag |= CF_ACGT_ONLY; break;      // omit sites where first base in REF is N (the new default)
            case 'A': args.aux.flag |= CALL_KEEPALT; break;
            case 'c': args.flag |= CF_CCALL; break;          // the original EM based calling method
            case 'i': args.flag |= CF_INS_MISSED; break;
            case 'v': args.aux.flag |= CALL_VARONLY; break;
            case 'o': args.output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args.output_type = FT_BCF_GZ; break;
                          case 'u': args.output_type = FT_BCF; break;
                          case 'z': args.output_type = FT_VCF_GZ; break;
                          case 'v': args.output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'C':
                      if ( !strcasecmp(optarg,"alleles") ) args.aux.flag |= CALL_CONSTR_ALLELES;
                      else if ( !strcasecmp(optarg,"trio") ) args.aux.flag |= CALL_CONSTR_TRIO;
                      else error("Unknown argument to -C: \"%s\"\n", optarg);
                      break;
            case 'V':
                      if ( !strcasecmp(optarg,"snps") ) args.flag |= CF_INDEL_ONLY;
                      else if ( !strcasecmp(optarg,"indels") ) args.flag |= CF_NO_INDEL;
                      else error("Unknown skip category \"%s\" (-S argument must be \"snps\" or \"indels\")\n", optarg);
                      break;
            case 'm': args.flag |= CF_MCALL; break;         // multiallelic calling method
            case 'p':
                args.aux.pref = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --pval-threshold %s\n", optarg);
                break;
            case 'P': args.aux.theta = strtod(optarg,&tmp);
                      if ( *tmp ) error("Could not parse, expected float argument: -P %s\n", optarg);
                      break;
            case 'n': parse_novel_rate(&args,optarg); break;
            case 'r': args.regions = optarg; break;
            case 'R': args.regions = optarg; args.regions_is_file = 1; break;
            case 't': args.targets = optarg; break;
            case 'T': args.targets = optarg; args.targets_is_file = 1; break;
            case 's': args.samples_fname = optarg; break;
            case 'S': args.samples_fname = optarg; args.samples_is_file = 1; break;
            default: usage(&args);
        }
    }
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args.bcf_fname = "-";  // reading from stdin
        else usage(&args);
    }
    else args.bcf_fname = argv[optind++];

    // Sanity check options and initialize
    if ( ploidy_fname ) args.ploidy = ploidy_init(ploidy_fname, 2);
    else
    {
        args.ploidy = ploidy_init_string(
                "X 1 60000 M 1\n"
                "X 2699521 154931043 M 1\n"
                "Y 1 59373566 M 1\n"
                "Y 1 59373566 F 0\n"
                "MT 1 16569 M 1\n"
                "MT 1 16569 F 1\n", 2);
    }
    if ( !args.ploidy ) error("Could not initialize ploidy\n");
    if ( (args.flag & CF_CCALL ? 1 : 0) + (args.flag & CF_MCALL ? 1 : 0) + (args.flag & CF_QCALL ? 1 : 0) > 1 ) error("Only one of -c or -m options can be given\n");
    if ( !(args.flag & CF_CCALL) && !(args.flag & CF_MCALL) && !(args.flag & CF_QCALL) ) error("Expected -c or -m option\n");
    if ( args.aux.n_perm && args.aux.ngrp1_samples<=0 ) error("Expected -1 with -U\n");    // not sure about this, please fix
    if ( args.aux.flag & CALL_CONSTR_ALLELES )
    {
        if ( !args.targets ) error("Expected -t or -T with \"-C alleles\"\n");
        if ( !(args.flag & CF_MCALL) ) error("The \"-C alleles\" mode requires -m\n");
    }
    if ( args.flag & CF_INS_MISSED && !(args.aux.flag&CALL_CONSTR_ALLELES) ) error("The -i option requires -C alleles\n");
    init_data(&args);

    while ( bcf_sr_next_line(args.aux.srs) )
    {
        bcf1_t *bcf_rec = args.aux.srs->readers[0].buffer[0];
        if ( args.samples_map ) bcf_subset(args.aux.hdr, bcf_rec, args.nsamples, args.samples_map);
        bcf_unpack(bcf_rec, BCF_UN_STR);

        // Skip unwanted sites
        int i, is_indel = bcf_is_snp(bcf_rec) ? 0 : 1;
        if ( (args.flag & CF_INDEL_ONLY) && !is_indel ) continue;
        if ( (args.flag & CF_NO_INDEL) && is_indel ) continue;
        if ( (args.flag & CF_ACGT_ONLY) && (bcf_rec->d.allele[0][0]=='N' || bcf_rec->d.allele[0][0]=='n') ) continue;   // REF[0] is 'N'

        // Which allele is symbolic? All SNPs should have it, but not indels
        args.aux.unseen = 0;
        for (i=1; i<bcf_rec->n_allele; i++)
        {
            if ( bcf_rec->d.allele[i][0]=='X' ) { args.aux.unseen = i; break; }  // old X
            if ( bcf_rec->d.allele[i][0]=='<' )
            {
                if ( bcf_rec->d.allele[i][1]=='X' && bcf_rec->d.allele[i][2]=='>' ) { args.aux.unseen = i; break; } // old <X>
                if ( bcf_rec->d.allele[i][1]=='*' && bcf_rec->d.allele[i][2]=='>' ) { args.aux.unseen = i; break; } // new <*>
            }
        }
        int is_ref = (bcf_rec->n_allele==1 || (bcf_rec->n_allele==2 && args.aux.unseen>0)) ? 1 : 0;

        if ( !args.aux.unseen && !is_indel )
        {
            // No symbolic <*> allele, the site cannot be called. Print it as it is or skip if -v was given.
            if ( !(args.aux.flag & CALL_VARONLY) || !is_ref ) bcf_write1(args.out_fh, args.aux.hdr, bcf_rec);
            continue;
        }
        else if ( is_ref )
        {
            if ( args.flag & CF_GVCF )
            {
                bcf_rec = gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, bcf_rec, 1);
                if ( !bcf_rec || args.aux.flag&CALL_VARONLY ) continue;
            }
            else if ( args.aux.flag&CALL_VARONLY )  continue;
        }

        bcf_unpack(bcf_rec, BCF_UN_ALL);
        if ( args.nsex ) set_ploidy(&args, bcf_rec);

        // Various output modes: QCall output (todo)
        if ( args.flag & CF_QCALL )
        {
            qcall(&args.aux, bcf_rec);
            continue;
        }

        // Calling modes which output VCFs
        int ret;
        if ( args.flag & CF_MCALL )
            ret = mcall(&args.aux, bcf_rec);
        else
            ret = ccall(&args.aux, bcf_rec);
        if ( ret==-1 ) error("Something is wrong\n");

        // gVCF output. If is_ref is set, the record has been already processed
        if ( args.flag & CF_GVCF && !is_ref )
        {
            bcf_rec = gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, bcf_rec, ret?0:1);
            if ( !bcf_rec ) continue;
        }

        // Normal output
        if ( (args.aux.flag & CALL_VARONLY) && ret==0 ) continue;     // not a variant
        bcf_write1(args.out_fh, args.aux.hdr, bcf_rec);
    }
    if ( args.flag & CF_GVCF ) gvcf_write(args.out_fh, &args.gvcf, args.aux.hdr, NULL, 0);
    if ( args.flag & CF_INS_MISSED ) bcf_sr_regions_flush(args.aux.srs->targets);
    destroy_data(&args);
    return 0;
}

