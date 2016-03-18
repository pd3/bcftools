/* 
    Copyright (C) 2014 Genome Research Ltd.

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
    THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdint.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <inttypes.h>
#include <unistd.h>
#include "bcftools.h"

#define GUESS_GT 1
#define GUESS_PL 2
#define GUESS_GL 3

typedef struct
{
    uint64_t ncount;
    double phap, pdip;
}
count_t;

typedef struct
{
    char *chr;
    uint32_t start, end;
    count_t *counts;    // per-sample counts: counts[isample]
}
stats_t;

typedef struct
{
    int argc;
    char **argv, *af_tag;
    double af_dflt;
    stats_t stats;
    int nsample, verbose, tag, include_indels;
    int *counts, ncounts;       // number of observed GTs with given ploidy, used when -g is not given
    double *tmpf, *pl2p, gt_err_prob;
    float *af;
    int maf;
    int32_t *arr, narr;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
}
args_t;

const char *about(void)
{
    return "Determine sample sex by checking genotype likelihoods in haploid regions.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Determine sample sex by checking genotype likelihoods (GL,PL) or genotypes (GT)\n"
        "       in the non-PAR region of chrX. The HWE is assumed, so given the alternate allele\n" 
        "       frequency fA and the genotype likelihoods pRR,pRA,pAA, the probabilities are\n"
        "       calculated as\n"
        "           P(dip) = pRR*(1-fA)^2 + pAA*fA^2 + 2*pRA*(1-fA)*fA\n"
        "           P(hap) = pRR*(1-fA) + pAA*fA\n"
        "       When genotype likelihoods are not available, the -e option is used to account\n"
        "       for genotyping errors with -t GT. The alternate allele frequency fA is estimated\n"
        "       directly from the data (the default) or can be provided by an INFO tag.\n"
        "       The results can be visualized using the accompanied guess-ploidy.py script.\n"
        "       Note that this plugin is intended to replace the former vcf2sex plugin.\n"
        "\n"
        "Usage: bcftools +guess-ploidy <file.vcf.gz> [Plugin Options]\n"
        "Plugin options:\n"
        "       --AF-dflt <float>           the default alternate allele frequency [0.5]\n"
        "       --AF-tag <TAG>              use TAG for allele frequency\n"
        "   -e, --error-rate <float>        probability of GT being wrong (with -t GT) [1e-3]\n"
        "   -i, --include-indels            do not skip indel sites\n"
        "   -r, --regions <chr:beg-end>     [X:2699521-154931043]\n"
        "   -R, --regions-file <file>       regions listed in a file\n"
        "   -t, --tag <tag>                 genotype or genotype likelihoods: GT, PL, GL [PL]\n"
        "   -v, --verbose                   verbose output (specify twice to increase verbosity)\n"
        "\n"
        "Examples:\n"
        "   bcftools +guess-ploidy in.vcf.gz\n"
        "   bcftools +guess-ploidy in.vcf.gz -t GL -r chrX:2699521-154931043\n"
        "   bcftools view file.vcf.gz -r chrX:2699521-154931043 | bcftools +guess-ploidy\n"
        "   bcftools +guess-ploidy in.bcf -v > ploidy.txt && guess-ploidy.py ploidy.txt img\n"
        "\n";
}

void process_region_guess(args_t *args)
{
    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( !args->include_indels && !(bcf_get_variant_types(rec)&VCF_SNP) )  continue;

        double freq[2] = {0,0}, sum;
        int ismpl,i;
        if ( args->tag & GUESS_GT )   // use GTs to guess the ploidy, considering only one ALT
        {
            int ngt = bcf_get_genotypes(args->hdr,rec,&args->arr,&args->narr);
            if ( ngt<=0 ) continue;
            ngt /= args->nsample;
            for (ismpl=0; ismpl<args->nsample; ismpl++)
            {
                int32_t *ptr = args->arr + ismpl*ngt;
                double *tmp = args->tmpf + ismpl*3;

                if ( ptr[0]==bcf_gt_missing ) 
                {
                    tmp[0] = -1;
                    continue;
                }
                if ( ptr[1]==bcf_int32_vector_end )
                {
                    if ( bcf_gt_allele(ptr[0])==0 ) // haploid R
                    {
                        tmp[0] = 1 - 2*args->gt_err_prob;
                        tmp[1] = tmp[2] = args->gt_err_prob;
                    }
                    else    // haploid A
                    {
                        tmp[0] = tmp[1] = args->gt_err_prob;
                        tmp[2] = 1 - 2*args->gt_err_prob;
                    }
                    continue;
                }
                if ( bcf_gt_allele(ptr[0])==0 && bcf_gt_allele(ptr[1])==0 ) // RR
                {
                    tmp[0] = 1 - 2*args->gt_err_prob;
                    tmp[1] = tmp[2] = args->gt_err_prob;
                }
                else if ( bcf_gt_allele(ptr[0])==bcf_gt_allele(ptr[1]) ) // AA
                {
                    tmp[0] = tmp[1] = args->gt_err_prob;
                    tmp[2] = 1 - 2*args->gt_err_prob;
                }
                else  // RA or hetAA, treating as RA
                {
                    tmp[1] = 1 - 2*args->gt_err_prob;
                    tmp[0] = tmp[2] = args->gt_err_prob;
                }
                freq[0] += 2*tmp[0]+tmp[1];
                freq[1] += tmp[1]+2*tmp[2];
            }
        }
        else
        {
            // use PL or GL to guess the ploidy, restrict to first ALT allele
            int npl = bcf_get_format_int32(args->hdr,rec,args->tag&GUESS_PL?"PL":"GL",&args->arr,&args->narr);
            if ( npl<=0 ) continue;
            npl /= args->nsample;
            for (ismpl=0; ismpl<args->nsample; ismpl++)
            {
                int32_t *ptr = args->arr + ismpl*npl;
                double *tmp = args->tmpf + ismpl*3;

                // restrict to first ALT
                if ( ptr[0]==bcf_int32_missing || ptr[1]==bcf_int32_missing || ptr[2]==bcf_int32_missing ) 
                {
                    tmp[0] = -1;
                    continue;
                }
                if ( ptr[0]==ptr[1] && ptr[0]==ptr[2] )
                {
                    tmp[0] = -1;
                    continue;
                }
                if ( args->tag&GUESS_PL )
                {
                    for (i=0; i<3; i++)
                        tmp[i] = (ptr[i]<0 || ptr[i]>=256) ? args->pl2p[255] : args->pl2p[ptr[i]];
                }
                else
                    for (i=0; i<3; i++) tmp[i] = pow(10.,ptr[i]);   // GL

                sum = 0;
                for (i=0; i<3; i++) sum += tmp[i];
                for (i=0; i<3; i++) tmp[i] /= sum;
                freq[0] += 2*tmp[0]+tmp[1];
                freq[1] += tmp[1]+2*tmp[2];
            }
        }
        if ( args->af_tag )
        {
            int ret = bcf_get_info_float(args->hdr,rec,args->af_tag,&args->af, &args->maf);
            if ( ret>0 ) { freq[0] = 1 - args->af[0]; freq[1] = args->af[0]; }
        }

        if ( !freq[0] && !freq[1] ) { freq[0] = 1 - args->af_dflt; freq[1] = args->af_dflt; }
        sum = freq[0] + freq[1];
        freq[0] /= sum;
        freq[1] /= sum;
        for (ismpl=0; ismpl<args->nsample; ismpl++)
        {
            count_t *counts = &args->stats.counts[ismpl];
            double *tmp = args->tmpf + ismpl*3;
            if ( tmp[0] < 0 ) continue;
            double phap = freq[0]*tmp[0] + freq[1]*tmp[2];
            double pdip = freq[0]*freq[0]*tmp[0] + 2*freq[0]*freq[1]*tmp[1] + freq[1]*freq[1]*tmp[2];
            counts->phap += log(phap);
            counts->pdip += log(pdip);
            counts->ncount++;
            if ( args->verbose>1 )
                printf("DBG\t%s\t%d\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n", bcf_seqname(args->hdr,rec),rec->pos+1,bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,ismpl),
                    freq[1],tmp[0],tmp[1],tmp[2],phap,pdip);
        }
    }
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->tag    = GUESS_PL;
    args->argc   = argc; args->argv = argv;
    args->gt_err_prob = 1e-3;
    args->af_dflt = 0.5;
    char *region  = "X:2699521-154931043";
    int region_is_file = 0;
    static struct option loptions[] =
    {
        {"AF-tag",1,0,0},
        {"AF-dflt",1,0,1},
        {"verbose",0,0,'v'},
        {"include-indels",0,0,'i'},
        {"error-rate",1,0,'e'},
        {"tag",1,0,'t'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"background",1,0,'b'},
        {0,0,0,0}
    };
    char c, *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:t:e:i",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 0: args->af_tag = optarg; break;
            case 1: 
                    args->af_dflt = strtod(optarg,&tmp);
                    if ( *tmp ) error("Could not parse: --AF-dflt %s\n", optarg);
                    break;
            case 'i': args->include_indels = 1; break;
            case 'e':
                args->gt_err_prob = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -e %s\n", optarg);
                if ( args->gt_err_prob<0 || args->gt_err_prob>1 ) error("Expected value from the interval [0,1]: -e %s\n", optarg);
                break;
            case 'R': region_is_file = 1; break; 
            case 'r': region = optarg; break; 
            case 'v': args->verbose++; break; 
            case 't':
                if ( !strcasecmp(optarg,"GT") ) args->tag = GUESS_GT;
                else if ( !strcasecmp(optarg,"PL") ) args->tag = GUESS_PL;
                else if ( !strcasecmp(optarg,"GL") ) args->tag = GUESS_GL;
                else error("The argument not recognised, expected --tag GT, PL or GL: %s\n", optarg);
                break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    char *fname = NULL;
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else { error(usage()); }
    }
    else if ( optind+1!=argc ) error(usage());
    else fname = argv[optind];

    args->sr = bcf_sr_init();
    if ( strcmp("-",fname) )
    {
        args->sr->require_index = 1;
        if ( region )
        {
            if ( bcf_sr_set_regions(args->sr, region, region_is_file)<0 )
                error("Failed to read the regions: %s\n",region);
        }
    }
    if ( !bcf_sr_add_reader(args->sr,fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = args->sr->readers[0].header;
    args->nsample = bcf_hdr_nsamples(args->hdr);
    args->stats.counts = (count_t*) calloc(args->nsample,sizeof(count_t));

    if ( args->af_tag && !bcf_hdr_idinfo_exists(args->hdr,BCF_HL_INFO,bcf_hdr_id2int(args->hdr,BCF_DT_ID,args->af_tag)) )
        error("No such INFO tag: %s\n", args->af_tag);

    int i;
    if ( args->tag&GUESS_PL )
    {
        args->pl2p = (double*) calloc(256,sizeof(double));
        for (i=0; i<256; i++) args->pl2p[i] = pow(10., -i/10.);
    }
    if ( args->tag&GUESS_PL || args->tag&GUESS_GL )
        args->tmpf = (double*) malloc(sizeof(*args->tmpf)*3*args->nsample);

    if ( args->verbose )
    {
        printf("# This file was produced by: bcftools +guess-ploidy(%s+htslib-%s)\n", bcftools_version(),hts_version());
        printf("# The command line was:\tbcftools +%s", args->argv[0]);
        for (i=1; i<args->argc; i++)
            printf(" %s",args->argv[i]);
        printf("\n");
        printf("# [1]SEX\t[2]Sample\t[3]Predicted sex\t[4]log P(Haploid)/nSites\t[5]log P(Diploid)/nSites\t[6]nSites\t[7]Score: F < 0 < M ($4-$5)\n");
        if ( args->verbose>1 )
            printf("# [1]DBG\t[2]Chr\t[3]Pos\t[4]Sample\t[5]AF\t[6]pRR\t[7]pRA\t[8]pAA\t[9]P(Haploid)\t[10]P(Diploid)\n");
    }

    process_region_guess(args);

    for (i=0; i<args->nsample; i++)
    {
        double phap = args->stats.counts[i].ncount ? args->stats.counts[i].phap / args->stats.counts[i].ncount : 0.5;
        double pdip = args->stats.counts[i].ncount ? args->stats.counts[i].pdip / args->stats.counts[i].ncount : 0.5;
        if ( args->verbose )
        {
            printf("SEX\t%s\t%s\t%f\t%f\t%"PRId64"\t%f\n", args->hdr->samples[i],phap>pdip?"M":"F",
                    phap,pdip,args->stats.counts[i].ncount,phap-pdip);
        }
        else
            printf("%s\t%s\n", args->hdr->samples[i],phap>pdip?"M":"F");
    }
   
    bcf_sr_destroy(args->sr);
    free(args->pl2p);
    free(args->tmpf);
    free(args->counts);
    free(args->stats.counts);
    free(args->arr);
    free(args->af);
    free(args);
    return 0;
}



