/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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

/*
    Trio haplotypes: mother (A,B), father (C,D), child (E,F)
    Modeling the following states:
        01|23|02 
        01|23|03
        01|23|12
        01|23|13
        01|23|20
        01|23|30
        01|23|21
        01|23|31
    with the likelihoods of two haplotypes A,B segments sharing an allele:
        P(01|A==B)  .. e (P of error)
        P(00|A==B)  .. 1-e
    and
        P(ab,cd,ef|E=A,F=C) = P(ea|E=A)*P(fc|F=C)


    Unrelated samples: (A,B) and (C,D)
    Modeling the states:
        xxxx .. A!=C,A!=D,B!=C,B!=D
        0x0x .. A=C,B!=D
        0xx0 .. A=D,B!=C
        x00x .. B=C,A!=D
        x0x0 .. B=D,A!=C
        0101 .. A=C,B=D
        0110 .. A=D,B=C
    with the likelihoods
        P(01|A!=B)  .. f*(1-f)
        P(00|A!=B)  .. (1-f)*(1-f)
        P(11|A!=B)  .. f*f

    Assuming 2x30 crossovers, P=2e-8.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <errno.h>
#include "bcftools.h"
#include "HMM.h"

#define C_TRIO 1
#define C_UNRL 2

// states for unrelated samples
#define UNRL_xxxx  0
#define UNRL_0x0x  1
#define UNRL_0xx0  2
#define UNRL_x00x  3
#define UNRL_x0x0  4
#define UNRL_0101  5
#define UNRL_0110  6

// trio states
#define TRIO_AC 0
#define TRIO_AD 1
#define TRIO_BC 2
#define TRIO_BD 3
#define TRIO_CA 4
#define TRIO_DA 5
#define TRIO_CB 6
#define TRIO_DB 7

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    hmm_t *hmm;
    double *eprob, *tprob, pij, pgt_err;
    uint32_t *sites;
    int32_t *gt_arr;
    int nsites, msites, ngt_arr, prev_rid;
    int mode, nstates;
    int imother,ifather,ichild, isample,jsample;
    void (*set_observed_prob) (bcf1_t *rec);
    char *prefix;
    FILE *fp;
}
args_t;

static args_t args;

static void set_observed_prob_trio(bcf1_t *rec);
static void set_observed_prob_unrelated(bcf1_t *rec);
static void init_hmm_trio(args_t *args);
static void init_hmm_unrelated(args_t *args);


const char *about(void)
{
    return "Color shared chromosomal segments, requires phased GTs.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Color shared chromosomal segments, requires phased GTs.\n"
        "Usage: bcftools +color-chrs [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --prefix <path>     output files prefix\n"
        "   -t, --trio <m,f,c>      names of mother, father and the child\n"
        "   -u, --unrelated <a,b>   names of two unrelated samples\n"
        "\n"
        "Example:\n"
        "   bcftools +color-chrs in.vcf --\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char *trio_samples = NULL, *unrelated_samples = NULL;
    memset(&args,0,sizeof(args_t));
    args.prev_rid = -1;
    args.hdr = in;
    args.pij = 2e-8;
    args.pgt_err = 1e-9;

    static struct option loptions[] =
    {
        {"prefix",1,0,'p'},
        {"trio",1,0,'t'},
        {"unrelated",1,0,'u'},
        {0,0,0,0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "?ht:u:p:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': args.prefix = optarg; break;
            case 't': trio_samples = optarg; break;
            case 'u': unrelated_samples = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc ) error(usage());
    if ( trio_samples && unrelated_samples ) error("Expected only one of the -t/-u options\n");
    if ( !trio_samples && !unrelated_samples ) error("Expected one of the -t/-u options\n");
    if ( !args.prefix ) error("Expected the -p option\n");

    int ret = bcf_hdr_set_samples(args.hdr, trio_samples ? trio_samples : unrelated_samples, 0);
    if ( ret<0 ) error("Could not parse samples: %s\n", trio_samples ? trio_samples : unrelated_samples);
    else if ( ret>0 ) error("%d-th sample not found: %s\n", ret,trio_samples ? trio_samples : unrelated_samples);

    if ( trio_samples )
    {
        int i,n = 0;
        char **list = hts_readlist(trio_samples, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        args.imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        args.ichild  = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[2]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
        args.set_observed_prob = set_observed_prob_trio;
        args.mode = C_TRIO;
        init_hmm_trio(&args);
    }
    else
    {
        int i,n = 0;
        char **list = hts_readlist(unrelated_samples, 0, &n);
        if ( n!=2 ) error("Expected two sample names with -u\n");
        args.isample = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.jsample = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
        args.set_observed_prob = set_observed_prob_unrelated;
        args.mode = C_UNRL;
        init_hmm_unrelated(&args);
    }
    return 1;
}

static void init_hmm_trio(args_t *args)
{
    int i,j;
    args->nstates = 8;
    args->tprob   = (double*) malloc(sizeof(double)*args->nstates*args->nstates);

    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            MAT(args->tprob,args->nstates,i,j) = 0;
    }

    MAT(args->tprob,args->nstates,TRIO_AD,TRIO_AC) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_BC,TRIO_AC) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_BD,TRIO_AC) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_AC,TRIO_AD) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_BC,TRIO_AD) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_BD,TRIO_AD) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_AC,TRIO_BC) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_AD,TRIO_BC) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_BD,TRIO_BC) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_AC,TRIO_BD) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_AD,TRIO_BD) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_BC,TRIO_BD) = args->pij;

    MAT(args->tprob,args->nstates,TRIO_DA,TRIO_CA) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_CB,TRIO_CA) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_DB,TRIO_CA) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_CA,TRIO_DA) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_CB,TRIO_DA) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_DB,TRIO_DA) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_CA,TRIO_CB) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_DA,TRIO_CB) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_DB,TRIO_CB) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_CA,TRIO_DB) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,TRIO_DA,TRIO_DB) = args->pij;
    MAT(args->tprob,args->nstates,TRIO_CB,TRIO_DB) = args->pij;

    for (i=0; i<args->nstates; i++)
    {
        double sum = 0;
        for (j=0; j<args->nstates; j++)
            if ( i!=j ) sum += MAT(args->tprob,args->nstates,i,j);
        MAT(args->tprob,args->nstates,i,i) = 1 - sum;
    }

    #if 0
    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            fprintf(stderr,"\t%e",MAT(args->tprob,args->nstates,j,i));
        fprintf(stderr,"\n");
    }
    #endif

    args->hmm = hmm_init(args->nstates, args->tprob, 10000);
}
static void init_hmm_unrelated(args_t *args)
{
    int i,j;
    args->nstates = 7;
    args->tprob   = (double*) malloc(sizeof(double)*args->nstates*args->nstates);

    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            MAT(args->tprob,args->nstates,i,j) = args->pij;
    }
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_xxxx) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_xxxx) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_x0x0,UNRL_0x0x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_0x0x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_x00x,UNRL_0xx0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_0xx0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_x00x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_x0x0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_0101) = args->pij*args->pij;

    for (i=0; i<args->nstates; i++)
    {
        for (j=i+1; j<args->nstates; j++)
            MAT(args->tprob,args->nstates,i,j) = MAT(args->tprob,args->nstates,j,i);
    }
    for (i=0; i<args->nstates; i++)
    {
        double sum = 0;
        for (j=0; j<args->nstates; j++)
            if ( i!=j ) sum += MAT(args->tprob,args->nstates,i,j);
        MAT(args->tprob,args->nstates,i,i) = 1 - sum;
    }

    #if 0
    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            fprintf(stderr,"\t%e",MAT(args->tprob,args->nstates,j,i));
        fprintf(stderr,"\n");
    }
    #endif

    args->hmm = hmm_init(args->nstates, args->tprob, 10000);
}
static inline double prob_shared(float af, int a, int b)
{
    return a==b ? 1-args.pgt_err : args.pgt_err;
}
static inline double prob_not_shared(float af, int a, int b)
{
    if ( a!=b ) return af*(1-af);
    else if ( a==0 ) return (1-af)*(1-af);
    else return af*af;
}
static void set_observed_prob_unrelated(bcf1_t *rec)
{
    float af = 0.5;  // alternate allele frequency

    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return;
    if ( ngt!=4 ) return;   // chrX

    int32_t a,b,c,d;
    a = args.gt_arr[2*args.isample];
    b = args.gt_arr[2*args.isample+1];
    c = args.gt_arr[2*args.jsample];
    d = args.gt_arr[2*args.jsample+1];
    if ( bcf_gt_is_missing(a) || bcf_gt_is_missing(b) ) return;
    if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) return;
    if ( !bcf_gt_is_phased(a) && !bcf_gt_is_phased(b) ) return; // only the second allele should be set when phased
    if ( !bcf_gt_is_phased(c) && !bcf_gt_is_phased(d) ) return;
    a = bcf_gt_allele(a);
    b = bcf_gt_allele(b);
    c = bcf_gt_allele(c);
    d = bcf_gt_allele(d);

    int m = args.msites;
    args.nsites++;
    hts_expand(uint32_t,args.nsites,args.msites,args.sites);
    if ( m!=args.msites ) args.eprob = (double*) realloc(args.eprob, sizeof(double)*args.msites*args.nstates);

    args.sites[args.nsites-1] = rec->pos;
    double *prob = args.eprob + args.nstates*(args.nsites-1);
    prob[UNRL_xxxx] = prob_not_shared(af,a,c) * prob_not_shared(af,a,d) * prob_not_shared(af,b,c) * prob_not_shared(af,b,d);
    prob[UNRL_0x0x] = prob_shared(af,a,c) * prob_not_shared(af,b,d);
    prob[UNRL_0xx0] = prob_shared(af,a,d) * prob_not_shared(af,b,c);
    prob[UNRL_x00x] = prob_shared(af,b,c) * prob_not_shared(af,a,d);
    prob[UNRL_x0x0] = prob_shared(af,b,d) * prob_not_shared(af,a,c);
    prob[UNRL_0101] = prob_shared(af,a,c) * prob_shared(af,b,d);
    prob[UNRL_0110] = prob_shared(af,a,d) * prob_shared(af,b,c);

#if 0
    static int x = 0;
    if ( !x++)
    {
        printf("p(0==0) .. %f\n", prob_shared(af,0,0));
        printf("p(0!=0) .. %f\n", prob_not_shared(af,0,0));
        printf("p(0==1) .. %f\n", prob_shared(af,0,1));
        printf("p(0!=1) .. %f\n", prob_not_shared(af,0,1));
    }
    printf("%d|%d %d|%d  x:%f 11:%f 12:%f 21:%f 22:%f 11,22:%f 12,21:%f  %d\n", a,b,c,d,
            prob[UNRL_xxxx], prob[UNRL_0x0x], prob[UNRL_0xx0], prob[UNRL_x00x], prob[UNRL_x0x0], prob[UNRL_0101], prob[UNRL_0110], rec->pos+1);
#endif
}
static void set_observed_prob_trio(bcf1_t *rec)
{
    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return;
    if ( ngt!=6 ) return;   // chrX

    int32_t a,b,c,d,e,f;
    a = args.gt_arr[2*args.imother];
    b = args.gt_arr[2*args.imother+1];
    c = args.gt_arr[2*args.ifather];
    d = args.gt_arr[2*args.ifather+1];
    e = args.gt_arr[2*args.ichild];
    f = args.gt_arr[2*args.ichild+1];
    if ( bcf_gt_is_missing(a) || bcf_gt_is_missing(b) ) return;
    if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) return;
    if ( bcf_gt_is_missing(e) || bcf_gt_is_missing(f) ) return;
    if ( !bcf_gt_is_phased(a) && !bcf_gt_is_phased(b) ) return; // only the second allele should be set when phased
    if ( !bcf_gt_is_phased(c) && !bcf_gt_is_phased(d) ) return;
    if ( !bcf_gt_is_phased(e) && !bcf_gt_is_phased(f) ) return;
    a = bcf_gt_allele(a);
    b = bcf_gt_allele(b);
    c = bcf_gt_allele(c);
    d = bcf_gt_allele(d);
    e = bcf_gt_allele(e);
    f = bcf_gt_allele(f);

    int m = args.msites;
    args.nsites++;
    hts_expand(uint32_t,args.nsites,args.msites,args.sites);
    if ( m!=args.msites ) args.eprob = (double*) realloc(args.eprob, sizeof(double)*args.msites*args.nstates);

    args.sites[args.nsites-1] = rec->pos;
    double *prob = args.eprob + args.nstates*(args.nsites-1);
    prob[TRIO_AC] = prob_shared(0,e,a) * prob_shared(0,f,c);
    prob[TRIO_AD] = prob_shared(0,e,a) * prob_shared(0,f,d);
    prob[TRIO_BC] = prob_shared(0,e,b) * prob_shared(0,f,c);
    prob[TRIO_BD] = prob_shared(0,e,b) * prob_shared(0,f,d);
    prob[TRIO_CA] = prob_shared(0,e,c) * prob_shared(0,f,a);
    prob[TRIO_DA] = prob_shared(0,e,d) * prob_shared(0,f,a);
    prob[TRIO_CB] = prob_shared(0,e,c) * prob_shared(0,f,b);
    prob[TRIO_DB] = prob_shared(0,e,d) * prob_shared(0,f,b);
}

void flush_viterbi(args_t *args)
{
    const char *s1, *s2, *s3 = NULL;
    if ( args->mode==C_UNRL )
    {
        s1 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->isample);
        s2 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->jsample);
    }
    else if ( args->mode==C_TRIO )
    {
        s1 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->imother);
        s3 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->ifather);
        s2 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->ichild);
    }

    if ( !args->fp )
    {
        kstring_t str = {0,0,0};
        kputs(args->prefix, &str);
        kputs(".dat", &str);
        args->fp = fopen(str.s,"w");
        if ( !args->fp ) error("%s: %s\n", str.s,strerror(errno));
        free(str.s);
        fprintf(args->fp,"# SG, shared segment\t[2]Chromosome\t[3]Start\t[4]End\t[5]%s:1\t[6]%s:2\n",s2,s2);
    }

    hmm_run_viterbi(args->hmm,args->nsites,args->eprob,args->sites);
    uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
    int i, iprev = -1, prev_state = -1, nstates = hmm_get_nstates(args->hmm);
    for (i=0; i<args->nsites; i++)
    {
        int state = vpath[i*nstates];
        if ( state!=prev_state || i+1==args->nsites )
        {
            uint32_t start = iprev>=0 ? args->sites[iprev]+1 : 1, end = i>0 ? args->sites[i-1] : 1;
            const char *chr = bcf_hdr_id2name(args->hdr,args->prev_rid);
            if ( args->mode==C_UNRL )
            {
                switch (prev_state)
                {
                    case UNRL_0x0x:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t-\n", chr,start,end,s1); break;
                    case UNRL_0xx0:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t-\t%s:1\n", chr,start,end,s1); break;
                    case UNRL_x00x:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t-\n", chr,start,end,s1); break;
                    case UNRL_x0x0:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t-\t%s:2\n", chr,start,end,s1); break;
                    case UNRL_0101:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s1,s1); break;
                    case UNRL_0110:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s1,s1); break;
                }
            }
            else if ( args->mode==C_TRIO )
            {
                switch (prev_state)
                {
                    case TRIO_AC:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr,start,end,s1,s3); break;
                    case TRIO_AD:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s1,s3); break;
                    case TRIO_BC:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s1,s3); break;
                    case TRIO_BD:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr,start,end,s1,s3); break;
                    case TRIO_CA:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr,start,end,s3,s1); break;
                    case TRIO_DA:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s3,s1); break;
                    case TRIO_CB:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s3,s1); break;
                    case TRIO_DB:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr,start,end,s3,s1); break;
                }
            }
            iprev = i-1;
        }
        prev_state = state;
    }
    args->nsites = 0;
}
    
bcf1_t *process(bcf1_t *rec)
{
    if ( args.prev_rid==-1 ) args.prev_rid = rec->rid;
    if ( args.prev_rid!=rec->rid ) flush_viterbi(&args);
    args.prev_rid = rec->rid;
    args.set_observed_prob(rec);
    return NULL;
}

void create_plots(args_t *args)
{
    fclose(args->fp);

    kstring_t str = {0,0,0};
    kputs(args->prefix, &str);
    kputs(".py", &str);
    FILE *fp = fopen(str.s,"w");
    if ( !fp ) error("%s: %s\n", str.s,strerror(errno));
    fprintf(fp,
            "import matplotlib as mpl\n"
            "mpl.use('Agg')\n"
            "import matplotlib.pyplot as plt\n"
            "import csv, os, urllib, gzip, itertools, sys\n"
            "import matplotlib.patches as ptch\n"
            "\n"
            "csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)\n"
            "\n"
            "# Download and parse karyotype file\n"
            "ktype_fname  = 'karyotype.gz'\n"
            "url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz'\n"
            "if not os.path.exists(ktype_fname):\n"
            "    print 'Downloading %%s to local file: %%s' %% (url, ktype_fname)\n"
            "    with open(ktype_fname, 'w') as k_file:\n"
            "        f = urllib.urlopen(url)\n"
            "        k_file.write(f.read())\n"
            "\n"
            "colors = itertools.cycle(['r','g','b','c','y','m','k'])\n"
            "chrs = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']\n"
            "max_len = 0\n"
            "offset  = 0\n"
            "ktype = {}\n"
            "with gzip.open(ktype_fname, 'rb') as f:\n"
            "    reader = csv.reader(f, 'tab')\n"
            "    for row in reader:\n"
            "        chr = row[0]\n"
            "        if chr[0:3]=='chr': chr = chr[3:]\n"
            "        if chr not in ktype:\n"
            "            ktype[chr] = []\n"
            "            offset = float(row[1])\n"
            "        ktype[chr].append([float(row[1])-offset,float(row[2])-offset,row[4]])\n"
            "        if max_len < float(row[2])-offset: max_len = float(row[2])-offset\n"
            "    f.close()\n"
            "\n"
            "def get_chr_length(ktype, chr):\n"
            "    return max([x[1] for x in ktype[chr]])\n"
            "\n"
            "def get_chr_centromere(ktype, chr):\n"
            "    acen = []\n"
            "    for x in ktype[chr]:\n"
            "        if x[2]=='acen':\n"
            "            acen.append(x[0])\n"
            "            acen.append(x[1])\n"
            "    return [acen[0],acen[1],acen[3]]\n"
            "\n"
            "def read_dat(fname, dat, hap):\n"
            "    with open(fname, 'rb') as f:\n"
            "        reader = csv.reader(f, 'tab')\n"
            "        for row in reader:\n"
            "            if row[0]=='SG':\n"
            "                chr   = row[1]\n"
            "                start = float(row[2])\n"
            "                end   = float(row[3])\n"
            "                hap1  = row[4]\n"
            "                hap2  = row[5]\n"
            "                if hap1 not in hap: \n"
            "                    smpl = hap1[2:] \n"
            "                    hap['1:'+smpl] = colors.next()\n"
            "                    hap['2:'+smpl] = colors.next()\n"
            "                    if smpl not in dat: dat[smpl] = {}\n"
            "                    dat[smpl]['1'] = []\n"
            "                    dat[smpl]['2'] = []\n"
            "                smpl2 = hap2[2:]\n"
            "                hap2  = hap2[0:1] \n"
            "                if smpl2 not in dat:\n"
            "                    dat[smpl2] = {}\n"
            "                    dat[smpl2]['1'] = []\n"
            "                    dat[smpl2]['2'] = []\n"
            "                if hap2 not in dat[smpl2]: dat[smpl2][hap2] = []\n"
            "                dat[smpl2][hap2].append([chr,start,end,hap[hap1]])\n"
            "        f.close()\n"
            "\n"
            "dat = {}\n"
            "hap_color = {}\n"
            "read_dat('%s.dat',dat,hap_color)\n"
            "for i in range(1,len(sys.argv)): read_dat(sys.argv[i],dat,hap_color)\n"
            "\n"
            "n = float(len(chrs))\n"
            "width = 0.2/n\n"
            "pad   = 0.7*width\n"
            "\n"
            "def mirror_coors(hap,x,coor):\n"
            "    if hap==1: return coor\n"
            "    xm = x + 2*width + pad\n"
            "    for id in range(len(coor)): coor[id] = (xm - coor[id][0] + x,coor[id][1])\n"
            "    return coor\n"
            "\n"
            "def chr2polygon(ax,chr,hap,x,color='gray',reg=None):\n"
            "    ystart = 0\n"
            "    yend   = get_chr_length(ktype,chr)/max_len\n"
            "    if reg:\n"
            "        ystart = reg[1]/max_len\n"
            "        yend   = reg[2]/max_len\n"
            "    cmer = get_chr_centromere(ktype,chr)\n"
            "    yc1 = cmer[0]/max_len\n"
            "    yc2 = cmer[1]/max_len\n"
            "    yc3 = cmer[2]/max_len\n"
            "    args = {'facecolor':color,'edgecolor':'none','aa':True,'lod':True,'lw':0.01}\n"
            "    if hap==1:\n"
            "        rectangle = ptch.Rectangle([x,yc1],width+0.5*pad,yc3-yc1,**args)\n"
            "    else:\n"
            "        rectangle = ptch.Rectangle([x+width+0.5*pad,yc1],width+0.5*pad,yc3-yc1,**args)\n"
            "    ax.add_patch(rectangle)\n"
            "    if ystart<=yc1:\n"
            "        coor = []\n"
            "        coor.append([x,ystart])\n"
            "        coor.append([x+width,ystart])\n"
            "        y = yend\n"
            "        if y>yc1: y = yc1\n"
            "        coor.append([x+width,y])\n"
            "        coor.append([x,y])\n"
            "        coor = mirror_coors(hap,x,coor)\n"
            "        polygon = ptch.Polygon(coor, **args)\n"
            "        ax.add_patch(polygon)\n"
            "    if yend>=yc3:\n"
            "        coor = []\n"
            "        y = ystart\n"
            "        if y<yc3: y = yc3\n"
            "        coor.append([x,y])\n"
            "        coor.append([x+width,y])\n"
            "        coor.append([x+width,yend])\n"
            "        coor.append([x,yend])\n"
            "        coor.append([x,y])\n"
            "        coor = mirror_coors(hap,x,coor)\n"
            "        polygon = ptch.Polygon(coor, **args)\n"
            "        ax.add_patch(polygon)\n"
            "\n"
            "\n"
            "for smpl in dat:\n"
            "    hap1 = '1:'+smpl\n"
            "    hap2 = '2:'+smpl\n"
            "    fig,ax = plt.subplots(1,1,figsize=(10,5))\n"
            "    for i in range(len(chrs)):\n"
            "        chr = chrs[i]\n"
            "        x = i/n\n"
            "        col = 'gray'\n"
            "        if hap1 in hap_color: col = hap_color[hap1]\n"
            "        chr2polygon(ax,chr,1,x,col)\n"
            "        for sg in dat[smpl]['1']:\n"
            "            if sg[0]==chr:\n"
            "                chr2polygon(ax,chr,1,x,sg[3],sg)\n"
            "        col = 'gray'\n"
            "        if hap2 in hap_color: col = hap_color[hap2]\n"
            "        chr2polygon(ax,chr,2,x,col)\n"
            "        for sg in dat[smpl]['2']:\n"
            "            if sg[0]==chr:\n"
            "                chr2polygon(ax,chr,2,x,sg[3],sg)\n"
            "        ax.annotate(chrs[i],(x+width+0.5*pad,-2*width),fontsize='x-small',va='top',ha='center')\n"
            "    ax.set_xlim(-width,1)\n"
            "    ax.set_ylim(-0.1,1+width)\n"
            "    ax.axes.get_yaxis().set_visible(False)\n"
            "    ax.axes.get_xaxis().set_visible(False)\n"
            "    ax.set_yticks([])\n"
            "    ax.set_frame_on(False)\n"
            "    plt.subplots_adjust(left=0,right=1,top=1,bottom=0)\n"
            "    plt.savefig('%s-'+smpl+'.png',dpi=100)\n"
            "    plt.close()\n",
            args->prefix,args->prefix
    );
    fclose(fp);

return;
    str.l = 0;
    ksprintf(&str,"python %s.py", args->prefix);
    int ret = system(str.s);
    if ( ret) fprintf(stderr, "The command returned non-zero status %d: %s\n", ret, str.s);
    free(str.s);
}

void destroy(void)
{
    flush_viterbi(&args);
    create_plots(&args);
    free(args.gt_arr);
    free(args.tprob);
    free(args.sites);
    free(args.eprob);
    hmm_destroy(args.hmm);
}



