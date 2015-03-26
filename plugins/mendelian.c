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

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <errno.h>
#include "bcftools.h"

#define MODE_COUNT     1
#define MODE_LIST_GOOD 2
#define MODE_LIST_BAD  4

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    int32_t *gt_arr;
    int mode;
    int ngt_arr;
    int nok, nbad, nrec;
    int imother,ifather,ichild;
}
args_t;

static args_t args;

const char *about(void)
{
    return "Count Mendelian consistent / inconsistent genotypes.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Count Mendelian consistent / inconsistent genotypes.\n"
        "Usage: bcftools +mendelian [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -c, --count             count the number of consistent sites\n"
        "   -l, --list [+x]         list consistent (+) or inconsistent (x) sites\n"
        "   -t, --trio <m,f,c>      names of mother, father and the child\n"
        "\n"
        "Example:\n"
        "   bcftools +mendelian in.vcf --\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char *trio_samples = NULL;
    memset(&args,0,sizeof(args_t));
    args.hdr  = in;
    args.mode = 0;

    static struct option loptions[] =
    {
        {"trio",1,0,'t'},
        {"list",1,0,'l'},
        {"count",0,0,'c'},
        {0,0,0,0}
    };
    char c;
    while ((c = getopt_long(argc, argv, "?ht:l:c",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'c': args.mode = MODE_COUNT; break;
            case 'l': 
                if ( !strcmp("+",optarg) ) args.mode = MODE_LIST_GOOD; 
                else if ( !strcmp("x",optarg) ) args.mode = MODE_LIST_BAD; 
                else error("The argument not recognised: --list %s\n", optarg);
                break;
            case 't': trio_samples = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc ) error(usage());
    if ( !trio_samples ) error("Expected the -t option\n");
    if ( !args.mode ) error("Expected one of the -c or -l options\n");

    int ret = bcf_hdr_set_samples(args.hdr, trio_samples, 0);
    if ( ret<0 ) error("Could not parse samples: %s\n", trio_samples);
    else if ( ret>0 ) error("No such sample: %s\n", argv[optind+ret-1]);

    int i,n = 0;
    char **list = hts_readlist(trio_samples, 0, &n);
    if ( n!=3 ) error("Expected three sample names with -t\n");
    args.imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
    args.ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
    args.ichild  = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[2]);
    for (i=0; i<n; i++) free(list[i]);
    free(list);

    return args.mode&(MODE_LIST_GOOD|MODE_LIST_BAD) ? 0 : 1;
}

bcf1_t *process(bcf1_t *rec)
{
    int mother,father,child;
    int32_t a,b,c,d,e,f;

    args.nrec++;

    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) goto not_applicable;
    if ( ngt!=6 ) goto not_applicable;   // chrX

    a = args.gt_arr[2*args.imother];
    b = args.gt_arr[2*args.imother+1];
    c = args.gt_arr[2*args.ifather];
    d = args.gt_arr[2*args.ifather+1];
    e = args.gt_arr[2*args.ichild];
    f = args.gt_arr[2*args.ichild+1];
    if ( bcf_gt_is_missing(a) || bcf_gt_is_missing(b) ) goto not_applicable;
    if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) goto not_applicable;
    if ( bcf_gt_is_missing(e) || bcf_gt_is_missing(f) ) goto not_applicable;

    mother = (1<<bcf_gt_allele(a)) | (1<<bcf_gt_allele(b));
    father = (1<<bcf_gt_allele(c)) | (1<<bcf_gt_allele(d));
    child  = (1<<bcf_gt_allele(e)) | (1<<bcf_gt_allele(f));

    if ( (mother&child) && (father&child) ) 
    {
        args.nok++;
        if ( args.mode&MODE_LIST_GOOD ) return rec;
    }
    else
    {
        args.nbad++; 
        if ( args.mode&MODE_LIST_BAD ) return rec;
    }
    return NULL;

not_applicable:
    return args.mode&MODE_LIST_GOOD ? rec : NULL;
}

void destroy(void)
{
    fprintf(stderr,"nok/nbad/nskipped:\t%d / %d / %d\n", args.nok,args.nbad,args.nrec-(args.nok+args.nbad));
    free(args.gt_arr);
}



