/*  bam_plcmd.c -- mpileup subcommand.

    Copyright (C) 2008-2014 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <htslib/regidx.h>
#include "bcftools.h"
#include "gvcf.h"

static inline int printw(int c, FILE *fp)
{
    char buf[16];
    int l, x;
    if (c == 0) return fputc('0', fp);
    for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
    if (c < 0) buf[l++] = '-';
    buf[l] = 0;
    for (x = 0; x < l/2; ++x) {
        int y = buf[x]; buf[x] = buf[l-1-x]; buf[l-1-x] = y;
    }
    fputs(buf, fp);
    return 0;
}

#include <assert.h>
#include "bam2bcf.h"
#include "sample.h"

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)

typedef struct {
    int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter, output_type;
    int openQ, extQ, tandemQ, min_support; // for indels
    double min_frac; // for indels
    char *reg, *pl_list, *fai_fname, *output_fname, *samples_list;
    faidx_t *fai;
    void *rghash_exc;
    regidx_t *bed;
    gvcf_t *gvcf;
    int argc;
    char **argv;
} mplp_conf_t;

typedef struct {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

typedef struct {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
    void *rghash_inc;   // whitelist for readgroups on per-bam basis
} mplp_aux_t;

typedef struct {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
} mplp_pileup_t;

static int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    if (r->ref[1])
        free(r->ref[1]);
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                0x7fffffff,
                                &r->ref_len[0]);

    if (!r->ref) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

static int mplp_func(void *data, bam1_t *b)
{
    char *ref;
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, skip = 0, ref_len;
    do {
        int has_ref;
        ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
            skip = 1;
            continue;
        }
        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
        if (ma->conf->bed) { // test overlap
            skip = regidx_overlap(ma->conf->bed, ma->h->target_name[b->core.tid],b->core.pos, bam_endpos(b)-1, NULL) ? 0 : 1;
            if (skip) continue;
        }
        if ( ma->rghash_inc )
        {
            // skip unless read group is listed
            uint8_t *rg = bam_aux_get(b, "RG");
            if ( !rg ) { skip = 1; continue; }
            if ( !khash_str2int_has_key(ma->rghash_inc, (const char*)(rg+1)) ) { skip = 1; continue; }
        }
        if (ma->conf->rghash_exc) { // exclude read groups
            uint8_t *rg = bam_aux_get(b, "RG");
            skip = (rg && khash_str2int_get(ma->conf->rghash_exc, (const char*)(rg+1), NULL)==0);
            if (skip) continue;
        }
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }
        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
                        __func__, b->core.pos, ref_len, b->core.tid);
                skip = 1;
                continue;
            }
        } else {
            has_ref = 0;
        }
        skip = 0;
        if (has_ref && (ma->conf->flag&MPLP_REALN)) bam_prob_realn_core(b, ref, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = bam_cap_mapQ(b, ref, ma->conf->capQ_thres);
            if (q < 0) skip = 1;
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) skip = 1;
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;
    } while (skip);
    return ret;
}

static void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
                       int n, char *const*fn, int *n_plp, const bam_pileup1_t **plp, int ignore_rg)
{
    int i, j;
    memset(m->n_plp, 0, m->n * sizeof(int));
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n_plp[i]; ++j) {
            const bam_pileup1_t *p = plp[i] + j;
            uint8_t *q;
            int id = -1;
            q = ignore_rg? 0 : bam_aux_get(p->b, "RG");
            if (q) id = bam_smpl_rg2smid(sm, fn[i], (char*)q+1, buf);
            if (id < 0) id = bam_smpl_rg2smid(sm, fn[i], 0, buf);
            if (id < 0 || id >= m->n) {
                assert(q); // otherwise a bug
                fprintf(stderr, "[%s] Read group %s used in file %s but absent from the header or an alignment missing read group.\n", __func__, (char*)q+1, fn[i]);
                exit(1);
            }
            if (m->n_plp[id] == m->m_plp[id]) {
                m->m_plp[id] = m->m_plp[id]? m->m_plp[id]<<1 : 8;
                m->plp[id] = (bam_pileup1_t*) realloc(m->plp[id], sizeof(bam_pileup1_t) * m->m_plp[id]);
            }
            m->plp[id][m->n_plp[id]++] = *p;
        }
    }
}

static void flush_bcf_records(mplp_conf_t *conf, htsFile *fp, bcf_hdr_t *hdr, bcf1_t *rec)
{
    if ( !conf->gvcf )
    {
        if ( rec ) bcf_write1(fp, hdr, rec);
        return;
    }

    if ( !rec )
    {
        gvcf_write(conf->gvcf, fp, hdr, NULL, 0);
        return;
    }

    int is_ref = 0;
    if ( rec->n_allele==1 ) is_ref = 1;
    else if ( rec->n_allele==2 )
    {
        // second allele is mpileup's X, not a variant
        if ( rec->d.allele[1][0]=='<' && rec->d.allele[1][1]=='*' && rec->d.allele[1][2]=='>' ) is_ref = 1;
    }
    rec = gvcf_write(conf->gvcf, fp, hdr, rec, is_ref);
    if ( rec ) bcf_write1(fp,hdr,rec);
}

/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 */
static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
    extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
    extern void bcf_call_del_rghash(void *rghash);
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, beg0 = 0, end0 = 1u<<29, ref_len, max_depth, max_indel_depth;
    const bam_pileup1_t **plp;
    mplp_ref_t mp_ref = MPLP_REF_INIT;
    bam_mplp_t iter;
    bam_hdr_t *h = NULL; /* header of first file in input list */
    char *ref;
    void *rghash_exc = NULL;

    bcf_callaux_t *bca = NULL;
    bcf_callret1_t *bcr = NULL;
    bcf_call_t bc;
    htsFile *bcf_fp = NULL;
    bcf_hdr_t *bcf_hdr = NULL;

    bam_sample_t *sm = NULL;
    kstring_t buf;
    mplp_pileup_t gplp;

    memset(&gplp, 0, sizeof(mplp_pileup_t));
    memset(&buf, 0, sizeof(kstring_t));
    memset(&bc, 0, sizeof(bcf_call_t));
    data = (mplp_aux_t**) calloc(n, sizeof(mplp_aux_t*));
    plp = (const bam_pileup1_t**) calloc(n, sizeof(bam_pileup1_t*));
    n_plp = (int*) calloc(n, sizeof(int));
    sm = bam_smpl_init();

    if (n == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(1);
    }

    void *samples_inc = NULL;
    if ( conf->samples_list )
    {
        int nsamples = 0;
        char **samples = hts_readlist(conf->samples_list, 0, &nsamples);
        if ( nsamples ) samples_inc = khash_str2int_init();
        for (i=0; i<nsamples; i++)
        {
            khash_str2int_inc(samples_inc, strdup(samples[i]));
            free(samples[i]);
        }
        free(samples);
    }

    // read the header of each file in the list and initialize data
    for (i = 0; i < n; ++i) {
        bam_hdr_t *h_tmp;
        data[i] = (mplp_aux_t*) calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open(fn[i], "rb");
        if ( !data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
            exit(1);
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(1);
        }
        if (conf->fai_fname && hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(1);
        }
        data[i]->conf = conf;
        data[i]->ref = &mp_ref;
        h_tmp = sam_hdr_read(data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
            exit(1);
        }
        data[i]->h = i? h : h_tmp; // for i==0, "h" has not been set yet
        if ( samples_inc ) data[i]->rghash_inc = khash_str2int_init();
        bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text, samples_inc, data[i]->rghash_inc);
        // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
        rghash_exc = bcf_call_add_rg(rghash_exc, h_tmp->text, conf->pl_list);
        if (conf->reg) {
            hts_idx_t *idx = sam_index_load(data[i]->fp, fn[i]);
            if (idx == 0) {
                fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
                exit(1);
            }
            if ( (data[i]->iter=sam_itr_querys(idx, data[i]->h, conf->reg)) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, conf->reg);
                exit(1);
            }
            if (i == 0) beg0 = data[i]->iter->beg, end0 = data[i]->iter->end;
            hts_idx_destroy(idx);
        }
        else
            data[i]->iter = NULL;

        if (i == 0) h = h_tmp; /* save the header of first file in list */
        else {
            // FIXME: to check consistency  and h_tmp
            bam_hdr_destroy(h_tmp);
            
            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            data[i]->h = h;
        }
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    gplp.n = sm->n;
    gplp.n_plp = (int*) calloc(sm->n, sizeof(int));
    gplp.m_plp = (int*) calloc(sm->n, sizeof(int));
    gplp.plp = (bam_pileup1_t**) calloc(sm->n, sizeof(bam_pileup1_t*));

    fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
    // write the VCF header
    bcf_fp = hts_open(conf->output_fname?conf->output_fname:"-", hts_bcf_wmode(conf->output_type));
    if (bcf_fp == NULL) {
        fprintf(stderr, "[%s] failed to write to %s: %s\n", __func__, conf->output_fname? conf->output_fname : "standard output", strerror(errno));
        exit(1);
    }

    bcf_hdr = bcf_hdr_init("w");
    kstring_t str = {0,0,0};

    ksprintf(&str, "##mpileupVersion=%s+htslib-%s\n",bcftools_version(),hts_version());
    bcf_hdr_append(bcf_hdr, str.s);

    str.l = 0;
    ksprintf(&str, "##mpileupCommand=bcftools mpileup");
    for (i=1; i<conf->argc; i++) ksprintf(&str, " %s", conf->argv[i]);
    kputc('\n', &str);
    bcf_hdr_append(bcf_hdr, str.s);

    if (conf->fai_fname)
    {
        str.l = 0;
        ksprintf(&str, "##reference=file://%s\n", conf->fai_fname);
        bcf_hdr_append(bcf_hdr, str.s);
    }

    // todo: use/write new BAM header manipulation routines, fill also UR, M5
    for (i=0; i<h->n_targets; i++)
    {
        str.l = 0;
        ksprintf(&str, "##contig=<ID=%s,length=%d>", h->target_name[i], h->target_len[i]);
        bcf_hdr_append(bcf_hdr, str.s);
    }
    free(str.s);
    bcf_hdr_append(bcf_hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">");
#if CDF_MWU_TESTS
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");
#endif
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
    bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if ( conf->fmt_flag&B2B_FMT_DP )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
    if ( conf->fmt_flag&B2B_FMT_DV )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
    if ( conf->fmt_flag&B2B_FMT_DPR )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_INFO_DPR )
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_FMT_DP4 )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
    if ( conf->fmt_flag&B2B_FMT_SP )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
    if ( conf->fmt_flag&B2B_FMT_AD )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");
    if ( conf->fmt_flag&B2B_FMT_ADF )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">");
    if ( conf->fmt_flag&B2B_FMT_ADR )
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
    bcf_hdr_append(bcf_hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
    if ( conf->gvcf )
        gvcf_update_header(conf->gvcf, bcf_hdr);

    for (i=0; i<sm->n; i++)
        bcf_hdr_add_sample(bcf_hdr, sm->smpl[i]);
    bcf_hdr_add_sample(bcf_hdr, NULL);
    bcf_hdr_write(bcf_fp, bcf_hdr);

    bca = bcf_call_init(-1., conf->min_baseQ);
    bcr = (bcf_callret1_t*) calloc(sm->n, sizeof(bcf_callret1_t));
    bca->openQ = conf->openQ, bca->extQ = conf->extQ, bca->tandemQ = conf->tandemQ;
    bca->min_frac = conf->min_frac;
    bca->min_support = conf->min_support;
    bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;

    bc.bcf_hdr = bcf_hdr;
    bc.n = sm->n;
    bc.PL = (int32_t*) malloc(15 * sm->n * sizeof(*bc.PL));
    if (conf->fmt_flag)
    {
        assert( sizeof(float)==sizeof(int32_t) );
        bc.DP4 = (int32_t*) malloc(sm->n * sizeof(int32_t) * 4);
        bc.fmt_arr = (uint8_t*) malloc(sm->n * sizeof(float)); // all fmt_flag fields, float and int32
        if ( conf->fmt_flag&(B2B_INFO_DPR|B2B_FMT_DPR|B2B_INFO_AD|B2B_INFO_ADF|B2B_INFO_ADR|B2B_FMT_AD|B2B_FMT_ADF|B2B_FMT_ADR) )
        {
            // first B2B_MAX_ALLELES fields for total numbers, the rest per-sample
            bc.ADR = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            bc.ADF = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            for (i=0; i<sm->n; i++)
            {
                bcr[i].ADR = bc.ADR + (i+1)*B2B_MAX_ALLELES;
                bcr[i].ADF = bc.ADF + (i+1)*B2B_MAX_ALLELES;
            }
        }
    }

    // begin pileup
    iter = bam_mplp_init(n, mplp_func, (void**)data);
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
    max_depth = conf->max_depth;
    if (max_depth * sm->n > 1<<20)
        fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
    if (max_depth * sm->n < 8000) {
        max_depth = 8000 / sm->n;
        fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
    }
    max_indel_depth = conf->max_indel_depth * sm->n;
    bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;
    while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
        if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        if (conf->bed && tid >= 0 && !regidx_overlap(conf->bed, h->target_name[tid], pos, pos, NULL)) continue;
        mplp_get_ref(data[0], tid, &ref, &ref_len);
        int total_depth, _ref0, ref16;
        for (i = total_depth = 0; i < n; ++i) total_depth += n_plp[i];
        group_smpl(&gplp, sm, &buf, n, fn, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
        _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
        ref16 = seq_nt16_table[_ref0];
        bcf_callaux_clean(bca, &bc);
        for (i = 0; i < gplp.n; ++i)
            bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
        bc.tid = tid; bc.pos = pos;
        bcf_call_combine(gplp.n, bcr, bca, ref16, &bc);
        bcf_clear1(bcf_rec);
        bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, 0, 0);
        flush_bcf_records(conf, bcf_fp, bcf_hdr, bcf_rec);
        // call indels; todo: subsampling with total_depth>max_indel_depth instead of ignoring?
        // check me: rghash in bcf_call_gap_prep() should have no effect, reads mplp_func already excludes them
        if (!(conf->flag&MPLP_NO_INDEL) && total_depth < max_indel_depth && bcf_call_gap_prep(gplp.n, gplp.n_plp, gplp.plp, pos, bca, ref, NULL) >= 0)
        {
            bcf_callaux_clean(bca, &bc);
            for (i = 0; i < gplp.n; ++i)
                bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], -1, bca, bcr + i);
            if (bcf_call_combine(gplp.n, bcr, bca, -1, &bc) >= 0) {
                bcf_clear1(bcf_rec);
                bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, bca, ref);
                flush_bcf_records(conf, bcf_fp, bcf_hdr, bcf_rec);
            }
        }
    }
    flush_bcf_records(conf, bcf_fp, bcf_hdr, NULL);

    // clean up
    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
    if (bcf_fp)
    {
        hts_close(bcf_fp);
        bcf_hdr_destroy(bcf_hdr);
        bcf_call_destroy(bca);
        free(bc.PL);
        free(bc.DP4);
        free(bc.ADR);
        free(bc.ADF);
        free(bc.fmt_arr);
        free(bcr);
    }
    if ( conf->gvcf ) gvcf_destroy(conf->gvcf);
    bam_smpl_destroy(sm); free(buf.s);
    for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
    free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
    bcf_call_del_rghash(rghash_exc);
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);
    for (i = 0; i < n; ++i) {
        sam_close(data[i]->fp);
        if (data[i]->iter) hts_itr_destroy(data[i]->iter);
        if ( data[i]->rghash_inc ) khash_str2int_destroy_free(data[i]->rghash_inc);
        free(data[i]);
    }
    free(data); free(plp); free(n_plp);
    if (mp_ref.ref[0]) free(mp_ref.ref[0]);
    if (mp_ref.ref[1]) free(mp_ref.ref[1]);
    if ( samples_inc ) khash_str2int_destroy_free(samples_inc);
    return ret;
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = (char**) calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) )
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (stat(buf, &sb) != 0)
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; }
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = (char**) realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    fclose(fh);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}
#undef MAX_PATH_LEN

int parse_format_flag(const char *str)
{
    int i, flag = 0, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        if ( !strcasecmp(tags[i],"DP") ) flag |= B2B_FMT_DP;
        else if ( !strcasecmp(tags[i],"DV") ) flag |= B2B_FMT_DV;
        else if ( !strcasecmp(tags[i],"SP") ) flag |= B2B_FMT_SP;
        else if ( !strcasecmp(tags[i],"DP4") ) flag |= B2B_FMT_DP4;
        else if ( !strcasecmp(tags[i],"DPR") ) flag |= B2B_FMT_DPR;
        else if ( !strcasecmp(tags[i],"INFO/DPR") ) flag |= B2B_INFO_DPR;
        else if ( !strcasecmp(tags[i],"AD") ) flag |= B2B_FMT_AD;
        else if ( !strcasecmp(tags[i],"ADF") ) flag |= B2B_FMT_ADF;
        else if ( !strcasecmp(tags[i],"ADR") ) flag |= B2B_FMT_ADR;
        else if ( !strcasecmp(tags[i],"INFO/AD") ) flag |= B2B_INFO_AD;
        else if ( !strcasecmp(tags[i],"INFO/ADF") ) flag |= B2B_INFO_ADF;
        else if ( !strcasecmp(tags[i],"INFO/ADR") ) flag |= B2B_INFO_ADR;
        else
        {
            fprintf(stderr,"Could not parse tag \"%s\" in \"%s\"\n", tags[i], str);
            exit(1);
        }
        free(tags[i]);
    }
    if (n_tags) free(tags);
    return flag;
}

static void print_usage(FILE *fp, const mplp_conf_t *mplp)
{
    char *tmp_require = bam_flag2str(mplp->rflag_require);
    char *tmp_filter  = bam_flag2str(mplp->rflag_filter);

    fprintf(fp,
            "\n"
            "About: This is the original \"samtools mpileup\" command which was transferred\n"
            "       to bcftools in order to avoid errors resulting from use of incompatible\n"
            "       versions of samtools and bcftools.\n"
            "Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]\n"
            "\n"
            "Input options:\n"
            "  -6, --illumina1.3+      quality is in the Illumina-1.3+ encoding\n"
            "  -A, --count-orphans     do not discard anomalous read pairs\n"
            "  -b, --bam-list FILE     list of input BAM filenames, one per line\n"
            "  -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)\n"
            "  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]\n"
            "  -d, --max-depth INT     max per-BAM depth; avoids excessive memory usage [%d]\n", mplp->max_depth);
    fprintf(fp,
            "  -E, --redo-BAQ          recalculate BAQ on the fly, ignore existing BQs\n"
            "  -f, --fasta-ref FILE    faidx indexed reference sequence file\n"
            "  -G, --exclude-RG FILE   exclude read groups listed in FILE\n"
            "  -l, --positions FILE    skip unlisted positions (chr pos) or regions (BED)\n"
            "  -q, --min-MQ INT        skip alignments with mapQ smaller than INT [%d]\n", mplp->min_mq);
    fprintf(fp,
            "  -Q, --min-BQ INT        skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp->min_baseQ);
    fprintf(fp,
            "  -r, --region REG        region in which pileup is generated\n"
            "  -R, --ignore-RG         ignore RG tags (one BAM = one sample)\n"
            "  -s, --samples <list>    list of samples to include\n"
            "  --rf, --incl-flags STR|INT  required flags: skip reads with mask bits unset [%s]\n", tmp_require);
    fprintf(fp,
            "  --ff, --excl-flags STR|INT  filter flags: skip reads with mask bits set\n"
            "                                            [%s]\n", tmp_filter);
    fprintf(fp,
            "  -x, --ignore-overlaps   disable read-pair overlap detection\n"
            "\n"
            "Output options:\n"
            "      --gvcf INT[,...]    group non-variant sites into gVCF blocks according to minimum per-sample DP\n"
            "  -o, --output FILE       write output to FILE [standard output]\n"
            "  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n"
            "  -t, --output-tags LIST  optional tags to output:\n"
            "               DP,DPR,DV,DP4,INFO/DPR,SP,AD,ADF,ADR,INFO/AD,INFO/ADF,INFO/ADR []\n"
            "\n"
            "SNP/INDEL genotype likelihoods options (effective with -g/-v):\n"
            "  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [%d]\n", mplp->extQ);
    fprintf(fp,
            "  -F, --gap-frac FLOAT    minimum fraction of gapped reads [%g]\n", mplp->min_frac);
    fprintf(fp,
            "  -h, --tandem-qual INT   coefficient for homopolymer errors [%d]\n", mplp->tandemQ);
    fprintf(fp,
            "  -I, --skip-indels       do not perform indel calling\n"
            "  -L, --max-idepth INT    maximum per-sample depth for INDEL calling [%d]\n", mplp->max_indel_depth);
    fprintf(fp,
            "  -m, --min-ireads INT    minimum number gapped reads for indel candidates [%d]\n", mplp->min_support);
    fprintf(fp,
            "  -o, --open-prob INT     Phred-scaled gap open seq error probability [%d]\n", mplp->openQ);
    fprintf(fp,
            "  -p, --per-sample-mF     apply -m and -F per-sample for increased sensitivity\n"
            "  -P, --platforms STR     comma separated list of platforms for indels [all]\n"
            "\n"
            "Notes: Assuming diploid individuals.\n");

    free(tmp_require);
    free(tmp_filter);
}

int bam_mpileup(int argc, char *argv[])
{
    int c;
    const char *file_list = NULL;
    char **fn = NULL;
    int nfiles = 0, use_orphan = 0;
    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.min_baseQ = 13;
    mplp.capQ_thres = 0;
    mplp.max_depth = 250; mplp.max_indel_depth = 250;
    mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
    mplp.min_frac = 0.002; mplp.min_support = 1;
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    mplp.argc = argc; mplp.argv = argv;
    mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    mplp.output_type = FT_VCF;
    static const struct option lopts[] =
    {
        {"rf", required_argument, NULL, 1},   // require flag
        {"ff", required_argument, NULL, 2},   // filter flag
        {"incl-flags", required_argument, NULL, 1},
        {"excl-flags", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"open-prob", required_argument, NULL, 4},
        {"gvcf", required_argument, NULL, 5},
        {"illumina1.3+", no_argument, NULL, '6'},
        {"count-orphans", no_argument, NULL, 'A'},
        {"bam-list", required_argument, NULL, 'b'},
        {"no-BAQ", no_argument, NULL, 'B'},
        {"no-baq", no_argument, NULL, 'B'},
        {"adjust-MQ", required_argument, NULL, 'C'},
        {"adjust-mq", required_argument, NULL, 'C'},
        {"max-depth", required_argument, NULL, 'd'},
        {"redo-BAQ", no_argument, NULL, 'E'},
        {"redo-baq", no_argument, NULL, 'E'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"exclude-RG", required_argument, NULL, 'G'},
        {"exclude-rg", required_argument, NULL, 'G'},
        {"positions", required_argument, NULL, 'l'},
        {"region", required_argument, NULL, 'r'},
        {"ignore-RG", no_argument, NULL, 'R'},
        {"ignore-rg", no_argument, NULL, 'R'},
        {"samples", required_argument, NULL, 's'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"ignore-overlaps", no_argument, NULL, 'x'},
        {"BCF", no_argument, NULL, 'g'},
        {"bcf", no_argument, NULL, 'g'},
        {"VCF", no_argument, NULL, 'v'},
        {"vcf", no_argument, NULL, 'v'},
        {"output-type", no_argument, NULL, 'O'},
        {"output-bp", no_argument, NULL, 'O'},
        {"output-tags", required_argument, NULL, 't'},
        {"uncompressed", no_argument, NULL, 'u'},
        {"ext-prob", required_argument, NULL, 'e'},
        {"gap-frac", required_argument, NULL, 'F'},
        {"tandem-qual", required_argument, NULL, 'h'},
        {"skip-indels", no_argument, NULL, 'I'},
        {"max-idepth", required_argument, NULL, 'L'},
        {"min-ireads ", required_argument, NULL, 'm'},
        {"per-sample-mF", no_argument, NULL, 'p'},
        {"per-sample-mf", no_argument, NULL, 'p'},
        {"platforms", required_argument, NULL, 'P'},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "Agf:r:l:q:Q:uRC:BDSd:L:b:P:po:e:h:Im:F:EG:6O:s:Vvxt:",lopts,NULL)) >= 0) {
        switch (c) {
        case 'x': mplp.flag &= ~MPLP_SMART_OVERLAPS; break;
        case  1 :
            mplp.rflag_require = bam_str2flag(optarg);
            if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse --rf %s\n", optarg); return 1; }
            break;
        case  2 :
            mplp.rflag_filter = bam_str2flag(optarg);
            if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --ff %s\n", optarg); return 1; }
            break;
        case  3 : mplp.output_fname = optarg; break;
        case  4 : mplp.openQ = atoi(optarg); break;
        case  5 :
                mplp.gvcf = gvcf_init(optarg);
                if ( !mplp.gvcf ) error("Could not parse: --gvcf %s\n", optarg);
                break;
        case 'f':
            mplp.fai = fai_load(optarg);
            if (mplp.fai == 0) return 1;
            mplp.fai_fname = optarg;
            break;
        case 'd': mplp.max_depth = atoi(optarg); break;
        case 'r': mplp.reg = strdup(optarg); break;
        case 'l':
                  // In the original version the whole BAM was streamed which is inefficient
                  //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine
                  //  best strategy, that is streaming or jumping.
                  mplp.bed = regidx_init(optarg,NULL,NULL,0,NULL);
                  if (!mplp.bed) { fprintf(stderr, "Could not read file: \"%s\"\n", optarg); return 1; }
                  break;
        case 'P': mplp.pl_list = strdup(optarg); break;
        case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
        case 'g': fprintf(stderr,"[warning] bcftools mpileup `-g` is functional, but deprecated. Please swith to -O in future.\n"); mplp.flag |= MPLP_BCF; break;
        case 'v': fprintf(stderr,"[warning] bcftools mpileup `-v` is functional, but deprecated. Please swith to -O in future.\n"); mplp.flag |= MPLP_BCF | MPLP_VCF; break;
        case 'u': fprintf(stderr,"[warning] bcftools mpileup `-u` is functional, but deprecated. Please swith to -O in future.\n"); mplp.flag |= MPLP_NO_COMP | MPLP_BCF; break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'D': mplp.fmt_flag |= B2B_FMT_DP; fprintf(stderr, "[warning] bcftools mpileup option `-D` is functional, but deprecated. Please switch to `-t DP` in future.\n"); break;
        case 'S': mplp.fmt_flag |= B2B_FMT_SP; fprintf(stderr, "[warning] bcftools mpileup option `-S` is functional, but deprecated. Please switch to `-t SP` in future.\n"); break;
        case 'V': mplp.fmt_flag |= B2B_FMT_DV; fprintf(stderr, "[warning] bcftools mpileup option `-V` is functional, but deprecated. Please switch to `-t DV` in future.\n"); break;
        case 'I': mplp.flag |= MPLP_NO_INDEL; break;
        case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
        case '6': mplp.flag |= MPLP_ILLUMINA13; break;
        case 'R': mplp.flag |= MPLP_IGNORE_RG; break;
        case 's': mplp.samples_list = optarg; break;
        case 'O': 
                  switch (optarg[0]) {
                      case 'b': mplp.output_type = FT_BCF_GZ; break;
                      case 'u': mplp.output_type = FT_BCF; break;
                      case 'z': mplp.output_type = FT_VCF_GZ; break;
                      case 'v': mplp.output_type = FT_VCF; break;
                      default: error("[error] The option \"-O\" changed meaning, which one do you want: \"bcftools mpileup --output-type\" or \"samtools mpileup --output-BP\"?\n", optarg); 
                  }
            break;
        case 'C': mplp.capQ_thres = atoi(optarg); break;
        case 'q': mplp.min_mq = atoi(optarg); break;
        case 'Q': mplp.min_baseQ = atoi(optarg); break;
        case 'b': file_list = optarg; break;
        case 'o': {
                char *end;
                long value = strtol(optarg, &end, 10);
                // Distinguish between -o INT and -o FILE (a bit of a hack!)
                if (*end == '\0') mplp.openQ = value;
                else mplp.output_fname = optarg;
            }
            break;
        case 'e': mplp.extQ = atoi(optarg); break;
        case 'h': mplp.tandemQ = atoi(optarg); break;
        case 'A': use_orphan = 1; break;
        case 'F': mplp.min_frac = atof(optarg); break;
        case 'm': mplp.min_support = atoi(optarg); break;
        case 'L': mplp.max_indel_depth = atoi(optarg); break;
        case 'G': {
                FILE *fp_rg;
                char buf[1024];
                mplp.rghash_exc = khash_str2int_init();
                if ((fp_rg = fopen(optarg, "r")) == 0)
                    fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
                while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
                    khash_str2int_inc(mplp.rghash_exc, strdup(buf));
                fclose(fp_rg);
            }
            break;
        case 't': mplp.fmt_flag |= parse_format_flag(optarg); break;
        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            return 1;
        }
    }
    if ( mplp.gvcf && !(mplp.fmt_flag&B2B_FMT_DP) ) mplp.fmt_flag |= B2B_FMT_DP;
    if ( mplp.flag&(MPLP_BCF|MPLP_VCF|MPLP_NO_COMP) )
    {
        if ( mplp.flag&MPLP_VCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_VCF;
            else mplp.output_type = FT_VCF_GZ;
        }
        else if ( mplp.flag&MPLP_BCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_BCF;
            else mplp.output_type = FT_BCF_GZ;
        }
    }
    if ( !(mplp.flag&MPLP_REALN) && mplp.flag&MPLP_REDO_BAQ )
    {
        fprintf(stderr,"Error: The -B option cannot be combined with -E\n");
        return 1;
    }
    if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
    if (argc == 1)
    {
        print_usage(stderr, &mplp);
        return 1;
    }
    int ret;
    if (file_list) {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        ret = mpileup(&mplp,nfiles,fn);
        for (c=0; c<nfiles; c++) free(fn[c]);
        free(fn);
    }
    else
        ret = mpileup(&mplp, argc - optind, argv + optind);
    if (mplp.rghash_exc) khash_str2int_destroy_free(mplp.rghash_exc);
    free(mplp.reg); free(mplp.pl_list);
    if (mplp.fai) fai_destroy(mplp.fai);
    if (mplp.bed) regidx_destroy(mplp.bed);
    return ret;
}
