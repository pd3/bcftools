#!/usr/bin/env python
#
# Plot the output of "bcftools +guess-ploidy -v"
#

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random, math, sys
import csv
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)

if len(sys.argv) != 3:
    print >> sys.stderr, 'Usage: guess-ploidy.py <guess-ploidy.out> <image-prefix>'
    sys.exit()

prefix = sys.argv[2]
dpi = 150

def add_value(dat,key,x,y):
    if key not in dat:
        dat[key] = []
    dat[key].append([x,y])

smpl2sex = {}
dat = {}
with open(sys.argv[1]) as f:
    reader = csv.reader(f, 'tab')
    for row in reader:
        if row[0][0]=="#": continue
        if row[0]=="SEX":
            smpl  = row[1]
            sex   = row[2]
            phap  = float(row[3])
            pdip  = float(row[4])
            ndat  = float(row[5])
            score = float(row[6])
            smpl2sex[smpl] = sex
            add_value(dat,'score',smpl,score)
            add_value(dat,'phap',smpl,phap)
            add_value(dat,'pdip',smpl,pdip)
            add_value(dat,'ndat',smpl,ndat)

def sort_by_val(arr):
    for x in (sorted(arr, key=lambda x:x[1])):
        id = len(smpl2id)
        smpl2id[x[0]] = id
    arr = sorted(arr, key=lambda x:smpl2id[x[0]])
    return arr

def select_sex(arr,sex):
    out = []
    for x in arr:
        if smpl2sex[x[0]]==sex: out.append(x)
    return out

col = {}
col['blue']   = '#396ab1' 
col['orange'] = '#da7c30'
col['green']  = '#3e9651'
col['red']    = '#cc2529'
#col['grey']   = '#535154'   # too light
col['grey']   = '#3e3e3e'
col['purple'] = '#6b4c9a'
col['brown']  = '#922428'
#col['yellow'] = '#948b3d'   # too dark
col['yellow'] = '#ccc210'

#col['blue']   = '#7293cb'      # light colors for bars
#col['orange'] = '#e1974c'
#col['green']  = '#84ba5b'
#col['red']    = '#d35e60'
#col['grey']   = '#808585'
#col['purple'] = '#9067a7'
#col['brown']  = '#ab6857'
#col['yellow'] = '#ccc210'


if True:
    fig,ax1 = plt.subplots(1,1,figsize=(6,4))
    smpl2id = {}
    dat['score']  = sort_by_val(dat['score'])
    dat['scoreM'] = select_sex(dat['score'],'M')
    dat['scoreF'] = select_sex(dat['score'],'F')
    ax2 = ax1.twinx()
    plots  = ax1.plot([smpl2id[x[0]] for x in dat['phap']],[x[1] for x in dat['phap']],'.',color=col['blue'],ms=3,label='log P(haploid)')
    plots += ax1.plot([smpl2id[x[0]] for x in dat['pdip']],[x[1] for x in dat['pdip']],'.',color=col['yellow'],ms=3,label='log P(diploid)')
    plots += ax1.plot([smpl2id[x[0]] for x in dat['scoreM']],[x[1] for x in dat['scoreM']],'.',color=col['green'],label='Total score: Males')
    plots += ax1.plot([smpl2id[x[0]] for x in dat['scoreF']],[x[1] for x in dat['scoreF']],'.',color=col['red'],label='Total score: Females')
    plots += ax2.plot([smpl2id[x[0]] for x in dat['ndat']],[x[1] for x in dat['ndat']],'.',color=col['grey'],ms=3,label='Number of sites')
    labels = [l.get_label() for l in plots]
    ax1.legend(plots,labels,loc='best', frameon=False, numpoints=1, prop={'size':10})
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.patch.set_visible(False)
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Score')
    ax2.set_ylabel('Number of sites')
    ax2.set_yscale('log')
    #   ax1.set_yscale('log')
    ax1.ticklabel_format(style='sci', scilimits=(-3,4), axis='x')
    plt.subplots_adjust(left=0.13,right=0.89,bottom=0.13,top=0.9,hspace=0.1)
    plt.savefig(prefix+'.png',dpi=dpi)

plt.close()



