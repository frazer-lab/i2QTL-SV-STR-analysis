
# coding: utf-8

# In[2]:

from __future__ import division
import numpy as np
import os
import sys
import datetime
from subprocess import call
import subprocess
import glob
import djPyBio as DJ
from djPyBio import Common as CM
from djPyBio import mpltools as axtools


import pandas as pd
import csv
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy 
import pybedtools as pbt
import ciepy
import cardipspy as cpy
import itertools
import tempfile
import six
import networkx as nx
import scipy.stats as stats
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import  make_axes_locatable
import datetime

import argparse
from scipy.stats import mode


# In[5]:

# pedigree_280 = pd.read_pickle('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/private_output/pedigree_analysis_unrelateds/pedigree_280.pkl') # 2017_07_26_11AM

# unrelated_set = pedigree_280[pedigree_280.In_Unrelated_Set == 'Yes'].WGS_UUID.tolist()


# In[4]:

def safe_div(x, y, alt=0):
    try:
        return x/y
    except:
        return alt


# In[12]:

def calculate_minor_allele_freq(tran, uuids, chrom_dict, sex_dict):
    
    
    nnref_col = []
    ref_col = []
    ref_uuids_col = []
    nnref_uuids_col = []
    alt_allele_freq = []
    ref_allele_freq = []
    minor_allele_freq = []
    minor_allele = []
    missing_col = []
    alleles_dist_col = []

    for c in tran.columns:
        chrom = chrom_dict[c]
        l = tran[c].tolist()
        
        #uuids containing the 1 allele
        non_ref = 0
        non_ref_uuids = []
        
        #uuids containing the 0 allele 
        ref = 0
        ref_uuids = []
        
        ref_allele = 0
        alt_allele = 0
        missing = 0
        # for chroms that are diploid
        # correct for allelic probs on sex chroms
        gts_dict = {'0/0':0, '0/1':1, '1/1':2, './.':0}
        gts_dict_sex = {'0/0':0, '0/1':1, '1/1':1, './.':0}
        alleles_dist = {}
        
        count = 0
        for c in l:
            uuid = uuids[count]
            sex = sex_dict[uuid]
            
            c = copy.deepcopy(str(c))
            alleles_dist[c] = alleles_dist.get(c, 0) + 1
            
            if c == './.':
                missing +=1
                continue
            
            if chrom=='Y':
                if sex == 'M':
                    aa = gts_dict_sex[c]
                    ra = 1 - aa
                    
                    
            
            elif chrom=='X':
                if sex == 'M':
                    aa = gts_dict_sex[c]
                    ra = 1 - aa
                else:
                    aa = gts_dict[c]
                    ra = 2 - aa
                
                
            else:
                try:
                    aa = gts_dict[c]
                except:
                    print c
                    break
                ra = 2 - aa
#                 ref_allele += ra
#                 alt_allele += aa
            
            ref_allele += ra
            alt_allele += aa
                
           
            if c not in ['0/0','./.']:
                non_ref +=1
                non_ref_uuids.append(uuids[count])

            if c == '0/0':
                ref +=1
                ref_uuids.append(uuids[count])

        

                
            count +=1
        
        # fix alleles dist for any missing gt types
        
        allele_types = ['0/0', '0/1', '1/1', './.']
        for a in allele_types:
            alleles_dist[a] = alleles_dist.get(a, 0)
        
        
        alleles_dist_col.append(alleles_dist)
        missing_col.append(missing)
        nnref_col.append(non_ref)
        nnref_uuids_col.append(non_ref_uuids)
             
        ref_col.append(ref)
        ref_uuids_col.append(ref_uuids)
        
        
        
        tot_alleles = ref_allele + alt_allele
        
        
        ref_af = safe_div(ref_allele, tot_alleles, alt='All Missing')
        alt_af = safe_div(alt_allele, tot_alleles, alt='All Missing')
        
#         ref_af = ref_allele/tot_alleles
#         alt_af = alt_allele/tot_alleles
        
        afs = [ref_af, alt_af]
        maf = min(afs)
        
        if afs.index(maf) == 0:
            minor_allele.append('REF')
        else:
            minor_allele.append('ALT')
        
        minor_allele_freq.append(maf)
        ref_allele_freq.append(ref_af)
        alt_allele_freq.append(alt_af)
        
    
    out_names = ['NNREF', 'NNREF_UUIDs','ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'REF_UUIDS', 'Missing', 'Alleles_Dist']
    out_data = [nnref_col, nnref_uuids_col, alt_allele_freq, ref_allele_freq, minor_allele_freq, minor_allele, ref_col, ref_uuids_col, missing_col, alleles_dist_col]
    
    data_dict = dict(zip(out_names, out_data))
    
    return data_dict
    


# In[ ]:

def annotate_gencode_genes(cnv_info):
    
    
    tss_bt = pbt.BedTool('/publicdata/gencode_v19_20151104/tss_merged.bed')                     
    genes = pbt.BedTool('/publicdata/gencode_v19_20151104/genes.bed')
    exons = pbt.BedTool('/publicdata/gencode_v19_20151104/exons.bed')
    
    transcript_to_gene = '/publicdata/gencode_v19_20151104/transcript_to_gene.tsv'
    tg= pd.read_table(transcript_to_gene, index_col=0, header=None, squeeze=True)
    
    # need the chr prefix to intersect gencode 
    
    cnv_info.Chr = cnv_info.Chr.astype(str)
    cnv_info['chrom'] = cnv_info.Chr.apply(lambda x : 'chr' + x)
    
    cnv_info = cnv_info.sort_values(['Chr', 'Start', 'End'])
    
    
    cnv_bt = pbt.BedTool.from_dataframe(cnv_info[['chrom','Start', 'End', 'ID']])
    
    
    cnv_info.index = cnv_info.ID
    
    # Find genes that the CNV overlaps.
    res = cnv_bt.intersect(genes, wo=True)
    df = res.to_dataframe()
    df['gene'] = df.thickEnd
    gb = df[['name', 'gene']].groupby('name')
    se = pd.Series(dict(list(gb['gene'])))
    cnv_info['overlaps_gene'] = se.apply(lambda x: set(x))

    # Find genes that the CNV contains completely.
    df = df[df.blockSizes == df.thickStart - df.strand]
    gb = df[['name', 'gene']].groupby('name')
    se = pd.Series(dict(list(gb['gene'])))
    cnv_info['contains_gene'] = se.apply(lambda x: set(x))  

    # Annotate with genes where the CNV overlaps exonic regions.
    res = cnv_bt.intersect(exons, wo=True)
    df = res.to_dataframe()
    df['gene'] = df.thickEnd.apply(lambda x: tg[x])
    gb = df[['name', 'gene']].groupby('name')
    se = pd.Series(dict(list(gb['gene'])))
    cnv_info['overlaps_gene_exon'] = se.apply(lambda x: set(x))    

    # Distance to nearest TSS.
    res = cnv_bt.closest(tss_bt, D='b')
    df = res.to_dataframe()
    cnv_info.loc[df.name, 'nearest_tss_dist'] = df.thickEnd.values
    
    
    return cnv_info


# In[ ]:

def dist_lambda(x,dict_in):
    Chr = str(x.Chr)
    Start = int(x.Start)
    End = int(x.End)
    
    cent_start = dict_in[Chr][0]
    cent_end = dict_in[Chr][1]
    
    ## CNV before
    if Start < cent_start and End < cent_start:
        dist = End - cent_start
    ## CNV after 
    
    elif Start > cent_end:
        dist = Start- cent_end
    
    ## CNV overlaps right edge:
    elif Start < cent_start and End > cent_start and End < cent_end:
        dist = 0
    ## CNV overlaps left edge
    
    elif Start > cent_start and Start < cent_end and End > cent_end:
        dist=0
    
    ## CNV Overlaps entire region:
    elif Start < cent_start and End > cent_end:
        dist =0
    ## CNV_within centromere entirely
    elif Start > cent_start and Start < cent_end and End > cent_start and End < cent_end:
        dist = 0
    return dist


# In[ ]:

def annotate_filters(cnv_info):
    
    gs_BT = pbt.BedTool.from_dataframe(cnv_info[['Chr','Start', 'End', 'ID']]).sort()
    
    cent_BT = pbt.BedTool('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/output/publicdata/centromeres_merged.bed')
    lcrs_BT = pbt.BedTool('/frazer01/home/djakubosky/software/wham/data/LCR-hs37d5.bed')

    seg_dupes_BT = pbt.BedTool('/frazer01/home/djakubosky/masks/hg19.segdup.mod.bed')  
    
    vdj_BT = pbt.BedTool('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/output/publicdata/vdj.bed')
    mhc_BT = pbt.BedTool('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/output/publicdata/mhc.bed')
    
    pseudo_x_bt = pbt.BedTool('X\t60001\t2699520', from_string=True)
    
    
    centromeres = pd.read_pickle('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/output/publicdata/centromeres.pkl')
    cent_dict = {}
    for x,y,z in zip(centromeres.chrom, centromeres.start, centromeres.end):
        cent_dict[str(x)] = [int(y),int(z)]
        
        
    telomeres = pd.read_pickle('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/output/publicdata/telomeres.pkl')
        
    tel_start = telomeres[telomeres.ID.apply(lambda x: x.split('_')[1]=='Start')]
    tel_end = telomeres[telomeres.ID.apply(lambda x: x.split('_')[1]=='End')]
    tel_start_dict = {}
    for x, y, z in zip(tel_start.Chr, tel_start.Start, tel_start.End):
        tel_start_dict[x] = [int(y), int(z)]
    tel_end_dict = {}
    for x, y, z in zip(tel_end.Chr, tel_end.Start, tel_end.End):
        tel_end_dict[x] = [int(y), int(z)] 
    
    for label,bed in zip(['MHC', 'VDJ', 'centromere', 'seg_dupe', 'pseudo_auto'], [mhc_BT, vdj_BT, cent_BT, seg_dupes_BT, pseudo_x_bt]):
        try:

            df = gs_BT.intersect(bed, wa=True).to_dataframe()
            cnv_info[label]=False
            cnv_info.loc[df['name'].tolist(), label] = True
        except:
            cnv_info[label]=False    
  
    
    
    cnv_info['cent_dist']= cnv_info.apply(lambda x: dist_lambda(x, cent_dict), axis=1)
    cnv_info['tel_start_dist']= cnv_info.apply(lambda x: dist_lambda(x,tel_start_dict), axis=1)
    cnv_info['tel_end_dist']= cnv_info.apply(lambda x: dist_lambda(x,tel_end_dict), axis=1)
    
    return cnv_info


# In[17]:

def annotate_new_cols_to_info(data_out, info_df, suffix = False):
    cols = data_out.keys()
    
    if suffix != 'FALSE':
        cols_out = [ i + '_' + suffix for i in cols]
    else:
        cols_out = cols
    
    for c_name, k in zip(cols_out, cols):
        info_df[c_name] = data_out[k]
    
    return info_df


# In[3]:

def annotate_maf(info, gts, samples, sex_dict, suff):
    
    gt_trans = gts[samples].T
    chrom_dict = info['#CHROM'].to_dict()
    data_out = calculate_minor_allele_freq(gt_trans, samples, chrom_dict, sex_dict)
    info = annotate_new_cols_to_info(data_out, info, suffix= suff)
    return info


# In[2]:

def annotate_maf_sample_subsets(sample_lists, suffix_lists, info, gts, sex_dict):
    
    for samps, suff in zip(sample_lists, suffix_lists):
        info = annotate_maf(info, gts, samps, sex_dict, suff)
    
    return info


# # Command Line Argument handling

# In[ ]:

def add_arguments_to_parser(parser):
    

    parser.add_argument("-info", "--info", dest="info", metavar='<info>', help="lumpy info pickle for all variants (av)", required=True)

    
    parser.add_argument("-gt", "--lumpy_gt", dest="lumpy_gt", metavar='<gt>', help="lumpy gt pickle", required=True)
        
    
    parser.add_argument("-caller", "--caller", dest="caller", metavar='<variant caller>', help="MELT/lumpy", required=True)

    
    parser.add_argument("-gender_map", "--gender_map", dest="gender_map", metavar='<gender_map>', help="gender file UUID Sex, tab delimited", required=True)   
    
    parser.add_argument("-intersect", "--intersections", dest="intersect", action = 'store_true', help="intersect with MHC VDJ Centromeres, Telomeres")
    
    
    parser.add_argument("-s", "--samples", dest="samples", metavar='<fn_samples1,fn_samples2>', help="sets of samples to annotate", required=True)
    
    parser.add_argument("-ss_suff", "--suffix_subsets", dest="suffix_subsets", metavar='<suff1,suff2,suff3>', help="suffixes to use for naming columns of annotated subsets of samples ex: unrelated, related", required=True, default=False)
    
    
    
#     parser.add_argument("-intersect", "--intersections", dest="intersect", metavar='<True/False>', help="intersect with MHC VDJ Centromeres, Telomeres", required=False, default=True)
    
#     parser.add_argument("-somatic", "--somatic", dest="somatic", metavar='<Sample,region,svtype>', help="annotate presence of a known somatic variant, svtype indicates which variant type the somatic variant is.  useful if you know that you have a somatic variant with other evidence such as arrays from different cell types on the same individuals", required=False, default=False)
    
    
    parser.add_argument("-genes", "--genes", dest="genes", action = 'store_true', help="annotate intersection with gencode genes")
        
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    parser.add_argument("-pre", "--prefix", dest="prefix", metavar='<prefix>', help="prefix to name files", default = False)
    
    parser.add_argument("-suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to annotate extracted info and genotypes from lumpy/melt with various things such as gene intersections, minor allele frequency etc.')
    add_arguments_to_parser(parser)
    return parser


# In[ ]:

def run_from_args(args):
    info = pd.read_pickle(args.info)
    gts = pd.read_pickle(args.lumpy_gt)
    caller = args.caller
    output_dir = args.output_dir
    
    # sample sets
    samples = str(args.samples)
    sample_files = samples.split(',')
    
    sample_lists = []
    for fn in sample_files:
        samples = [line.rstrip() for line in open(fn)]
        sample_lists.append(samples)
    
    # suffixes for sample sets in info
    suffixes = args.suffix_subsets
    suffixes = suffixes.split(',')
    
    gender_file = args.gender_map
    
    sex_dict = {line.rstrip().split()[0]:line.rstrip().split()[1] for line in open(gender_file)}
   
    
    print 'starting annotation'
    print CM.datestring(hour=True, minute=True)
    
    # fix a few naming conventions fix irregularities from BND calls
    # dtypes will be screwed up if we don't make some helper columns
    
    info['Chr'] = info['#CHROM'].astype(str)
    info['Start'] = info['POS'].astype(int)
    
    try:
        
        inds = info[info.END=='Column_Not_Present'].index.tolist()
        info.loc[inds, 'END'] = info.loc[inds, 'Start']
    except:
        pass
    
    info['End'] = info['END'].astype(int)
    
    gts = gts.loc[info.index.tolist()]
    info = annotate_maf_sample_subsets(sample_lists, suffixes, info, gts, sex_dict)
    
    
    
    if args.intersect:
        print "annotating centromere/telomere distances, MHC, VDJ Regions"
        info = annotate_filters(info)
#         print gs_info
    
#     if args.somatic:
#         print "annotating somatic variants"
#         spl = args.somatic.split(',')
#         uuid = spl[0]
#         region = spl[1]
#         svtype = spl[2]
#         gs_info = annotate_somatic_var(gs_info, uuid, region, svtype)
        
    if args.genes == True:
        print "annotating gencode genes"
        info = annotate_gencode_genes(info)
        
        
    
    if args.suffix:
        fn_info = os.path.join(output_dir, '{}_info'.format(caller) + args.suffix)
        var_name_info = '{}_info'.format(caller) + args.suffix

    
    else:
        fn_info = os.path.join(output_dir, '{}_info'.format(caller))
        var_name_info = '{}_info'.format(caller)
        
    
    
    print 'data annotated'
    CM.save_dataframe(var_name_info, info, output_dir, print_vars_recorded_loc=False)
#     CM.save_dataframe(var_name_cns, gs_cns, output_location, print_vars_recorded_loc=False)


# In[ ]:

2


# In[ ]:

2


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))

