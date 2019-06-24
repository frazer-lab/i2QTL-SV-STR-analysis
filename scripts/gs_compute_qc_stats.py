
# coding: utf-8

# In[ ]:

from __future__ import division
import numpy as np
import os
import sys
import datetime
from subprocess import call
import subprocess
import glob
import djPyi2 as DJ
from djPyi2 import Common as CM
from djPyi2 import mpltools as axtools

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
dy_name = 'gs_qc_analysis'

private_out = os.path.join(DJ.root, 'private_output', dy_name)
if not os.path.exists(private_out):
    cpy.makedir(private_out)
import gzip


# # Chi Square Test of the Possible batch effects

# In[2]:

def nz_alleles(dist1, dist2):
    nz1 = [z for z in dist1.keys() if dist1[z] != 0]
    nz2 = [z for z in dist2.keys() if dist2[z] != 0]
    
    keys_to_test = list(set(nz1 + nz2))
    keys_to_test = sorted(keys_to_test)
    return keys_to_test

def test_chi(keys, d1, d2):
    
    def add_missing_keys(keys, d):
        d = d.copy()
        for k in keys:
            d[k] = d.get(k, 0)
        return d
    d1 = add_missing_keys(keys,d1)
    d2 = add_missing_keys(keys,d2)
    
    v1 = [d1[z] for z in keys]
    v2 = [d2[z] for z in keys]
    g, p, dof, expctd = stats.chi2_contingency([v1, v2])
    return p, v1, v2

def get_chi_p(d1, d2):
    keys = nz_alleles(d1, d2)
    p,v1, v2 = test_chi(keys, d1, d2)
    return p, keys, v1, v2


def compare_var_sites(comp, suffix1, suffix2):
    names = ['diploid_alleles', 'alleles_dist']
    n1 = ["_".join([i, suffix1]) for i in names]
    n2 = ["_".join([i, suffix2]) for i in names]
    print n1
    print n2
    da1 = comp[n1[0]].tolist()
    da2 = comp[n2[0]].tolist()
    dist1 = comp[n1[1]].tolist()
    dist2 = comp[n2[1]].tolist()
    inds = comp.index.tolist()
    
    data = []
    
    for d1, d2, ad1, ad2 in zip(da1, da2, dist1, dist2):
        s1 = set(d1)
        s2 = set(d2)
        
        same_gt = (ad1 == ad2)
        same_alleles = False
        s1_diff = set()
        s2_diff = set()
        if s1 == s2:
            same_alleles = True
        else:
            s1_diff = s1.difference(s2)
            s2_diff = s2.difference(s1)
            
 
        chi_p, keys, v1, v2 = get_chi_p(ad1, ad2)
      
        corr = np.corrcoef(v1, v2)[0][1]         
        data.append([same_alleles, s1_diff, s2_diff, corr, chi_p, keys, same_gt])
    
    df = pd.DataFrame(data, columns=['same_alleles', 'unique_alleles_1', 'unique_alleles_np', 'corr', 'chi_p', 'keys_compared', 'same_gts'], index=inds)
    
    df['ID'] = comp.index
    df['ad_ipscore'] = comp[n1[1]]
    df['ad_hipsci'] = comp[n2[1]]
    
    chi_p_thresh = 0.05 / df.shape[0]
    df['filter_chi_p'] = df.chi_p < chi_p_thresh
    
    
    return df
    

def prepare_allele_dist_batch(info):
    
    info = info.copy()
    
    data = []
    for d1, x1, d2, x2 in zip(info.alleles_dist_ipscore_unrel.tolist(), 
                              info.num_lq_ipscore_unrel.tolist(),info.alleles_dist_hipsci_fib.tolist(),
                             info.num_lq_hipsci_fib.tolist()):
        d1 = copy.deepcopy(d1)
        d1['LQ'] = x1
        d2 = copy.deepcopy(d2)
        d2['LQ'] = x2
        
        out = [d1, d2]
        data.append(out)
    
    
    df = pd.DataFrame(data, columns= ['alleles_dist_ipscore_unrel', 'alleles_dist_hipsci_fib'], index = info.index)
    
    cols = ['diploid_alleles_ipscore_unrel', 'diploid_alleles_hipsci_fib']
    df = df.join(info[cols])
    return df

def prepare_info_v1_batch(info):
    cols = ['alleles_distipscore_unr', 'alleles_disthipsci_fib', 'diploid_allelesipscore_unr','diploid_alleleshipsci_fib']
    cols_rename = ['alleles_dist_ipscore_unrel', 'alleles_dist_hipsci_fib', 'diploid_alleles_ipscore_unrel','diploid_alleles_hipsci_fib']
    
    df_batch = info_v1[cols].copy()
    df_batch.columns = cols_rename
    return df_batch

def filter_info_all_annot(df):
    df = df[(df.FILTER_GSCNQUAL == False) & (df.somatic == False) & (df.primary_site == True) & (df.stitch_constituent==False) & (df.cnv_class != 'Non_Bi')].copy()
    return df

def calculate_chi_and_annotate_info(info):
    df_batch = prepare_allele_dist_batch(info)
    comp = compare_var_sites(df_batch, 'ipscore_unrel', 'hipsci_fib')
    inds= comp[comp.filter_chi_p ==True].index.tolist()
    info['filter_batch'] = False
    info.loc[inds, 'filter_batch'] = True
    info['batch_p'] = comp.chi_p
    return comp, info
    
    


# # HWE 

# In[19]:

def hardy_weinberg_expectation(dist1):
    keys = ['0/0', '0/1', '1/1']
    vals1 = [dist1.get(i, 0) for i in keys]
    dist = {i:dist1.get(i, 0) for i in keys}
    tot = sum(vals1)
    
    try:
        p = (2*dist['0/0'] + dist['0/1'])/(2*tot)
    except:
        print dist
        return
    q = 1 - p
    
    exp_AA = (p**2) * tot
    exp_Aa = 2*p*q*tot
    exp_aa = (q**2)*tot
    
    hw_exp = [exp_AA, exp_Aa, exp_aa]
    return vals1, hw_exp

def convert_dist_to_alleles(x, cnv_class = 'DUP'):
    
    convert_dict_del = {2: '0/0', 1: '0/1', 0: '1/1'}
    convert_dict_dup = {2: '0/0', 3: '0/1', 4: '1/1'}
    dict_out = {}
    if cnv_class == 'DUP':
        
        for i in x.keys():
            k = convert_dict_dup.get(i, False)
            if not k:
                return False
            else:
                dict_out[k] = x[i]
        return dict_out
                
    if cnv_class == 'DEL':
        for i in x.keys():
            k = convert_dict_del.get(i, False)
            if not k:
                return False
            else:
                dict_out[k] = x[i]
        return dict_out

def prep_info_hwe(info):
    info = info.copy()
    info = info[info.cnv_class != 'Non_Bi']
    info = info[(info.cnv_class != 'mCNV')]
    info = info[~info.Chr.isin(['X', 'Y'])]
    return info

def get_hwe_chi_comparison_pval(dist):
    
    obs, hwe = hardy_weinberg_expectation(dist)
    inds = get_non_zero_inds(obs, hwe)
    obs_test, hwe_test = (subset_list(obs, inds), subset_list(hwe, inds))
    s = stats.chisquare(obs_test, f_exp=hwe_test)
    p_val = s.pvalue
    
    return p_val

def get_non_zero_inds(l1, l2):
    nz1 = [i for i,l in list(enumerate(l1)) if l !=0]
    nz2 = [i for i,l in list(enumerate(l2)) if l !=0]
    inds_to_test = list(set(nz1 + nz2))    
    return inds_to_test

def subset_list(l, inds):
    return [l[i] for i in inds]

def calculate_hwe_and_annotate_info(info_hwe, info):

#     info_hwe['allele_dist_hwe_corrected'] = info_hwe.alleles_dist.apply(lambda x: convert_dist_to_alleles(x,cnv_class))
    
    data = []
    for i, x in info_hwe.iterrows():
    
        dist_all = x['allele_dist_hwe_corrected']
        site = i
        p_val_all = get_hwe_chi_comparison_pval(dist_all)

        data.append(p_val_all)


    inds = info_hwe.index.tolist()
    cols =['hwe_p']
    hwe = pd.DataFrame(data, index=inds, columns=cols)
    hwe_p_thresh = 0.05 / hwe.shape[0]
    hwe['filter_hwe'] = hwe.hwe_p < hwe_p_thresh
    info_hwe['filter_hwe'] = hwe.filter_hwe
    
    inds = hwe[hwe.filter_hwe].index.tolist()
    info['filter_hwe'] = False
    info.loc[inds, 'filter_hwe'] = True
    info['hwe_p'] = hwe.hwe_p
    # these weren't tested because they were mCNV or non-simple DUP
    info['hwe_p'] = info.hwe_p.fillna(False)
    return info_hwe, hwe, info
    


# In[3]:

def run_hwe_annotations(info_v3):
    info_hwe = info_v3.pipe(prep_info_hwe)

    info_hwe['alleles_str'] = info_hwe.diploid_alleles.apply(lambda x: ",".join([str(i) for i in x]))

    info_hwe_del = info_hwe[info_hwe.cnv_class == 'DEL'].copy()
    info_hwe_del['allele_dist_hwe_corrected'] = info_hwe_del.alleles_dist.apply(lambda x: convert_dist_to_alleles(x, 'DEL'))

    info_hwe_dup = info_hwe[info_hwe.cnv_class == 'DUP'].copy()
    info_hwe_dup['allele_dist_hwe_corrected'] = info_hwe_dup.alleles_dist.apply(lambda x: convert_dist_to_alleles(x, 'DUP'))

    info_hwe_dup = info_hwe_dup[info_hwe_dup.allele_dist_hwe_corrected != False]

    info_hwe_testing = pd.concat([info_hwe_dup, info_hwe_del])
    info_hwe_testing, hwe, info_v3 = calculate_hwe_and_annotate_info(info_hwe_testing, info_v3)
    return hwe, info_v3


# In[ ]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-cns", "--cns", dest="cns", metavar='<cns>', help="genome strip cns pickle (Copy Number States)", required=True)
    
    parser.add_argument("-info", "--info", dest="info", metavar='<info>', help="genome strip info pickle", required=True)
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    

    parser.add_argument("-file_suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line tool to compute hwe and batch effects for genome strip variants')
    
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    gs_info = pd.read_pickle(args.info)
    print 'Testing batch and hwe of gs variants'
    print CM.datestring(hour=True, minute=True)
   
    comp, gs_info = calculate_chi_and_annotate_info(gs_info)
    hwe, gs_info = run_hwe_annotations(gs_info)
    
    comp['ID'] = comp.index
    hwe['ID'] = hwe.index
    
    stats = pd.merge(comp, hwe, how = 'outer')
    
    output_location = args.output_dir
    if args.suffix:
        fn_info = os.path.join(output_location, 'gs_info' + args.suffix)      
        var_name_info = 'gs_info' + args.suffix
        var_name_stats = 'gs_qc_stats' + args.suffix
       
       
    else:
        fn_info = os.path.join(output_location, 'gs_info')
        fn_cns = os.path.join(output_location, 'gs_cns')
        var_name_info = 'gs_info'
        var_name_stats = 'gs_qc_stats' 
        
    
    print 'data annotated'
    print CM.datestring(hour=True, minute=True)
    
    CM.save_dataframe(var_name_stats, stats, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_info, gs_info, output_location, print_vars_recorded_loc=False)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

