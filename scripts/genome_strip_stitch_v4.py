
# coding: utf-8

# In[35]:

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
import argparse

import pandas as pd
import csv
import copy 
import pybedtools as pbt
import ciepy
import cardipspy as cpy
import itertools
import tempfile
import six
import networkx as nx
from scipy.stats import mode
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import  make_axes_locatable
import datetime

from scipy.stats import mode
import scipy.stats as stats


# In[36]:

from collections import Counter


# In[37]:

def combine_cnvs(cnvs, info):
    """Combine the list of CNVs cnvs into a single CNV of form CNV_{chrom}_{start}_{end}."""
    return 'CNV_{}_{}_{}'.format(info.loc[cnvs[0], 'Chr'], 
                                 info.loc[cnvs, 'Start'].min(), 
                                 info.loc[cnvs, 'End'].max())

def merge_cnvs(a, b, info, cns):
    """Return boolean indicating whether CNVs a and b should be merged."""
    return ((abs(info.loc[a, 'order'] - info.loc[b, 'order']) == 1) and
            ((cns.loc[a] - cns.loc[b]).abs().mean() < 0.5))


# In[38]:

def compare_lists(l1, l2):
    count = 0
    for i1, i2 in zip(l1, l2):
        if i1 != i2:
            count +=1
    return count
def calculate_absolute_mean_diff(a1, a2):
    """calculate the absolute mean difference of elements in 2 equal length arrays"""
    diff = np.array(a1) - np.array(a2)
    diff_abs = np.absolute(diff)
    mean_diff = diff_abs.mean()
    return mean_diff


# In[76]:

def compare_sites(ind1, ind2, cns_t, samples,samples_nmode1, samples_nmode2, cnmode1, cnmode2, samples_lq1, samples_lq2, subtract_lq = True):
    """ calculate the pearson corr of sites in the intersection, calculate the number of differences among non-mode samples """
    data = []

      
    ids = [ind1, ind2]
    ids_mod = ["_".join(i.split('_')[:-1]) for i in ids]
    
    
    cns_1 = cns_t[ind1].to_dict()
    cns_2 = cns_t[ind2].to_dict()
    
        
    samples_nmode1 = [i for i in samples if cns_1[i] != cnmode1]
    samples_nmode2 = [i for i in samples if cns_2[i] != cnmode2]
    samples_to_compare = set(samples_nmode1 + samples_nmode2)
    if subtract_lq:
        # exclude LQ samps
        samples_to_exclude = set(samples_lq1 + samples_lq2)
        samples_to_compare_nmode = list(samples_to_compare.difference(samples_to_exclude))
        
        samples_to_compare_nmode = list(samples_to_compare.difference(samples_to_exclude))
        samples_to_compare_corr = list(set(samples).difference(samples_to_exclude))
    
    else: 
        samples_to_compare_corr = samples
        samples_nmode1 = [i for i in samples if cns_1[i] != cnmode1]
        samples_nmode2 = [i for i in samples if cns_2[i] != cnmode2]
        samples_to_compare_nmode = list(set(samples_nmode1 + samples_nmode2))
        samples_to_exclude = []
    
    cns_nmode_1 = [int(cns_1[i]) for i in samples_to_compare_nmode]
    cns_nmode_2 = [int(cns_2[i]) for i in samples_to_compare_nmode]

    cns_corr_1 = [int(cns_1[i]) for i in samples_to_compare_corr]
    cns_corr_2 = [int(cns_2[i]) for i in samples_to_compare_corr]
    
    allele_dist1 = dict(Counter(cns_corr_1))
    allele_dist2 = dict(Counter(cns_corr_2))
    
    corr_coef = stats.pearsonr(cns_corr_1, cns_corr_2)[0]
        
    nsamp = len(samples_to_compare_nmode)
    nsamp_pass = len(samples_to_compare_corr)

    num_diff = compare_lists(cns_nmode_1, cns_nmode_2)
    alleles = set(cns_corr_1 + cns_corr_2)
    
    alleles1 = set(cns_corr_1)
    num_alleles1 = len(alleles1)
    alleles2 = set(cns_corr_2)
    num_alleles2 = len(alleles2)

    mean_diff_all = calculate_absolute_mean_diff(cns_corr_1, cns_corr_2)
    mean_diff_nmode = calculate_absolute_mean_diff(cns_nmode_1, cns_nmode_2)
    
    exact_match = (cns_corr_1 == cns_corr_2)

    
    try:
        perc_diff = num_diff/nsamp
    except:
#         print nsamp, 'nsamp is zero'
        perc_diff = 0

    out = [ind1, ind2, corr_coef, num_diff, nsamp, 
           perc_diff, samples_to_compare_nmode, samples_to_compare_corr, 
           list(samples_to_exclude), nsamp_pass, exact_match, mean_diff_all, mean_diff_nmode, 
           alleles1, alleles2, num_alleles1, num_alleles2, allele_dist1, allele_dist2]

    return out


# In[40]:

def split_id_to_coord(ID):
    spl = ID.split('_')
    chrom = spl[1]
    start = int(spl[2])
    end = int(spl[3])
    return chrom, start, end


# In[41]:

def compute_dist(id1, id2):
    chrom1, start1, end1 = split_id_to_coord(id1)
    chrom2, start2, end2 = split_id_to_coord(id2)
    dist = start2 - end1
    return dist
    


# In[65]:

def collect_data_adjacent_sites(info, cns_t, samples, cn_mode_col = 'cn_mode', subtract_lq=True):
    
    data = []

    for chrom, df in info.groupby('Chr'):
        inds = df.index.tolist()
        diff_mode_uuids = info.diff_mode_uuids.to_dict()
        mode_cn_all = info[cn_mode_col].to_dict()
        lq_uuids = info.lq_samps.to_dict()
        cnv_classes = info.cnv_class.to_dict()


#         if chrom not in ['X', 'Y']:
        max_range = len(inds)-1 
        for i in range(0, max_range):
            if i < max_range:

                ind1 = inds[i]
                ind2 = inds[i+1]
                pair_ind = '-'.join([ind1, ind2])

                distance_between = compute_dist(ind1, ind2)
                if distance_between < 0:
                    absolute_dist = 0
                else:
                    absolute_dist = distance_between

                cnv_class1 = cnv_classes[ind1]
                cnv_class2 = cnv_classes[ind2]

                diff_uuids1 = diff_mode_uuids[ind1]
                diff_uuids2 = diff_mode_uuids[ind2]
                mode_cn1 = mode_cn_all[ind1]
                mode_cn2 = mode_cn_all[ind2]
                lq_uuids1 = lq_uuids[ind1]
                lq_uuids2 = lq_uuids[ind2]
                pair = [ind1,ind2]

                num_diff1 = len(diff_uuids1)
                num_diff2 = len(diff_uuids2)

                comp = compare_sites(ind1, ind2, cns_t, samples, 
                                     diff_uuids1, diff_uuids2, mode_cn1, mode_cn2, 
                                     lq_uuids1, lq_uuids2, subtract_lq=subtract_lq)

                comp = comp + [pair, mode_cn1, mode_cn2, cnv_class1, cnv_class2, distance_between, 
                               absolute_dist, pair_ind, num_diff1, num_diff2, diff_uuids1, diff_uuids2,
                              chrom]
                data.append(comp)

    df = pd.DataFrame(data, columns=['ID1', 'ID2', 'corr_coef',
                                     'num_diff', 'num_non_mode', 'percent_non_mode_diff', 
                                     'samps_to_compare_nmode', 'samps_to_compare_corr',(
                                     'samps_to_exclude', 'num_pass', 'exact_cn_match', 
                                     'mean_cn_diff_all', 'mean_cn_diff_nmode','alleles1', 'alleles2', 
                                     'num_alleles1', 'num_alleles2','allele_dist1', 'allele_dist2','pair', 
                                     'mode_cn1', 'mode_cn2', 'cnv_class1', 'cnv_class2', 
                                     'distance_between', 'distance_between_mod', 'cat_pair', 
                                     'num_diff1', 'num_diff2', 'diff_uuids1', 'diff_uuids2', 'chrom'])
    
    df.index = df.cat_pair
    return df


# In[31]:

def prep_info(df):
    df = df.copy()
    df = df.sort_values(['Chr', 'Start', "End"])
    return df


# In[32]:

def prep_cns(df, info):
    df = df.copy()
    df = df.reindex(info.index)
    return df 


# In[33]:

def annotate_passing(df, thresh):
    df = df.copy()
    inds = df[(df.corr_coef > 0.9) & (df.percent_non_mode_diff <= 0.2) & (df.mean_cn_diff_all < 0.5) &
              (df.mean_cn_diff_nmode < 0.5) & (df.distance_between_mod < thresh)].index.tolist()
    

    new_size = df.shape[0]
    
    
    df['passing_criteria'] = False
    df.loc[inds, 'passing_criteria'] = True
    return df


# In[51]:

def prep_adj_sites(adj_sites_prestich, thresh = 5000):  
    adj_sites_prestich = adj_sites_prestich.copy()
    adj_sites_prestich['matching_mode'] = (adj_sites_prestich.mode_cn1 == adj_sites_prestich.mode_cn2)
#     adj_sites_plot = adj_sites_prestich[(adj_sites_prestich.matching_mode == True) & 
#                                     (adj_sites_prestich.num_alleles1 > 1) & 
#                                         (adj_sites_prestich.num_alleles2 > 1)].copy()
    adj_sites_prestich['log_dist'] = np.log10(adj_sites_prestich.distance_between_mod + 1)
    adj_sites_prestich = adj_sites_prestich.pipe(annotate_passing, thresh)
    return adj_sites_prestich


# In[43]:

info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_combine_ipscore_hipsci/info_all_sites_rmdup_filt.pkl').pipe(prep_info)
cns = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/i2QTL_combined/cns_all.pkl').pipe(prep_cns, info)
cns.drop('old_index', axis =1, inplace=True)


# In[46]:

sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')

samples_discovery = sample_info[sample_info.CELL_TYPE != 'iPSC'].WGS_ID.tolist()

cns_t = cns.T.copy()
cns_t = cns_t.loc[samples_discovery]


# In[117]:

samples_ipscore_males = sample_info[(sample_info.STUDY == 'iPSCORE') & (sample_info.SEX == 'M')].WGS_ID.tolist()


# In[116]:

samples_females = sample_info[(sample_info.CELL_TYPE != 'iPSC') & (sample_info.SEX=='F')].WGS_ID.tolist()


# In[115]:

samples_males = sample_info[(sample_info.CELL_TYPE != 'iPSC') & (sample_info.SEX=='M')].WGS_ID.tolist()


# In[118]:

len(samples_males)


# In[119]:

len(samples_females)


# In[127]:

adj_sites = collect_data_adjacent_sites(info, cns_t, samples_discovery, subtract_lq=True).pipe(prep_adj_sites, thresh=30000)

passing_lq_correct = adj_sites[adj_sites.passing_criteria == True].index.tolist()


# In[128]:

adj_sites_all_males = collect_data_adjacent_sites(info, cns_t, samples_males, cn_mode_col='cn_mode_male', subtract_lq=True).pipe(prep_adj_sites, thresh=30000)


# In[129]:

adj_sites_all_females = collect_data_adjacent_sites(info, cns_t, samples_females, cn_mode_col='cn_mode_female', subtract_lq=True).pipe(prep_adj_sites, thresh=30000)


# In[130]:

adj_sites_male_ipscore = collect_data_adjacent_sites(info, cns_t, samples_males, cn_mode_col = 'cn_mode_male_ipscore_fb', subtract_lq=True).pipe(prep_adj_sites, thresh=30000)

# passing_lq_correct_male_ipscore = adj_sites[adj_sites.passing_criteria == True].index.tolist()


# In[144]:

comp_males_x = adj_sites_all_males[(adj_sites_all_males.chrom == 'X')  & (adj_sites_all_males.passing_criteria == True)].index.tolist() b


# In[175]:

def get_corr_at_samples(ind1, ind2, samples, cns_t):
    cns_1 = cns_t[ind1].to_dict()
    cns_2 = cns_t[ind2].to_dict()
    try:

        cns_corr1 = [int(cns_1[i]) for i in samples]
        cns_corr2 = [int(cns_2[i]) for i in samples]
        corr_coef = stats.pearsonr(cns_corr1, cns_corr2)[0]
    except:
        print samples
    return corr_coef


# In[176]:

def gather_consensus_X(adj_sites_all_males, adj_sites_all_females, adj_sites, cns_t):
    inds = adj_sites_all_males.index.tolist()
    assert (inds == adj_sites_all_females.index.tolist())
    
    num_non_mode_males = adj_sites_all_males.num_non_mode.tolist()
    num_non_mode_females = adj_sites_all_females.num_non_mode.tolist()
    
    ndiff_male = adj_sites_all_males.num_diff.tolist()
    ndiff_female = adj_sites_all_females.num_diff.tolist()
    
    passing_b_males = adj_sites_all_males.passing_criteria.tolist()
    passing_b_females = adj_sites_all_females.passing_criteria.tolist()
    
    samples_nm_males = adj_sites_all_males.samps_to_compare_corr.tolist()
    samples_nm_females = adj_sites_all_females.samps_to_compare_corr.tolist()
    union_samples = [list(set(i1 + i2)) for i1, i2 in zip(samples_nm_males, samples_nm_females)]
   
    return union_samples
#     corr_all = adj_sites.corr_coef.tolist()
    
    corr_males = adj_sites_all_males.corr_coef.tolist()
    corr_females = adj_sites_all_females.corr_coef.tolist()
    dist = adj_sites_all_males.distance_between_mod.tolist()
    
    out = []
    for ind, nm_male, nm_female, pass_b_male, pass_b_female, cm, cf, ndm, ndf, d, union_samps in zip(inds,
                                                                                            num_non_mode_males,
                                                                           num_non_mode_females,
                                                                           passing_b_males, passing_b_females, 
                                                                           corr_males, corr_females,
                                                                           ndiff_male, ndiff_female, 
                                                                           dist, union_samples):
        
        
        ind1, ind2 = ind.split('-')
        chrom = ind.split('_')[1]
        
        if chrom == 'X':

            if [nm_male, nm_female] == [0, 0]:
                out.append([ind, False])

            elif (nm_male == 0) & (nm_female > 0):
                out.append([ind, pass_b_female])

            elif (nm_male > 0) & (nm_female == 0):
                out.append([ind, pass_b_male])

            elif (nm_male > 0) & (nm_female > 0):
                if all([pass_b_male, pass_b_female]):
                    out.append([ind, True])
                else:
                    perc_diff = (nm_male + nm_female)/(ndm + ndf)
                    corr = get_corr_at_samples(ind1, ind2, union_samples, cns_t)
                    if all([(cm > 0.9), (cf > 0.9), (perc_diff <=0.2), (d < 30000), (corr > 0.9)]):
                        print 'yes'
                        out.append([ind, True])
                    else:
                        out.append([ind,False])

            else:
                print "didn't account for all scenarios"
                break
    return out    


# In[177]:

x_chrom_consensus = gather_consensus_X(adj_sites_all_males, adj_sites_all_females, adj_sites, cns_t)


# In[178]:

x_chrom_consensus


# In[145]:

comp_females_x = adj_sites_all_females[(adj_sites_all_females.chrom == 'X')  & (adj_sites_all_females.passing_criteria == True)].index.tolist()


# In[146]:

x_chrom_consensus = set(comp_males_x).intersection(set(comp_females_x))


# In[147]:

len(x_chrom_consensus)


# In[148]:

len(comp_females_x)


# In[149]:

len(comp_males_x)


# In[83]:

len(samples_ipscore_males)


# In[73]:

adj_sites_male_ipscore[adj_sites_male_ipscore.num_non_mode == 0]


# In[64]:

adj_sites[~adj_sites.chrom.isin(['X', 'Y'])]


# In[ ]:




# In[ ]:




# In[54]:

adj_sites.


# In[3]:

def combine_connected_component_length(cnvs, chrom, info, cns, cn_modes, cn_modes_males, cn_modes_females, all_samples, males, females, order_dict, length=5000, corr_nmode = 0.8):
    
    """
    If appropriate, combine CNVs from the list cnvs. The CNVs in cnvs should
    have highly correlated genotypes (for instance, cnvs might be a list of connected
    nodes from a graph built from the genotype correlations).
    
    Parameters
    ----------
    cnvs : list
        List of CNVs of the form CNV_2_400_800 (aka CNV_{chrom}_{start}_{end}).
        
    info : pandas.DataFrame
        Data frame whose index contains the values in cnvs. This dataframe should at least
        have chrom, start, and end columns.
        
    cns : pandas.DataFrame
        Data frame whose index contains the values in cnvs. This data frame should contain
        the diploid copy number estimates.
        
    cn_modes: dict of modes to use for a given variant
    
    males: samples who are male- for y chrom variants this is important
        
    length: length threshold between merged cnvs (use this to stop CNV merging for things that are very far apart)
    
    """
    cnv_tup = [(i, order_dict[i]) for i in cnvs]
    cnv_tup = sorted(cnv_tup, key=lambda x: x[1])
    cnvs = [i[0] for i in cnv_tup]
    out = dict()
    
    def check_pairs(cnvs, cn_modes, samples, sample_set):
        """ sample set refers to set of samples (males females, all)"""
        combine = []
        data = []
        for i in range(len(cnvs) - 1):

            cnv1 = cnvs[i]
            cnv2 = cnvs[i+1]
            cn_mode1 = cn_modes[cnv1]
            cn_mode2 = cn_modes[cnv2]

            length_check = (info.loc[cnv2, 'Start'] - info.loc[cnv1, 'End']) < length


            order_check = (abs(info.loc[cnv1, 'order'] - info.loc[cnv2, 'order']) == 1)

            mean_gt_diff_check = ((cns.loc[cnv1] - cns.loc[cnv2]).abs().mean() < 0.5)


            corr_nmode_bool, inf, cr, nsamp = check_corr_non_mode_CN(cnv1, cnv2, cns, cn_mode1, cn_mode2, samples, corr_thresh = corr_nmode)

#             data.append([cnv1, cnv2, corr_nmode_bool, inf, cr, nsamp, sample_set])
            status = True
            test_metrics = [length_check, order_check, mean_gt_diff_check, corr_nmode_bool]
            try:
                i = test_metrics.index(False)
                status =False
            except:
                pass
              
            data.append([cnv1, cnv2, corr_nmode_bool, inf, cr, nsamp, sample_set, status, 
                         length_check, order_check, mean_gt_diff_check, cnvs])
            combine.append(status)
            
        return combine, data
    
    
   
    
    if chrom ==  'X':
        combine_females, data_females = check_pairs(cnvs, cn_modes_females, females, 'females')
        combine_males, data_males = check_pairs(cnvs, cn_modes_males, males, 'males')
        combine = []
        data = []
        for b1, b2, d1, d2 in zip(combine_females, combine_males, data_females, data_males):
            nsamp_males = d2[5]
            nsamp_females = d1[5]
            pass_females = b1
            pass_males = b2
            
            if ((nsamp_males > 0) & (nsamp_females > 0)):          
                b_combined = ([pass_females, pass_males] == [True, True])
                combine.append(b_combined)
                # collect data from both
                data.append(d1)
                data.append(d2)
            
            elif (nsamp_females == 0) & (nsamp_males > 0):
                combine.append(pass_males)
                # collect from just males
                data.append(d2)
                
            elif (nsamp_females == 0) & (nsamp_males == 0):
                combine.append(False)
                data.append(d2)
            
            elif (nsamp_females > 0) & (nsamp_males == 0):
                combine.append(pass_females)
                data.append(d1)
            else:
                print cnvs, nsamp_males, nsamp_females, pass_males, pass_females, ((nsamp_males > 0) and (nsamp_females > 0)), b1, b2,  ([pass_females, pass_males] == [True, True])
                break
    
    elif chrom ==  'Y':
        combine, data = check_pairs(cnvs, cn_modes_males, males, 'males')
    
    else:
        combine, data = check_pairs(cnvs, cn_modes, all_samples, 'all_samples')
            
        
    i = 0    
    to_combine = [cnvs[0]]
    while i < len(cnvs) - 1:
        if combine[i]:
            to_combine.append(cnvs[i + 1])
        else:
            if len(to_combine) > 1:
                out[combine_cnvs(to_combine, info)] = to_combine
            to_combine = [cnvs[i + 1]]
        i += 1
        
    if len(to_combine) > 1:
        out[combine_cnvs(to_combine, info)] = to_combine

    other_info = pd.DataFrame(data, columns = ['cnv1', 'cnv2', 'stitch_b', 'inf', 'corr', 'nsamp', 'sample_set', 'status_stitch', 'length_b', 'order_b', 'mean_gt_b', 'original_cnv_clust'])
    
    
    other_info['m_clust'] = str(out)
    
    
    return out, other_info


# In[4]:

def characterize_cnv_classes(x):
    Types = []
    out = 'none'
    if ('DUP' in x) and ('mCNV' in x) and ('DEL' not in x):
        out = 'DUP,mCNV'
    if ('DUP' in x )and ('DEL' not in x) and ('mCNV' not in x):
        out = 'DUP,DUP'
    if ('DEL' in x) and ('DUP' not in x) and ('mCNV' not in x):
        out = 'DEL,DEL'
    if ('DEL' in x) and ('mCNV' in x) and ('DUP' not in x):
        out = 'DEL,mCNV'
    if ('DEL' in x) and ('mCNV' in x) and ('DUP' in x):
        out = 'DEL,DUP,mCNV'
    if ('DEL' in x) and ('DUP' in x) and ('mCNV' not in x):
        out = 'DUP,DEL'
    if ('mCNV' in x) and ('DUP' not in x) and ('DEL' not in x):
        out = 'mCNV,mCNV'
        
    return out


# In[2]:

def lambda_add_stitched_tag(x):
    if x.split('_')[-1] in ['iPSCORE', 'HipSci']:
        return x
    else:
        return x + '_Stitched'


# In[76]:

def stitch_cnvs(gs_info, gs_cns, sex_dict, samples, max_distance=5000, correl = 0.9, corr_for_nmode = 0.8,
               males_Y=False, females =False, cn_mode_col = 'cn_mode', cn_mode_male_Y_col = 'cn_mode_male', cn_mode_male_col = 'cn_mode_male'):
    
    
    
    all_males = [i for i in samples  if sex_dict[i]=='M']
    females = [i for i in samples if sex_dict[i]=='F']
    
    

    gs_info.Chr =gs_info.Chr.astype(str)
    gs_info.End = gs_info.End.astype(int)
    gs_info.Start = gs_info.Start.astype(int)
    gs_info.sort_values(by=['Chr', 'Start', 'End'], inplace=True)
    
   
        
    gs_cns=gs_cns[samples]
    
    
    gs_info['order'] = range(gs_info.shape[0])
    order_dict = gs_info['order'].to_dict()
    
    cn_modes = gs_info[cn_mode_col].to_dict()
    cn_modes_males = gs_info[cn_mode_male_col].to_dict()
    if cn_mode_col != cn_mode_male_Y_col:
        cn_modes_males_Y = gs_info[cn_mode_male_Y_col].to_dict()
    else:
        cn_modes_males_Y = cn_modes_males
        
    cn_modes_females = gs_info.cn_mode_female.to_dict()
    
    
    
    st = pd.DataFrame()
    combined= dict()
    for chrom in set(gs_info.Chr):
        if chrom == 'Y':
            # use the cnm for males only
            t = gs_info[gs_info.Chr == chrom]
            
            # use the CNM of the Y chrom male subset, if there is one (same by default as autosomes)
            cnm_males = cn_modes_males_Y
            
            if males_Y:
                # use only a subset of males for Y chrom calls
                males = males_Y
                corr = gs_cns[males].loc[t.index].T.corr()
            else:
                males = all_males
                corr = gs_cns[males].loc[t.index].T.corr()
                # if no alternative subset provided


        else:
            cnm_males = copy.deepcopy(cn_modes_males)
            # use all males for autosomes (X chrom calls consider both separately)
            males = all_males
            t = gs_info[gs_info.Chr == chrom]
            corr = gs_cns.loc[t.index].T.corr()
                
            

        g = corr > correl
        edges = []
        for i in g.columns:
            se = g[i]
            # select correlated values
            se = se[se]
            for j in se.index:
                edges.append((i, j))
        g = nx.Graph(edges)
        cc = nx.connected_components(g)
        

        while True:
            try:
                c = list(cc.next())
                if len(c) > 1:
                    d, odf = combine_connected_component_length(c, chrom, gs_info, gs_cns, cn_modes, cnm_males, 
                                                                cn_modes_females, samples, males, females, order_dict,
                                                                length=max_distance, corr_nmode = corr_for_nmode)
                    st = st.append(odf)
                        
                    for k in d.keys():
                        combined[k] = d[k]
            except StopIteration:
                break

                
                
    mapping_from = []
    mapping_to = []
    for k in combined.keys():
        for x in combined[k]:
            mapping_from.append(x)
            mapping_to.append(k)
    mapping = pd.Series(mapping_to, index=mapping_from)
#     mapping.to_csv(os.path.join(private_out, suffix + '_combined_mapping.tsv'), sep='\t')
    
    
    to_remove = []
    for k in combined.keys():
        to_remove += combined[k]
    print('{} CNVs combined into {} CNVs.'.format(len(to_remove), len(combined)))
    
        
    data = []
    for i in combined.keys():
        spl = i.split('_')
        chrom = spl[1]
        start,end = int(spl[2]), int(spl[3])
        merged_calls = combined[i]
        len_merge = len(merged_calls)
        classes_merged = []
        cluster_str = ",".join(merged_calls)
        
        for z in merged_calls:
            class_call = gs_info.loc[z].cnv_class
            classes_merged.append(class_call)
        data.append([chrom, start, end, i, len_merge, merged_calls, classes_merged, cluster_str])

    merged = pd.DataFrame(data, columns=['chrom', 'start', 'end', 'stitched_cnv_site_ID','num_cnvs_merged', 'cnvs_merged', 'cnv_classes', 'cluster_str'])
    
    merged.sort_values(['chrom', 'start', 'end'], inplace=True)
    
    merged['merge_class']= merged.cnv_classes.apply(lambda x: characterize_cnv_classes(x))
    merged['Length']= merged.end - merged.start
    merged.index = merged.stitched_cnv_site_ID
    
    distances = []
    for i in combined.keys():
        list_sites = combined[i]
        num_combined = len(list_sites)

        for l,z in enumerate(list_sites[:-1]):
            p1, p2 = list_sites[l],list_sites[l+1]
            end_1 = int(p1.split('_')[3])
            start_2 = int(p2.split('_')[2])

            pair = p1 + ',' + p2

            dist = start_2 - end_1
            CNV_Types = str(gs_info.loc[p1].cnv_class) + ',' + str(gs_info.loc[p2].cnv_class)

            distances.append([pair, i, CNV_Types, dist, num_combined])
            

    dist_frame = pd.DataFrame(distances, columns=['pair', 'stitched_cnv_site_ID' ,'cnv_types', 'distance', 'num_combined'])
    try:
        
        mean_dist = dist_frame.groupby('stitched_cnv_site_ID').distance.mean().to_frame()
        mean_dist['mean_distance_between']= mean_dist.distance
        merged = merged.join(mean_dist['mean_distance_between'])
        
    except:
        pass
    

    cluster_ID_dict = {}
    stitched_ID_dict = {}
    count = 1
    for x1, x2 in zip(merged.cnvs_merged.tolist(), merged['stitched_cnv_site_ID'].tolist()):
        for ID in x1:
            cluster_ID_dict[ID] = count
            stitched_ID_dict[ID] = x2

        cluster_ID_dict[x2] = count 
        count +=1

    merged['cluster_ID'] = merged.stitched_cnv_site_ID.apply(lambda x: cluster_ID_dict[x])
    
    
    tdf = pd.DataFrame([x.split('_') for x in combined.keys()], columns=['cnv', 'Chr', 'Start', 'End'],
                   index=combined.keys()).drop('cnv', axis=1)
   
    tdf['ID']= tdf.index
    tdf = tdf.join(merged[['cluster_ID', 'cluster_str', 'stitched_cnv_site_ID']])
    
    
    tdf.Start = tdf.Start.astype(int)
    tdf.End = tdf.End.astype(int)
    # site generated by stitching?
    tdf['stitch_breakpoint']=True
    # is this a site that is a constituent of a stitching cluster?
    tdf['stitch_constituent'] = False
    
    # what cluster of stitched site does this correspond to- stitched site ID also gets this numeric ID
    
    gs_info['cluster_ID'] = gs_info.ID.apply(lambda x: cluster_ID_dict.get(x, 0))
    # what stitch site is this site a constituent of, if any?
    gs_info['stitched_cnv_site_ID'] = gs_info.ID.apply(lambda x: stitched_ID_dict.get(x, False))
    
    # all the gs_info sites are before stitching- aren't the newly generated stitch breakpoint
    gs_info['stitch_breakpoint'] = False
    
    gs_info['stitch_constituent']= False
    gs_info.loc[to_remove, 'stitch_constituent'] = True
    
    # carry over original info col from the original sites that are unstitched
    gs_info_trunc = gs_info[['Chr','Start', 'End', 'ID','stitch_breakpoint', 
                             'stitch_constituent', 'cluster_ID', 'stitched_cnv_site_ID']]

    gs_combined_info_unannotated = pd.concat([tdf, gs_info_trunc])
    gs_combined_info_unannotated = gs_combined_info_unannotated.sort_values(by=['Chr', 'Start', 'End'])
    

    
    gs_combined_info_unannotated.ID =  gs_combined_info_unannotated.ID.apply(lambda x: lambda_add_stitched_tag(x))
    gs_combined_info_unannotated.index = gs_combined_info_unannotated.ID

#     gs_combined_info_unannotated = pd.concat([tdf, gs_info_trunc.drop(to_remove)])
    
    return merged, dist_frame, gs_combined_info_unannotated, to_remove, combined, st


# In[6]:

def check_corr_non_mode_CN(cnv1, cnv2, gs_cns, cn_mode1, cn_mode2, samples, corr_thresh=0.8):
    
    """additional validation at the level of the non-mode samples- if the samples have less than 4 different alleles, we require less than 10% of samples to be different, by a maximum of 1 CN in any given sample, with a mean difference of less than 0.5 on average, if there are 4 or more alleles also require correlation of 80% at the sample level, this is to ensure that the initial correlation isn't driven only by correlation between mode copy number in adjacent CNVs"""
    
    t1 = (gs_cns.loc[cnv1] != cn_mode1)
    t2 = (gs_cns.loc[cnv2] != cn_mode2)
   
    s1 = t1[t1==True].index.tolist()
    s2 = t2[t2==True].index.tolist()
  
    samples_to_compare = list(set(s1 + s2))
    samples_to_compare = [i for i in samples_to_compare if i in samples]
    
    nsamp = len(samples_to_compare)
  
    chrom = cnv1.split('_')[1]

    
    # no non mode samps
    if nsamp == 0:
        return False, 'none', 'none', nsampdif
    
 
    elif nsamp <= 2:
        
        t = (gs_cns.loc[cnv1, samples_to_compare] == gs_cns.loc[cnv2, samples_to_compare]).value_counts().to_dict()
        b = t.get(True, 0) == nsamp
        return b, 'none', 'none', nsamp
      
        
    else:
        g1 = gs_cns.loc[cnv1, samples_to_compare]
        g2 = gs_cns.loc[cnv2, samples_to_compare]
        t = (g1 == g2).value_counts().to_dict()
        
        # check if all NM are the same CN, if so return true
        b = t.get(True) == nsamp
        
        if b:
            return b, 'none', 'none', nsamp
        
        else:
            # if all are not the same, check number that are different
            # max difference between copy numbers
            # pecentage of samples different
            # mean CN difference
            # if more than 3 distinct nonmode alleles, also check correlation

            d = (g1 == g2)
            d = d.value_counts()

            num_diff = d.get(False, 0)
            diff = np.abs((g1 - g2))

            mdiff = max(diff)
            msum = sum(diff)
            rdiff = num_diff/nsamp
            alleles = set(g1.tolist() + g2.tolist())
            rdiff_b = (rdiff <= 0.2)
            mean_gt_diff_b = (diff.mean() < 0.5)

            diff_b = (mdiff < 2)
            pass_thresh = (rdiff_b, mean_gt_diff_b) == (True, True)
            return pass_thresh,  [diff_b, mean_gt_diff_b, diff_b, num_diff], 'none', nsamp




# In[8]:

def find_header_end(fn):
    """ find header end line number of vcf or gzipped vcf"""
    try:

        if fn.split('.').pop() == 'gz':
            import gzip
            F = gzip.open(fn, 'rU')
        else:
            F = open(fn, 'rU')
    except:
        return "vcf file doesn't have the right extension or isn't a VCF"


    count = 0
    for line in F:
        count +=1
        try: 
            spl = line.split('\t')
            spl0 = spl[0] 
            if spl[0]=="#CHROM":
                F.close()
                return count
            if count > 2000:
                F.close()
                return 'incomplete or missing header, or very long header'
        except:
            break

    F.close()


# In[71]:

def vcf_format_df(df, vcf_df):
    """ takes as input a df in the format of 'merged' df here from cnv 
    merging in order to produce a vcf formatted pandas df for use 
    in constructing a vcf file for re-genotyping"""
    
    data = []
    for l, i in df.iterrows():
        chrom = i.chrom
        POS = i.start
        END = i.end
        ID  = i.stitched_cnv_site_ID
        first_Merged = i.cnvs_merged[0]
        ALT = '<CNV>'
        QUAL='.'
        FILTER='PASS'

        
        SVTYPE='CNV'
            
        INFO = "END={};SVTYPE={}".format(END, SVTYPE)

        REF = vcf_df.loc[first_Merged].REF
        out = [chrom,POS, ID, REF,ALT, QUAL, FILTER, INFO, END]
        data.append(out)

    Header = ['#CHROM', 'POS', 'ID', 'REF', 
                 'ALT', 'QUAL', 'FILTER', 'INFO', 'END']

    df_out = pd.DataFrame(data, columns=Header)
    return df_out


# In[72]:

def split_seq(seq, num_pieces):
    start = 0
    for i in xrange(num_pieces):
        stop = start + len(seq[i::num_pieces])
        yield seq[start:stop]
        start = stop


# In[84]:

def write_genotyping_vcfs(vcf_formatted_df, run_directory,
                         header_location = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_pipeline/header.txt',
                         genotyping_job = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/notebooks/scripts/regenotype.sh'):
    
    # adapted for i2QTL
    
    # Make initial directories for text files etc
    
    text_file_dir = os.path.join(run_directory, 'text_files')
    in_vcf_dir = os.path.join(run_directory, 'in_vcfs')
    logs_dir = os.path.join(run_directory, 'logs')
   
    results_dir = os.path.join(run_directory, 'results')
    dirs = [run_directory, text_file_dir, in_vcf_dir, logs_dir, results_dir]
    
    for d in dirs:
        if not os.path.exists(d):
            os.mkdir(d)
#         call('mkdir ' + d, shell=True)
    
    jobs_all = []
    for i, df in vcf_formatted_df.groupby('#CHROM'):
    
        name = os.path.join(text_file_dir, i + '_merges.txt')
        df.sort_values(by=['#CHROM', 'POS', 'END'], inplace=True)
        df[['#CHROM', 'POS', 'ID', 'REF', 
            'ALT', 'QUAL', 'FILTER', 'INFO']].to_csv(name,
                                                     sep='\t', index=False)

        in_vcf_fn = os.path.join(in_vcf_dir, i + '_merges.vcf')

        #! cat {header_location} {name} > {vcf_name}
        
        #Make the in vcfs
        
        command_in = 'cat {} {} > {}'.format(header_location, name, in_vcf_fn)
        call(command_in, shell=True)


# In[69]:

def read_vcf_to_df(fn):
    skiprows = find_header_end(fn) - 1
    df = pd.read_table(fn, skiprows=skiprows)
    return df

def combine_original_vcfs(vcf1, vcf2):
    def prep_vcf(df, source):
        df = df.copy()
        df['ID'] = df['ID'].apply(lambda x: "{}_{}".format(x, source))
        df.index = df.ID
        return df
    
    df1 = read_vcf_to_df(vcf1).pipe(prep_vcf, 'iPSCORE')
    df2 = read_vcf_to_df(vcf2).pipe(prep_vcf, 'HipSci')
    combined_df = pd.concat([df1, df2])
    return combined_df


# In[12]:

# sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')

# info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_combine_ipscore_hipsci/info_all_sites_rmdup_filt.pkl')

# cns = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_combine_ipscore_hipsci/cns_all_filt.pkl')

# samples_discovery = sample_info[sample_info.CELL_TYPE != 'iPSC'].WGS_ID.tolist()

# vcf_ipscore = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/run_ipscore/results/gs_cnv.genotypes.vcf.gz'
# vcf_hipsci = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/run_hipsci/results/gs_cnv.genotypes.vcf.gz'

# run_dir = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/stitching'

# sex_dict = sample_info.SEX.to_dict()
# males = [k for k,v in sex_dict.iteritems() if v == 'M']
# females = [k for k,v in sex_dict.iteritems() if v == 'F']

# merged, dist_df, gs_combined_info_unannotated, to_remove, combined, st = stitch_cnvs(info, cns, sex_dict, samples_discovery, max_distance=5000, correl=0.9, corr_for_nmode = 0.7)

# merged, dist_df, gs_combined_info_unannotated, to_remove, combined, st = stitch_cnvs(info, cns, sex_dict, samples_discovery, max_distance=5000, correl=0.9, corr_for_nmode = 0.7)

# gs_combined_info_unannotated.shape

# tdf = combine_original_vcfs(vcf_ipscore, vcf_hipsci)

# vcf_formatted_df = vcf_format_df(merged, tdf)
# # make_genotyping_jobs2(merged_vcf, run_dir )

# make_genotyping_jobs2(vcf_formatted_df, run_dir,
#                          header_location = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_pipeline/header.txt',
#                          genotyping_job = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/notebooks/scripts/regenotype.sh')
    


# In[216]:

# fn = '/frazer01/home/djakubosky/BF_GS_Discovery3/results/gs_cnv.genotypes.vcf'
# merged_vcf = vcf_format_df(merged, fn)


# In[13]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to identify neighboring sites from genome strip calls to stitch together based on correlation and similarity of copy number and create jobs to run genotyping')
    add_arguments_to_parser(parser)
    return parser


# In[80]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-cns", "--cns", dest="cns", metavar='<cns>', help="genome strip cns pickle", required=True)
    parser.add_argument("-info", "--info", dest="info", metavar='<info>', help="genome strip info pickle", required=True)
    
    parser.add_argument("-s", "--samples", dest="samples", metavar='<samples>', help="set of samples to annotate", required=True)
    
    parser.add_argument("-gender_map", "--gender_map", dest="gender_map", metavar='<gender_map>', help="gender file UUID Sex, tab delimited", required=True)
    
    
    parser.add_argument("-males_Y", "--males_Y", dest="males_Y", metavar='<males_fn>', help="specific subset of males to use for Y chrom stitching, helpful if some genomes don't include the Y chrom", required=False, default = False)   
    
    parser.add_argument("-cn_mode_male_Y_col", "--cn_mode_male_Y_col", dest="cn_mode_male_Y_col", metavar='<males_Y_fn>', help="alternative column to use for cn_mode if some subsetting of male samples is needed", required=False, default = 'cn_mode_male')   
    
    parser.add_argument("-cn_mode_male_col", "--cn_mode_male_col", dest="cn_mode_male_col", metavar='<males_y_column>', help="alternative column to use for cn_mode if some subsetting of male samples is needed", required=False, default = 'cn_mode_male')   
    
    
    
    parser.add_argument("-nj", "--no_jobs", dest="no_jobs", action = 'store_true', help="don't make jobs/job dirs", required=False)   
    
       
     
    parser.add_argument("-run_dir", "--run_dir", dest="run_dir", metavar='<run_dir>', help="directory to run the genotyping jobs from, create if non-existant", required=False, default=False)   
    
    
    parser.add_argument('-corr', '--corr', dest ='corr',  metavar='<corr_overall>', help="overall correlation of genotypes required at a site", required=False, default = 0.9, type = float)
    
    parser.add_argument('-corr_nmode', '--corr_nmode', dest ='corr_nm',  metavar='<corr_diff_mode>', help="correlation of between non-mode individuals (when > 3 individuals)", required = False, default=0.9, type = float)
    
    parser.add_argument('-dist', '--dist', dest ='dist_thresh',  metavar='<length_thresh>', help="maximum length between calls to merge", required = False, default = 10000, type = int)
        
    
    parser.add_argument("-fs", "--filter_somatic", dest="filter_somatic", action = 'store_true', help="filter previously marked somatic variants from the stitching process, save the info data without these variants", required = False)
    
    parser.add_argument("-vcf", "--vcf", dest="vcf", metavar='<vcf_fn>', help="original VCF Path", required=False)
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    parser.add_argument("-pre", "--prefix", dest="prefix", metavar='<prefix>', help="prefix to name files", default = False)
    
    parser.add_argument("-suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)


# In[85]:

# sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')

# info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_combine_ipscore_hipsci/info_all_sites_rmdup_filt.pkl')

# cns = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_combine_ipscore_hipsci/cns_all_filt.pkl')

# samples_discovery = sample_info[sample_info.CELL_TYPE != 'iPSC'].WGS_ID.tolist()
# sex_dict = sample_info.SEX.to_dict()


# males = sample_info[(sample_info.STUDY == 'iPSCORE') & (sample_info.SEX=='M')].WGS_ID.tolist()

# merged, dist_df, gs_combined_info_unannotated, to_remove, combined, st = stitch_cnvs(info, cns, sex_dict, samples_discovery, max_distance=5000, correl=0.9, corr_for_nmode = 0.7, males_Y = males, cn_mode_male_Y_col = 'cn_mode_male_ipscore_fb')

# vcf_ipscore = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/run_ipscore/results/gs_cnv.genotypes.vcf.gz'
# vcf_hipsci = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/run_hipsci/results/gs_cnv.genotypes.vcf.gz'

# tdf = combine_original_vcfs(vcf_ipscore, vcf_hipsci)

# vcf_formatted_df = vcf_format_df(merged, tdf)
# # make_genotyping_jobs2(merged_vcf, run_dir )

# run_dir = '/frazer01/projects/hipsci/pipeline/WGS/GenomeSTRiP_V3/stitching'

# write_genotyping_vcfs(vcf_formatted_df, run_dir,
#                          header_location = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_pipeline/header.txt',
#                          genotyping_job = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/notebooks/scripts/regenotype.sh')



# In[83]:

def run_from_args(args):
    cns = pd.read_pickle(args.cns)
    info = pd.read_pickle(args.info)
    samples_file = args.samples
    run_dir = args.run_dir
    males_Y_file = args.males_Y
    max_dist = args.dist_thresh
    corr_thresh = args.corr
    
    cn_mode_male_Y_col = args.cn_mode_male_Y_col
    cn_mode_male_col = args.cn_mode_male_col

#     print cn_mode_male_Y_col
#     print cn_mode_male_col
    
    # at this point we choose to leave annotated CNV that are somatic behind- won't be present in other files
    
    
    if 'somatic' in info.columns.tolist():
        
        if args.filter_somatic:
            info = info[info.somatic == False]
            

    samples = [line.rstrip() for line in open(samples_file)]
    males_Y =  [line.rstrip() for line in open(males_Y_file)]
 
    gender_file = args.gender_map
    
    sex_dict = {line.rstrip().split()[0]:line.rstrip().split()[1] for line in open(gender_file)}
    
    print 'starting stitching'
    print CM.datestring(hour=True, minute=True)
    
    
    merged, dist_df, gs_combined_info_unannotated, to_remove, combined, st = stitch_cnvs(info, cns, sex_dict, samples, max_distance=max_dist, correl= corr_thresh, corr_for_nmode = 0.7, males_Y = males_Y, cn_mode_male_Y_col = cn_mode_male_Y_col, cn_mode_male_col=cn_mode_male_col)
    
#     merged, dist_df, gs_combined_info_unannotated, to_remove, combined, st = stitch_cnvs(info, cns, sex_dict, samples, max_distance=args.length_thresh, correl=args.corr, corr_for_nmode = args.corr_nm)
    
    removed_df = pd.DataFrame(to_remove, columns=['name'])
    
    
    if not args.no_jobs:
        
    
        merged_vcf = vcf_format_df(merged, args.vcf)
        make_genotyping_jobs2(merged_vcf, run_dir, args.bam_list)


    output_location = args.output_dir
    
        
    
    
    if args.suffix:
        fn_info = os.path.join(output_location, 'gs_info' + args.suffix)
        fn_cns = os.path.join(output_location, 'gs_cns' + args.suffix)
        
        var_name_info = 'gs_info' + args.suffix
        # no need for extra suffixes for these,place in unique folders
        var_name_merge = 'stitched_sites'
        var_name_dist = 'dist_df'
        var_name_removed = 'removed_sites'
        
    
    else:
        fn_info = os.path.join(output_location, 'gs_info')
        fn_cns = os.path.join(output_location, 'gs_cns')
        var_name_info = 'gs_info_stitched_ua'
        var_name_merge = 'stitched_sites_info'
        var_name_dist = 'dist_df'
        var_name_removed = 'removed_sites'
        
        
        
    
    CM.save_dataframe(var_name_info, gs_combined_info_unannotated, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_merge, merged, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_dist, dist_df, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_removed, removed_df, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe('stitch_site_pairs_info', st, output_location, print_vars_recorded_loc=False)


# In[260]:

# pd.DataFrame(to_remove)


# In[234]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))


# In[217]:

# run_dir = DJ.root + '/private_output/gs_processing_pipeline/rd50'
# metadir = '/frazer01/home/djakubosky/BF_GS_MD'
# input_bam_list = '/frazer01/home/djakubosky/input_bam_lists/Blood_Fibroblasts_WGS.list'
# make_genotyping_jobs(merged_vcf, run_dir, metadir, input_bam_list)


# merged, dist_frame, gs_combined_info_unannotated, to_remove, combined = stitch_cnvs(gs_info, gs_cns_50, males, females, wgs_uuid_273.keys(),'_stitched', max_distance=100000)

