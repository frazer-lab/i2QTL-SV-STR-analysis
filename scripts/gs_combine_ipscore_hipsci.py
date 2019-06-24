
# coding: utf-8

# In[1]:

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
import matplotlib.ticker as ticker
import six
import networkx as nx
import scipy.stats as stats
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import  make_axes_locatable
import datetime
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
import scipy.stats as stats
from scipy.stats import mode
import argparse

import gzip

from collections import Counter
from djPyi2 import pandas_methods as pm


# # Combining iPSCORE and HipSci discovery sites
# 
#     1) Collapse redundancies
#     We will choose a single site in cases where the two sets discovered exactly the same break point, or almost exactly the same breakpoint, and genotyped it in the same individuals.  
#     - As part of this process before intersection we will consider the genotypes of sites from hipsci in ipscore samples, filtering out sites that have no non-ref calls in ipscore, as these are likely specific to hipsci samples and not redundancies to be collapsed - this step might be able to be skipped as other metrics (correlation/similarity coef) will also catch this
#     - we will then: intersect variants requiring 50% RO- we will be stringent for this part of the process, because will be stitching closely correlated variants after anyway- sites that are highly correlated and adjacent or only partially overlapping are probably more likely to represent sites that were broken up by the algorithm rather than true duplicates that we are looking for for this step
#     - we will intersect the variants, annotate the intersection bed file with info about each site, and calculate the correlation/num of differences among non-mode samples across the 478 HipSci/iPSCORE samples used for discovery (excluding hipsci and ipscore iPSC genomes because they represent redundant samples)
#     - we will go through clusters of overlapping variants (sometimes multiple calls may overlap one site heavily), and select the highest weighted pair- by average overlap percent- and look at the correlation and number of differences at the site-.  If the correlation is < 0.9 and there are differences in more than 5% of non-mode samples the sites will not be considered redundant, else, the higher GSCNQUAL site between the two will be chosen and the other marked as a redundancy.  Note: the final VCF file will contain all the sites- but include flags for redundancies and their mate sites, along with FILTER column annotations to easily subset to a mostly non-redundant call set
#     2) Stitching: Stitching will be performed similarly, on the entire set, minus the filtered sites, with clusters of correlated neighboring sites being considered a single large site, these larger sites will go on to be genotyped by genome STRiP again in all samples. 
#     
#     

# In[45]:

def make_bedtool_from_info(info, cnv_class_col):
    info['Chr'] = info.Chr.astype(str)
    info[['Start', 'End']] = info[['Start', 'End']].applymap(int)
    info.sort_values(['Chr', 'Start', 'End'], inplace=True)
    ipscore_bed_prior =  pbt.BedTool.from_dataframe(info[['Chr', 'Start', 'End', 'ID', cnv_class_col]])
    return ipscore_bed_prior


def intersect_variant_sets(ipscore_bed, gtex_bed, ro = 0.5, type_match = True):
    """given two pybedtools in the format ['chrom' 'start' 'end' 'ID' 'SVTYPE'] intersect them with RO
    and annotate what ones have matching classes """
    
    def filter_appropriate_overlaps(intersect):
        svtype1 = intersect.SVTYPE_A.tolist()
        svtype2 = intersect.SVTYPE_B.tolist()


        combos = (list(itertools.combinations(['DUP', 'mCNV'], 2)) + 
                  list(itertools.combinations(['DEL', 'mCNV'], 2)))

        acceptable_match = map(set, combos)

        data = []
        for s1,s2 in zip(svtype1, svtype2):
            if s1 == s2:
                data.append(True)
            else:
                b = set([s1, s2]) in acceptable_match
                data.append(b)
        return data 
    
    t = ['chrom', 'start', 'end', 'ID', 'SVTYPE']

    cols = ["{}_{}".format(i, 'A') for i in t] + ["{}_{}".format(i, 'B') for i in t] + ['overlap']

    # intersections
    intersect = ipscore_bed.intersect(gtex_bed, f = ro, F = ro, wo=True).to_dataframe(names = cols)
    
    # add amount of overlap (RO on each)
    try:
        
        intersect['Length_A'] = intersect.end_A.astype(int) - intersect.start_A.astype(int)
        intersect['Length_B'] = intersect.end_B.astype(int) - intersect.start_B.astype(int)

    except:
        return intersect

    intersect['RO_A'] = intersect['overlap'].astype(int)/ intersect.Length_A
    intersect['RO_B'] = intersect['overlap'].astype(int)/ intersect.Length_B
    intersect['average_RO'] = intersect[['RO_A', 'RO_B']].mean(axis = 1)
    
    if type_match:
        intersect['matching_svtypes'] = filter_appropriate_overlaps(intersect)
    else:
        intersect['matching_svtypes'] = True
    
    return intersect


def add_maf_to_intersect(intersect, info1, info2, col1, col2):
    
    intersect = intersect.merge(info2['percent_diff_from_mode'].to_frame('percent_diff_from_mode_B'), left_on='ID_B', right_index=True)
        
    intersect = intersect.merge(info1['percent_diff_from_mode'].to_frame('percent_diff_from_mode_A'), left_on='ID_A', right_index=True)

    return intersect

def intersect_variants(info1, cnv_class_col1, info2, cnv_class_col2, type_match = True, ro = 0.5, add_percent=True):
    bed1 = make_bedtool_from_info(info1, cnv_class_col1)
    bed2 = make_bedtool_from_info(info2, cnv_class_col2)
    intersect = intersect_variant_sets(bed1, bed2, ro = ro, type_match = type_match)
    
    if add_percent:
        intersect = add_maf_to_intersect(intersect, info1, info2, 'percent_diff_from_mode',
                                         'percent_diff_from_mode')
    
    return intersect

def intersect_and_get_venn(info1, cnv_class_col1, info2, cnv_class_col2, ro = 0.5, type_match=True, add_percent=True):
    intersect = intersect_variants(info1, cnv_class_col1, info2, cnv_class_col2, ro=ro, type_match = type_match, add_percent=add_percent)

    
    venn = get_sets_for_venn(intersect, info1, info2)
    return intersect, venn


def intersect_plot_venn(ax, info1, info2, cnv_class1, cnv_class2, labels, ro = 0.5, type_match = True):
    # intersect, venn sets
    int12, venn_sets = intersect_and_get_venn(info1, cnv_class1, info2, cnv_class2, ro=ro, type_match = type_match,
                                             add_percent=False)
    
    venn2(venn_sets, set_labels=labels, ax= ax)
    return int12

def get_sets_for_venn(df, info_A, info_B):
    
    ol_with_variants_A = df.ID_A.unique().shape[0]
    ol_with_variants_B = df.ID_B.unique().shape[0]
    middle = round(((ol_with_variants_A + ol_with_variants_B)/2), 0)
    left = info_A.shape[0] - middle
    right = info_B.shape[0] - middle
    
    out_dict = {'10':left, '11': middle, '01': right}
    return out_dict


# In[46]:

def get_info_cols_for_intersect_sites(intersect, info):
    """ add some info columns to the intersection df so we don't have to expensively look them up when processing"""
    intersect = intersect.copy()
    
    inds = intersect.index.tolist()
    
    
    ID_A_mod = intersect.ID_A.tolist()
    ID_A_mod = ["{}_{}".format(i, 'iPSCORE') for i in ID_A_mod]
    
    ID_B_mod = intersect.ID_B.tolist()
    ID_B_mod = ["{}_{}".format(i, 'HipSci') for i in ID_B_mod]
    
    columns_desired = ['alleles_dist_ipscore_fb', 'diploid_alleles_ipscore_fb', 'GSCNQUAL_ipscore', 'GSCNQUAL_hipsci_fib', 'diff_mode_uuids', 'percent_diff_from_mode_ipscore_fb', 'percent_diff_from_mode_hipsci_fib', 'cn_mode', 'lq_samps']
    
    name_mod_A = ["{}_{}".format(i, 'iPSCORE') for i in columns_desired]
    name_mod_B = ["{}_{}".format(i, 'HipSci') for i in columns_desired]
    
    t1 = info.loc[ID_A_mod, columns_desired]
    t1.index = inds
    t1.columns = name_mod_A
    
    
    t2 = info.loc[ID_B_mod, columns_desired]
    t2.index = inds
    t2.columns = name_mod_B
    
    tdf = pd.concat([t1, t2], axis = 1)
    
    order=['ID_A', 'ID_B', 'Length_A', 'Length_B', 'RO_A', 'RO_B', 'average_RO'] + tdf.columns.tolist()
    
    tdf['ID_A'] = ID_A_mod
    tdf['ID_B'] = ID_B_mod
    
    # pull over a few columns from intersect
    tdf = tdf.join(intersect[['Length_A', 'Length_B', 'RO_A', 'RO_B', 'average_RO']])
    tdf = tdf[order]
    
    qual_score_hip = 'GSCNQUAL_ipscore_HipSci	GSCNQUAL_hipsci_fib_HipSci'.split()
    qual_score_ipscore = 'GSCNQUAL_ipscore_iPSCORE	GSCNQUAL_hipsci_fib_iPSCORE'.split()
    tdf['GSCNQUAL_sum_HipSci_site'] = tdf[qual_score_hip].sum(axis =1)
    tdf['GSCNQUAL_sum_iPSCORE_site'] = tdf[qual_score_ipscore].sum(axis = 1)
    
    tdf['GSCNQUAL_max_HipSci_site'] = tdf[qual_score_hip].max(axis =1)
    tdf['GSCNQUAL_max_iPSCORE_site'] = tdf[qual_score_ipscore].max(axis = 1)
    
    return tdf, ID_A_mod, ID_B_mod
    


# In[47]:

def compare_lists(l1, l2):
    count = 0
    for i1, i2 in zip(l1, l2):
        if i1 != i2:
            count +=1
    return count


# In[48]:

def compare_sites(ind1, ind2, cns_t, samples, samples_nmode1, samples_nmode2, cnmode1, cnmode2, samples_lq1, samples_lq2, subtract_lq = True):
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


# In[49]:

def split_id_to_coord(ID):
    spl = ID.split('_')
    chrom = spl[1]
    start = int(spl[2])
    end = int(spl[3])
    return chrom, start, end

def compute_dist(id1, id2):
    chrom1, start1, end1 = split_id_to_coord(id1)
    chrom2, start2, end2 = split_id_to_coord(id2)
    dist = start2 - end1
    return dist
    


# In[50]:


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


# In[51]:

def annotate_redundant(df):
    
    df = df.copy()
    df['redundant'] = False
    inds = df[((df.num_diff==0) | (df.percent_non_mode_diff < 0.05)) & ((df.corr_coef > 0.95) & (df.mean_cn_diff_nmode < 0.2) & (df.mean_cn_diff_all < 0.2))].index.tolist()
    
    df['redundant'] = False
    if len(inds) > 0:
        df.loc[inds, 'redundant'] = True
    
    return df


# In[52]:

def compare_sites_intersect_v2(intersect_info, cns_t, samples, info, cn_mode_col= 'cn_mode', subtract_lq=True, chroms = CM.normal_Chrs):
    """ calculate the pearson corr of sites in the intersection, calculate the number of differences among non-mode samples """
    data = []
#     inds = intersect_info.index.tolist()
    cns_t = cns_t.loc[samples]
    diff_mode_uuids = info.diff_mode_uuids.to_dict()
    mode_cn_all = info[cn_mode_col].to_dict()
    lq_uuids = info.lq_samps.to_dict()
    cnv_classes = info.cnv_class.to_dict()

    inds = []
    for ind, x in intersect_info.iterrows():
        
        ind1 = x['ID_A']
        # hipsci
        ind2 = x['ID_B']
        RO = x['average_RO']
        
        ind_ipscore = ind1
        ind_hipsci = ind2
        chrom, start, end = split_id_to_coord(ind1)
        
        if chrom in chroms:
            inds.append(ind)
            

            qual_col = 'GSCNQUAL_sum_HipSci_site	GSCNQUAL_sum_iPSCORE_site'.split()

            qual_hipsci = x[qual_col[0]]
            qual_ipscore = x[qual_col[1]]

            quals = [qual_ipscore, qual_hipsci]
            max_qual = max(quals)

            # ipscore
            ind1 = x['ID_A']
            # hipsci
            ind2 = x['ID_B']



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

            best_qual_site = ind_ipscore
            worst_qual_site = ind_hipsci
            if quals.index(max_qual) != 0:
                best_qual_site = ind_hipsci
                worst_qual_site = ind_ipscore            





            comp = comp + [pair, mode_cn1, mode_cn2, cnv_class1, cnv_class2, distance_between, 
                           absolute_dist, pair_ind, num_diff1, num_diff2, diff_uuids1, diff_uuids2,
                          chrom, best_qual_site, worst_qual_site, RO, qual_ipscore, qual_hipsci]
            data.append(comp)
    
    
    df = pd.DataFrame(data, columns=['ID_iPSCORE', 'ID_HipSci', 'corr_coef',
                             'num_diff', 'num_non_mode', 'percent_non_mode_diff', 
                             'samps_to_compare_nmode', 'samps_to_compare_corr',
                             'samps_to_exclude', 'num_pass', 'exact_cn_match', 
                             'mean_cn_diff_all', 'mean_cn_diff_nmode','alleles1', 'alleles2', 
                             'num_alleles1', 'num_alleles2','allele_dist1', 'allele_dist2','pair', 
                             'mode_cn1', 'mode_cn2', 'cnv_class1', 'cnv_class2', 
                             'distance_between', 'distance_between_mod', 'cat_pair', 
                             'num_diff1', 'num_diff2', 'diff_uuids1', 'diff_uuids2', 'chrom', 'hq_site', 
                                     'lq_site', 'average_RO', 'GSCNQUAL_sum_iPSCORE_site','GSCNQUAL_sum_HipSci_site'],
                      index = inds)
    df = df.pipe(annotate_redundant)
    return df


# In[53]:

def gather_consensus_X(adj_sites_all_males, adj_sites_all_females, cns_t, bool_col ='redundant'):
    inds = adj_sites_all_males.index.tolist()
    assert (inds == adj_sites_all_females.index.tolist())
    
    num_non_mode_males = adj_sites_all_males.num_non_mode.tolist()
    num_non_mode_females = adj_sites_all_females.num_non_mode.tolist()
    
    ndiff_male = adj_sites_all_males.num_diff.tolist()
    ndiff_female = adj_sites_all_females.num_diff.tolist()
    
    passing_b_males = adj_sites_all_males[bool_col].tolist()
    passing_b_females = adj_sites_all_females[bool_col].tolist()
    
    samples_nm_males = adj_sites_all_males.samps_to_compare_corr.tolist()
    samples_nm_females = adj_sites_all_females.samps_to_compare_corr.tolist()
    union_samples = [list(set(i1 + i2)) for i1, i2 in zip(samples_nm_males, samples_nm_females)]
   
    corr_males = adj_sites_all_males.corr_coef.tolist()
    corr_females = adj_sites_all_females.corr_coef.tolist()
    dist = adj_sites_all_males.distance_between_mod.tolist()
    
    out = []
    pairs = adj_sites_all_males.cat_pair.tolist()
    for ind, nm_male, nm_female, pass_b_male, pass_b_female, cm, cf, ndm, ndf, d, us, p in zip(inds,
                                                                                            num_non_mode_males,
                                                                           num_non_mode_females,
                                                                           passing_b_males, passing_b_females, 
                                                                           corr_males, corr_females,
                                                                           ndiff_male, ndiff_female, 
                                                                           dist, union_samples, pairs):
        
        
        ind1, ind2 = p.split('-')
        chrom = ind1.split('_')[1]
        
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
                    out.append([ind,False])
                    
            else:
                print "didn't account for all scenarios"
                break
                
    df = pd.DataFrame(out, columns=['numeric_ind', 'redundant_consensus'])
    df.index = df.numeric_ind
    return df


# In[54]:

def annotate_consensus_X(int_info, int_info_males, int_info_females, int_info_males_ipscore, x_consensus):
    int_info['redundant_adjusted'] = int_info.redundant
    inds_x = int_info[int_info.chrom=='X'].index.tolist()
    inds_y = int_info[int_info.chrom=='Y'].index.tolist()
    int_info_y_males = int_info_males_ipscore.copy()

    int_info.loc[inds_x, 'redundant_adjusted'] = False
    int_info.loc[inds_y, 'redundant_adjusted'] = False

    int_info.loc[inds_x, 'redundant_adjusted'] = x_consensus.redundant_consensus
    
    if int_info_males_ipscore.shape[0] > 0:
        int_info.loc[inds_y, 'redundant_adjusted'] = int_info_y_males.redundant
    return int_info
    


# In[55]:

def get_intersect_info(intersect, cns_t, info_all_sites, samples):
    intersect_info, ID_A_mod, ID_B_mod = intersect.pipe(get_info_cols_for_intersect_sites, info_all_sites)
    intersect_info = intersect_info.pipe(compare_sites_intersect, cns_t, samples)
    return intersect_info, ID_A_mod, ID_B_mod


# In[56]:

def annotate_clusters_exact_ol(intersect_info, info_all_sites):
    # mark sites that exactly overlap
    info_all_sites['exact_overlap'] = False
    
    data_clusters = []
    inds = []
    for x, df in intersect_info.groupby('ID_iPSCORE'):
        cluster = df.ID_HipSci.tolist()
        cluster.append(x)
        cluster = ",".join(sorted(cluster))
        
        data_clusters.append(cluster)
        data_clusters.append(cluster)
        inds.append(x)
        inds.append(df.ID_HipSci.values[0])
    
    cluster_dict_exact = dict(zip(inds, data_clusters))

    info_all_sites.loc[inds, 'intersect_cluster'] = data_clusters
    info_all_sites.loc[inds, 'intersect_cluster_size'] = 2
    info_all_sites.loc[inds,'exact_overlap'] = True
    return info_all_sites, cluster_dict_exact


# In[57]:

def filter_non_exact_ol(other_sites):
    # pick a variant from each cluster to represent it
    
    data = []
    non_primary_sites = []
    primary_site = []
    
    Z = nx.Graph()
    
    
    for i, x in other_sites.iterrows():
        id_a = x['ID_iPSCORE']
        id_b = x['ID_HipSci']
        corr_coef = x['corr_coef']
        RO = x['average_RO']
        perc_diff = x['percent_non_mode_diff']
        hq_site = x['hq_site']
        lq_site = x['lq_site']
        redundant = x['redundant_adjusted']
        qual_ipscore = x['GSCNQUAL_sum_iPSCORE_site']
        qual_hipsci = x['GSCNQUAL_sum_HipSci_site']

        Z.add_edge(id_a, id_b, RO = RO, correlation = corr_coef, perc_diff = perc_diff,
                   hq_site = hq_site, lq_site=lq_site, qual_ipscore = qual_ipscore,qual_hipsci=qual_hipsci)
        
    
    
    t = nx.connected_component_subgraphs(Z)
    cluster_dict = {}
    for g in t:
        node_pairs = g.edges()
        if len(node_pairs) == 1:
         
            max_pair = sorted(node_pairs[0])
            p = node_pairs[0]
            
            ed = g.get_edge_data(p[0], p[1])
            hq_site = ed['hq_site']
            lq_site = ed['lq_site']
            
            cluster = ",".join(max_pair)
            cluster_dict[hq_site] = cluster
            cluster_dict[lq_site] = cluster
            data.append([cluster, cluster, hq_site, lq_site, ed, 2])
        else:
            max_weight = 0
            max_pair = ()
            for np in node_pairs:
                ed = g.get_edge_data(np[0], np[1])
                ro = ed['RO']
                if ro > max_weight:
                    max_weight = ro
                    max_pair = np
            ed= g.get_edge_data(max_pair[0], max_pair[1])
            hq_site = ed['hq_site']
            lq_site = ed['lq_site']
            
            cluster_nodes = sorted(list(set(CM.flatten_list(node_pairs))))
            
            
            
            cluster_size = len(cluster_nodes)
            cluster = ",".join(cluster_nodes)
            
            for s in cluster_nodes:
                cluster_dict[s] = cluster
                
            
           
            other_lq_sites = CM.flatten_list([i for i in node_pairs if i != max_pair])
            other_lq_sites = [i for i in other_lq_sites if i != hq_site] + [lq_site]
            other_lq_sites = list(set(other_lq_sites))
            other_lq_sites = ",".join(other_lq_sites)
            
            max_pair = ",".join(sorted(max_pair))
            data.append([cluster, max_pair, hq_site, other_lq_sites, ed, cluster_size])
    
    df = pd.DataFrame(data, columns=['cluster', 'max_weighted_pair', 'hq_site', 'lq_sites', 'ed', 'cluster_size'])
    df = df.sort_values('cluster')
    return df, cluster_dict
    


# In[58]:

def lambda_add_clusters(x, cluster_dict):
    l = x['intersect_cluster']
    ind = x['ID']
    if l == False:
        t = cluster_dict.get(ind, False)
        return t
    else:
        return l


# In[59]:

def lambda_add_clusters_size(x):
    if x != False:
        s = len(x.split(','))
        return s
    else:
        return 0


# In[60]:

# def add_cluster_IDs(info_all_sites):
#     info_all_sites = info_all_sites.copy()
#     info_all_sites = info_all_sites.sort_values(['Chr', 'Start', 'End'])

# #     info_all_sites['intersect_cluster_ID'] = 0
#     has_cluster = info_all_sites[info_all_sites.intersect_cluster != False]

#     cluster_id_dict = {}
#     clusters = []
#     clust_number = 0
#     for ind, x in has_cluster.iterrows():
#         clust = x['intersect_cluster']

#         if clust not in clusters:
#             clusters.append(clust)
#             clust_number += 1
#             cluster_id_dict[clust] = clust_number

#     info_all_sites['intersect_cluster_ID'] = info_all_sites.intersect_cluster.apply(lambda x: cluster_id_dict.get(x, 0))
    
    
    
#     primary_site_dict = {}
#     for ind,x in has_cluster[has_cluster.primary_site==True].iterrows():
#         clust = x['intersect_cluster']
#         spl = clust.split(',')
#         for x in spl:
#             primary_site_dict[x] = ind
#         primary_site_dict[x] = x
    

#     info_all_sites['intersect_cluster_primary_cnv_ID'] = info_all_sites.ID.apply(lambda x: primary_site_dict.get(x, False))
    
#     return info_all_sites


# In[ ]:

def add_cluster_IDs(info):
    info = info.copy()
    has_cluster = info[info.intersect_cluster != False]
    has_cluster = has_cluster.sort_values(['Chr', 'Start', 'End'])
    cluster_id_dict = {}
    clusters = []
    clust_number = 0
    for ind, x in has_cluster.iterrows():
        clust = x['intersect_cluster']

        if clust not in clusters:
            clusters.append(clust)
            clust_number += 1
            cluster_id_dict[clust] = clust_number

    primary_site_dict = {}
    for ind,x in has_cluster[has_cluster.primary_site==True].iterrows():
        clust = x['intersect_cluster']
        spl = clust.split(',')
        for x in spl:
            primary_site_dict[x] = ind
        primary_site_dict[x] = x


    info['intersect_cluster_primary_cnv_ID'] = info.ID.apply(lambda x: primary_site_dict.get(x, False))

    info['intersect_cluster_ID'] = info.intersect_cluster.apply(lambda x: cluster_id_dict.get(x, 0))

    int_cluster_dict = has_cluster.intersect_cluster.to_dict()

    info['intersect_cluster_primary'] = info.ID.apply(lambda x: int_cluster_dict.get(x, False))
    return info


# In[61]:

def filter_redundancies(intersect_info, info_all_sites, info_filt):
    
    
    intersect_info['ID_A'] = intersect_info.ID_iPSCORE.apply(lambda x: x.replace('_iPSCORE', ''))
    intersect_info['ID_B'] = intersect_info.ID_HipSci.apply(lambda x: x.replace('_HipSci', ''))

    exact_overlaps = intersect_info[intersect_info.ID_A == intersect_info.ID_B]
    
    
    # first annotate 
    info_all_sites['intersect_cluster'] = False
    info_all_sites['intersect_cluster_size'] = 0
    info_all_sites, cluster_dict_exact = annotate_clusters_exact_ol(exact_overlaps, info_all_sites)

    primary_sites = exact_overlaps.hq_site.tolist()
    non_primary_sites = exact_overlaps.lq_site.tolist()
    
    
    inds = info_filt.index.tolist()
    info_all_sites['primary_site'] = False
    info_all_sites.loc[inds, 'primary_site'] = True
    info_all_sites.loc[non_primary_sites, 'primary_site'] = False



    sites_A = exact_overlaps.ID_A.tolist()
    sites_B = exact_overlaps.ID_B.tolist()
    # sites that are exactly overlapped will be removed from other redundancy comparisons
    # remaining are non-exact overlaps that need to be checked for redundancy

    other_sites = intersect_info[~((intersect_info.ID_A.isin(sites_A))| (intersect_info.ID_B.isin(sites_B)))]
    other_sites = other_sites[other_sites.redundant_adjusted == True]
    
    tdf, cluster_dict = filter_non_exact_ol(other_sites)
    primary_sites = tdf['hq_site'].tolist()
    non_primary_sites = tdf.lq_sites.tolist()

    non_primary_sites = list(set(CM.flatten_list([i.split(',') for i in non_primary_sites])))
    info_all_sites.loc[non_primary_sites, 'primary_site'] = False
    info_all_sites['intersect_cluster'] = info_all_sites.apply(lambda x: lambda_add_clusters(x, cluster_dict), axis =1)
    info_all_sites['intersect_cluster_size'] = info_all_sites.intersect_cluster.apply(lambda x: lambda_add_clusters_size(x))
    
    info_all_sites = info_all_sites.pipe(add_cluster_IDs)
    return info_all_sites, tdf, exact_overlaps


# In[12]:

# info_ipscore,  filt_ipscore_ind, unfilt_ipscore = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/ipscore/gs_info_annot.pkl').pipe(prep_info).pipe(filter_info, qual_thresh_dup=12, qual_thresh_del=2)


# In[13]:

# info_hipsci_ipscore, filt_h_ipscore_ind, unfilt_h_ipscore = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/ipscore_h_rg/gs_info_annot.pkl').pipe(prep_info).pipe(filter_info, 15,2,20)


# In[62]:

def prep_info_callset(df, info):
    """make sure that info dfs have only sites that haven't been filtered in the info_all"""
    df = df.copy()
    df = df[df.ID.isin(info.ID_original.tolist())]
    return df


# In[63]:

def prep_info(df):
    df = df.copy()
    df = df[df.cnv_class != 'Non_Bi']
    df = df[df.FILTER_NUM_LQ == False]
    return df


# In[64]:

def filter_info(df, qual_thresh_dup=12, qual_thresh_del=2, qual_thresh_mcnv=14, num_lq_thresh=1000, qual_col = 'GSCNQUAL', prep_info_pipe=True):
    df = df.copy()
    tdf =df.copy()
    
    if prep_info_pipe:
        
        df = df.pipe(prep_info)
    
    df = df[df.somatic == False]
    inds_prefilt = copy.deepcopy(df.index.tolist())
    df[qual_col] = df[qual_col].astype(float)
    df = df[df.cnv_class != 'Non_Bi']
    
    
    df = df[~((df.cnv_class == 'DUP') & (df[qual_col] < qual_thresh_dup))]
    df = df[~((df.cnv_class == 'DEL') & (df[qual_col]< qual_thresh_del))]
    df = df[~((df.cnv_class == 'mCNV') & (df[qual_col] < qual_thresh_mcnv))]
#     df = df[df.num_union_bf <= num_lq_thresh]
    
    
    inds = df.ID.tolist()
    s1 = set(inds)
    s2 = set(inds_prefilt)
    filtered = list(s2.difference(s1))
    tdf['FILTER_GSCNQUAL'] = False
    tdf.loc[filtered, 'FILTER_GSCNQUAL'] = True
    
    
    return df, filtered, tdf


# In[65]:

def prep_cns_all_t(cns_all, info_all_sites):
    """ prep and filter transposed cns matrix"""
    
    cns_all = cns_all.copy()
    cns_all.drop('old_index', axis =1, inplace=True)
    cns_all = cns_all.reindex(info_all_sites.index)
    
    cns_t = cns_all.T.copy()
    return cns_t
    


# In[66]:

def prep_info_all_sites(info_all_sites):
    """ prepares info for intersect- can be run directly on the read file to produce the filtered and full info dfs"""
    
    info_all_sites = info_all_sites.copy()
    info_all_sites['GSCNQUAL_max'] = info_all_sites[['GSCNQUAL_ipscore','GSCNQUAL_hipsci_fib']].max(axis =1)

    info_filt,  filt_all_ind, info_all_sites = info_all_sites.pipe(filter_info, qual_thresh_dup=12, qual_thresh_del=2, qual_thresh_mcnv=14, qual_col = 'GSCNQUAL_max')
    
    
    return info_all_sites, info_filt


# In[3]:

def run_intersect_combine_pipeline(info_ipscore, info_hipsci_ipscore, info_all_sites, info_filt, sample_info, cns_t):

    samples_discovery = sample_info[sample_info.CELL_TYPE != 'iPSC'].WGS_ID.tolist()
    samples_ipscore_males = sample_info[(sample_info.STUDY == 'iPSCORE') & (sample_info.SEX == 'M')].WGS_ID.tolist()
    samples_females = sample_info[(sample_info.CELL_TYPE != 'iPSC') & (sample_info.SEX=='F')].WGS_ID.tolist()
    samples_males = sample_info[(sample_info.CELL_TYPE != 'iPSC') & (sample_info.SEX=='M')].WGS_ID.tolist()
    
    intersect = intersect_variants(info_ipscore, 'cnv_class', info_hipsci_ipscore, 'cnv_class', type_match=False, add_percent=True)
    
    
    intersect_info, ID_A_mod, ID_B_mod = intersect.pipe(get_info_cols_for_intersect_sites, info_all_sites)
    int_info = compare_sites_intersect_v2(intersect_info, cns_t, samples_discovery, info_all_sites)
    
    
    int_info_males = compare_sites_intersect_v2(intersect_info, cns_t, samples_males,
                                            info_all_sites, chroms =['X'],
                                            cn_mode_col='cn_mode_male', subtract_lq=True)
    int_info_females = compare_sites_intersect_v2(intersect_info, cns_t, samples_males,
                                              info_all_sites, chroms=['X'], 
                                              cn_mode_col='cn_mode_female', subtract_lq=True)
    int_info_males_ipscore = compare_sites_intersect_v2(intersect_info, cns_t, samples_males, info_all_sites, 
                                                    chroms=['Y'], cn_mode_col='cn_mode_male_ipscore_fb',
                                                    subtract_lq=True)
    
    x_consensus = gather_consensus_X(int_info_males, int_info_females, cns_t)
    int_info = annotate_consensus_X(int_info, int_info_males, int_info_females, int_info_males_ipscore, x_consensus)
    info_all_sites, cluster_info, exact_overlaps = filter_redundancies(int_info, info_all_sites, info_filt)
    
    int_info.index = int_info.cat_pair

    int_info['exact_overlap_pair'] = False
    
    inds = int_info.loc[exact_overlaps.cat_pair.tolist(), 'exact_overlap_pair'] = True
    
    return info_all_sites, int_info, cluster_info


# In[2]:

# def fix_cluster_annotations(info):
#     info = info.copy()
#     has_cluster = info[info.intersect_cluster != False]
#     has_cluster = has_cluster.sort_values(['Chr', 'Start', 'End'])
#     cluster_id_dict = {}
#     clusters = []
#     clust_number = 0
#     for ind, x in has_cluster.iterrows():
#         clust = x['intersect_cluster']

#         if clust not in clusters:
#             clusters.append(clust)
#             clust_number += 1
#             cluster_id_dict[clust] = clust_number

#     primary_site_dict = {}
#     for ind,x in has_cluster[has_cluster.primary_site==True].iterrows():
#         clust = x['intersect_cluster']
#         spl = clust.split(',')
#         for x in spl:
#             primary_site_dict[x] = ind
#         primary_site_dict[x] = x


#     info['intersect_cluster_primary_cnv_ID'] = info.ID.apply(lambda x: primary_site_dict.get(x, False))

#     info['intersect_cluster_ID'] = info.intersect_cluster.apply(lambda x: cluster_id_dict.get(x, 0))

#     int_cluster_dict = has_cluster.intersect_cluster.to_dict()

#     info['intersect_cluster_primary'] = info.ID.apply(lambda x: int_cluster_dict.get(x, False))
#     return info


# In[97]:

def filter_all_sites(info_all_sites):
    info_all_sites = info_all_sites.copy()
    
    
    info_all_sites = info_all_sites[(info_all_sites.somatic == False) & (info_all_sites.cnv_class != 'Non_Bi') & (info_all_sites.FILTER_GSCNQUAL == False) & (info_all_sites.primary_site ==True) & (info_all_sites.FILTER_NUM_LQ==False)]
    return info_all_sites


# In[ ]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-cns", "--cns", dest="cns", metavar='<cns_all>', help="genome strip cns pickle (Copy Number States)", required=True)
    
    parser.add_argument("-info", "--info", dest="info", metavar='<info_all>', help="genome strip info pickle", required=True)
    
    parser.add_argument("-info_ipscore", "--info_ipscore", dest="info_ipscore", metavar='<info_ipscore>', help="genome strip info pickle", required=True)
    
    parser.add_argument("-info_hipsci_ipscore", "--info_hipsci_ipscore", dest="info_hipsci_ipscore", metavar='<info_hipsci_ipscore', help="genome strip info pickle", required=True)
    
    parser.add_argument("-sample_info", "--sample_info", dest="sample_info", metavar='<sample_info>', help="sample info file", required=True)

    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    

    parser.add_argument("-file_suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line tool to combine ipscore and hipsci genome strip variants')
    
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    
    info_all_sites, info_filt = pd.read_pickle(args.info).pipe(prep_info_all_sites)
    info_ipscore = pd.read_pickle(args.info_ipscore).pipe(prep_info_callset, info_filt)
    info_hipsci_ipscore= pd.read_pickle(args.info_hipsci_ipscore).pipe(prep_info_callset, info_filt)
    sample_info = pd.read_pickle(args.sample_info)
    
    cns_t  = pd.read_pickle(args.cns).pipe(prep_cns_all_t, info_filt)
    
    print 'Intersecting and Marking Duplicate Sites in iPSCORE and HipSci'
    print CM.datestring(hour=True, minute=True)
   
    
    info_all_sites, int_info, cluster_info = run_intersect_combine_pipeline(info_ipscore, info_hipsci_ipscore, 
                                                             info_all_sites, info_filt, sample_info, cns_t)
    
#     info_all_sites = info_all_sites.pipe(fix_cluster_annotations)
    info_filt = info_all_sites.pipe(filter_all_sites)
#     int_info.redundant_adjusted.value_counts()
    
  
    combined_exact = int_info[int_info.exact_overlap_pair ==True].shape[0]
    exact_ols = combined_exact * 2
    other_sites = cluster_info.cluster_size.sum()
    combined_info = cluster_info.shape[0]
    
    print "{} exact overlaps combined into {}".format(exact_ols, combined_exact)
    print "{} non exact overlaps combined into {}".format(other_sites, combined_info)
    
    total_combined = exact_ols + other_sites
    total_out = combined_info + combined_exact
    
    print "{} total sites combined into {}".format(total_combined, total_out)
    
    
    
    
    output_location = args.output_dir
    if args.suffix:
        fn_info = os.path.join(output_location, 'gs_info' + args.suffix)      
        var_name_info = 'gs_info' + args.suffix
        var_name_info_filt = 'gs_info_filt' + args.suffix
        var_name_intersect = 'intersect_info'
        var_name_clust_intersect = 'intersect_cluster_info'
       
       
    else:
        fn_info = os.path.join(output_location, 'gs_info')
        fn_cns = os.path.join(output_location, 'gs_cns')
        var_name_info = 'gs_info' 
        var_name_info_filt = 'gs_info_filt' 
        var_name_intersect = 'intersect_info'
        var_name_clust_intersect = 'intersect_cluster_info'
       
    
    print 'data annotated'
    print CM.datestring(hour=True, minute=True)
    CM.save_dataframe(var_name_info, info_all_sites, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_info_filt, info_filt, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_intersect, int_info, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_clust_intersect, cluster_info, output_location, print_vars_recorded_loc=False)
    

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))


# In[96]:

# info_all_sites, info_filt = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/i2QTL_combined/gs_info_combined_annot_lq.pkl').pipe(prep_info_all_sites)


# info_ipscore = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/ipscore/gs_info_annot.pkl').pipe(prep_info_callset, info_filt)

# info_hipsci_ipscore = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/ipscore_h_rg/gs_info_annot.pkl').pipe(prep_info_callset, info_filt)


# sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')


# cns_t = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/i2QTL_combined/cns_all.pkl').pipe(prep_cns_all_t, info_filt)

# info_all_sites, int_info= run_intersect_combine_pipeline(info_ipscore, info_hipsci_ipscore, info_all_sites, info_filt, sample_info, cns_t)

# info_filt = info_all_sites.pipe(filter_all_sites)

# outdir = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/i2QTL_combined'
# # DJ.makedir(outdir)


# In[103]:

# CM.save_dataframe('info_all_sites_rmdup', info_all_sites, outdir)

# CM.save_dataframe('info_all_sites_rmdup_filt', info_filt, outdir)

# cns_filt = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V4/i2QTL_combined/cns_all.pkl').reindex(info_filt.index.tolist())

# CM.save_dataframe('cns_all_filt', cns_filt, outdir)

