
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
from collections import Counter

from scipy.stats import mode
from djPyBio import Common as CM



# In[2]:

def get_mode_col(tran, uuids):
    modes = []
    nnrefs = []

    for c in tran.columns:
        l = tran[c].tolist()
        m = int(mode(l).mode[0])
        
        nnref = 0
        count = 0
        nref_uuids = []
        for c in l:
            
            if c <> m:
                nnref +=1
                nref_uuids.append(uuids[count])
            count +=1
        
        nnrefs.append(nnref)
        modes.append(m)
        
        
    return modes, nnref, nref_uuids
    


# In[3]:

def safe_div(x, y, alt=0):
    try:
        return x/y
    except:
        return alt


# In[4]:

def transform_to_biallelic(l_dict, uuids_to_test, chrom, sex_dict, svtype = 'DEL'):
    
    if chrom not in ['X', 'Y']:
        
        del_dict = {0:'1/1', 1: '0/1', 2: '0/0'}
        dup_dict = {2:'0/0', 3: '0/1' , 4: '1/1'}
        trans_dict = {}
        if svtype == 'DUP':
            for s in uuids_to_test:
                trans_dict[s] = dup_dict[int(l_dict[s])]
        if svtype == 'DEL':
            for s in uuids_to_test:
                trans_dict[s] = del_dict[int(l_dict[s])]
        
    
    else:
        
        conv_dict = {'DEL': {'F': {0:'1/1', 1: '0/1', 2: '1/1'}, 'M': {0:'0/1', 1: '0/0'}}, 
                    'DUP': {'F': {2:'0/0', 3: '0/1' , 4: '1/1'}, 'M':  {1:'0/0', 2: '0/1'}}}
        
        if chrom == 'X':
            for s in uuids_to_test:
                sex = sex_dict[s]
                trans_dict[s] = conv_dict[svtype][sex][int(l_dict[s])]
                
        
        if chrom == 'Y':
            for s in uuids_to_test:
                sex = sex_dict[s]
                if sex == 'M':
                    trans_dict[s] = conv_dict[svtype][sex][int(l_dict[s])]
                else:
                    trans_dict[s] = './.'
                      
    return trans_dict


# In[5]:

def compute_maf_lightweight(gts_dict, chrom, sex_dict):
    
    c = Counter(gts_dict.values())
    if chrom not in ['X', 'Y']:
    
#         gts_dict_convert = {'0/0':0, '0/1':1, '1/1':2, './.':0}

        ref_allele = (c.get('0/0', 0) * 2) + (c.get('0/1', 0) * 1)
        alt_allele = (c.get('0/1', 0) * 1) + (c.get('1/1', 0) * 2) 
        tot_alleles = ref_allele + alt_allele

#         ref_af = safe_div(ref_allele, tot_alleles, alt='All Missing')
#         alt_af = safe_div(alt_allele, tot_alleles, alt='All Missing')
#         afs = [ref_af, alt_af]

            
    else:
        male_gts = [gts_dict[s] for s in gts_dict.keys() if sex_dict[s] == 'M']
        female_gts =  [gts_dict[s] for s in gts_dict.keys() if sex_dict[s] == 'F']
    
        c_male = Counter(male_gts)
        c_female = Counter(female_gts)
        
        
        
        ref_allele_f = (c_female.get('0/0', 0) * 2) + (c_female.get('0/1', 0) * 1) 
        alt_allele_f = (c_female.get('0/1', 0) * 1) + (c_female.get('1/1', 0) * 2) 
        
        ref_allele_m = (c_male.get('0/0', 0) * 1) + (c_male.get('0/1', 0) * 1) 
        alt_allele_m = (c_male.get('0/1', 0) * 1)
        
        ref_allele = ref_allele_f + ref_allele_m 
        alt_allele = alt_allele_f + alt_allele_m 
        
        tot_alleles = ref_allele + alt_allele
    
    
    ref_af = safe_div(ref_allele, tot_alleles, alt='All Missing')
    alt_af = safe_div(alt_allele, tot_alleles, alt='All Missing')
    afs = [ref_af, alt_af]


    maf = 0.0
    minor_allele = 'ALT'
    if not tot_alleles == 0:
        maf = min(afs)
        if afs.index(maf) == 0:
            minor_allele = 'REF'
        else:
            minor_allele = 'ALT'


    if minor_allele == 'ALT':
        non_ref = c.get('0/1', 0) + c.get('1/1', 0)
        ref = c.get('0/0', 0)
    if minor_allele == 'REF':
        non_ref = c.get('0/1', 0) + c.get('0/0', 0)
        ref = c.get('1/1', 0)      
        
    out_names = ['NNREF', 'ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF']
    out_data = [non_ref, alt_af, ref_af, maf, minor_allele, ref]
    data_dict = dict(zip(out_names, out_data))
    
    return out_data, data_dict


# In[6]:

def estimate_variant_class(min_allele, max_allele, ref):
    
    if min_allele == max_allele:
        return "Non_Bi"
    
    if (min_allele < ref) and (max_allele > ref):
        return 'MIXED'
    
    if (max_allele > ref) and (min_allele ==ref):
        return 'DUP'
    if (min_allele < ref) and( max_allele == ref):
        return 'DEL'
    
    if (min_allele > ref) and (max_allele > ref) and (min_allele != max_allele):
        
        return 'DUP'
    
    if (min_allele < ref) and (max_allele < ref):
        return 'DEL'


def minimum_predicted_alleles(min_in, max_in, diploid = True):
    """compute minimum predicted haploid 
    alleles that can explain a range of CN between the min and max diploid CN 
    that could exist in a population """
    
    if diploid == True:
        
        if min_in==max_in:
            num_alleles_predicted = 1

        elif min_in >= 1:
            num_alleles_predicted = 1 + round(max_in/2) - round((min_in-1)/2)

        else:
            num_alleles_predicted = round(max_in/2) - round((min_in-1)/2)
            
    else:
        # if a haploid situation, assume we have all the individual alleles in the range
        if min_in==max_in:
            num_alleles_predicted = 1
        else:
            num_alleles_predicted = (max_in - min_in) + 1
    
    return num_alleles_predicted


def get_variant_type_info(min_allele, max_allele, min_allele_males, max_allele_males, min_allele_females, 
             max_allele_females, chrom):
    
    ref = 2
    ref_males = 1
    
    if chrom =='Y':
        variant_class = estimate_variant_class(min_allele_males, max_allele_males, ref_males)
        min_pred_alleles= minimum_predicted_alleles(min_allele_males, max_allele_males, diploid= False)
        
    elif chrom == 'X':
        variant_class_males = estimate_variant_class(min_allele_males, max_allele_males, ref_males)
        variant_class_females = estimate_variant_class(min_allele_females, max_allele_females, ref)
        
        min_pred_alleles_males = minimum_predicted_alleles(min_allele_males, max_allele_males, diploid= False)
        min_pred_alleles_females = minimum_predicted_alleles(min_allele_females, max_allele_males)
        
        min_pred_alleles = max(min_pred_alleles_males, min_pred_alleles_females)
        
        vts = set([variant_class_males, variant_class_females])
        if len(vts) > 1:
            MIXED= set(['DUP', 'DEL'])
            if 'MIXED' in vts:
                variant_class = 'MIXED'
            elif vts == MIXED:
                variant_class ='MIXED'
            
            else:
                vts.remove('Non_Bi')
                variant_class = vts.pop()
            
        else:
            variant_class = vts.pop()
       
    else:
        variant_class = estimate_variant_class(min_allele, max_allele, ref)
        min_pred_alleles = minimum_predicted_alleles(min_allele, max_allele)
    
    
    cnv_class = copy.deepcopy(variant_class)
    
    if (min_pred_alleles >= 3):
        cnv_class = 'mCNV'
    if variant_class == 'MIXED':
        cnv_class = 'mCNV'
    
    
    
    return variant_class, cnv_class, min_pred_alleles


# In[7]:

def check_biallelic_status(unique_gts_site, unique_gts_males, unique_gts_females, cnv_class, chrom):
    if chrom not in ['X', 'Y']:
        bi_class = is_biallelic_DUP_DEL_auto(unique_gts_site, cnv_class)
    else:
        bi_class = is_biallelic_DUP_DEL_xy(unique_gts_males, unique_gts_females, chrom)
        
    return bi_class

def is_biallelic_DUP_DEL_xy(unique_gts_males, unique_gts_females, chrom):
    
    m_svtype = 'UNKNOWN'
    f_svtype = 'UNKNOWN'

    DUP_alleles_F = [2,3,4]
    DEL_alleles_F = [0,1,2]
    num_unique_F = len(unique_gts_females)

    DUP_alleles_M = [1, 2]
    DEL_alleles_M = [0,1]
    num_unique_M = len(unique_gts_males)


    # is it biallelic in M/F
    if num_unique_M > 2:
        m_svtype = 'mCNV'

    elif num_unique_F > 3:
        f_svtype = 'mCNV'

    intersect_M_dup  = set(unique_gts_males).intersection(DUP_alleles_M)
    intersect_M_del  = set(unique_gts_males).intersection(DEL_alleles_M)


    intersect_F_dup  = set(unique_gts_females).intersection(DUP_alleles_F)
    intersect_F_del  = set(unique_gts_females).intersection(DEL_alleles_F)
    
    if chrom == 'X':

        if len(intersect_M_dup) == num_unique_M:
            m_svtype == 'DUP'

        if len(intersect_M_del) == num_unique_M:
            m_svtype == 'DEL'

        if len(intersect_F_dup) == num_unique_F:
            f_svtype == 'DUP'

        if len(intersect_F_del) == num_unique_F:
            f_svtype == 'DEL'

        if f_svtype == m_svtype:
            return f_svtype

        else:
            return 'UNKNOWN'
        
    if chrom == 'Y':
        
        if len(intersect_M_dup) == num_unique_M:
            m_svtype == 'DUP'
            
        if len(intersect_M_del) == num_unique_M:
            m_svtype == 'DEL'
        
        return m_svtype

def is_biallelic_DUP_DEL_auto(unique_gts_site, cnv_class):
    """ for autosomes """
    DUP_alleles = [2,3,4]
    DEL_alleles = [0,1,2]
    num_unique = len(unique_gts_site)
    
    
    if num_unique > 3:
        return cnv_class
    
    if cnv_class == 'DUP':
        intersect = set(unique_gts_site).intersection(DUP_alleles)
        if len(intersect) == num_unique:
            return 'DUP'
    
    elif cnv_class == 'DEL':
        intersect = set(unique_gts_site).intersection(DEL_alleles)
        if len(intersect) == num_unique:
            return 'DEL'  
    else:
        return cnv_class
    


# In[8]:

def gs_annotations(tran, uuids, sex_dict, lq_dict = False):
    
    modes_col = []
    modes_male_col = []
    modes_female_col = []
    nnrefs_col = []
    unique_gts_col = []
    alleles_dist_col = []
    alleles_dist_males_col = []
    alleles_dist_females_col = []
    
    
    variants_discovered_col = []
    variant_allele_count_col = []
    count_alleles_col = []
    min_allele_col = []
    max_allele_col = []
    
    min_allele_males_col = []
    max_allele_males_col = []
    
    min_allele_females_col = []
    max_allele_females_col = []
    
    is_mixed_col = []
    percent_diff_from_mode_col = []
    
    nref_uuids_col = []
    
    min_alleles_pred_col = []
    cnv_class_col = []
    cnv_subclass_col = []
    lq_ind_col = []
    all_lq_col = []
    lq_samps_col = []
    num_lq_col = []
    num_non_lq_col = []
    non_lq_samps_col  = []
    bi_class_col = []
    non_ref_af_col = []
    
    unique_gts_males_col = []
    unique_gts_females_col = []
    
    maf_data = []
    
    maf_col = []
    minor_allele_col = []
    
    try:
        males = [i for i in uuids if sex_dict[i]=='M']
    except:
        males = []
    
    try:
        females = [i for i in uuids if sex_dict[i]=='F']
    except:
        females = []
    
 
    
    
    for c in tran.columns:
        chrom = c.split('_')[1]
        
        try:
            l = tran[c].tolist()
            l_dict = tran[c].to_dict()
            
        except:
            print c
            
        
        lq_inds = []
        lq_samps = []
        if lq_dict:
            lq = lq_dict[c]
            if lq:
                # filter to those in the current sample subset
                lq_samps = [i for i in lq if i in uuids]

                if len(lq_samps) > 0:
                    uuids_to_test = set(uuids).difference(lq_samps)
                    
                else:
                    uuids_to_test = copy.deepcopy(uuids)

            else:
                uuids_to_test = copy.deepcopy(uuids)
        else:
            uuids_to_test = copy.deepcopy(uuids)
        
        l_all = [l_dict[s] for s in uuids]
        l_males_all = [l_dict[s] for s in males]
        l_females_all = [l_dict[s] for s in females]
       
        l_filt = [l_dict[s] for s in uuids_to_test]
        
      
        non_lq_samps_col.append(uuids_to_test)
        lq_samps_col.append(lq_samps)
        
        lq_len = len(lq_samps)
        num_lq_col.append(lq_len)
        
        num_pass = len(l_filt)
        num_non_lq_col.append(num_pass)
          
        # check if all are filtered- to continue processing- 
        # we will continue instead as if they aren't LQ, and adjust later using the all_lq column
        
        all_lq = False
        if len(l_filt) == 0:
            all_lq = True
            m = int(mode(l_all).mode[0])
        else:
            try:
                m = int(mode(l_filt).mode[0])
            except:
                print 'THIS', l_all, c
                break
    
        # sex specific modes      
        no_males = False
        no_females = False
        
        male_gts_filt = copy.deepcopy(l_males_all)
        if len(l_males_all) > 0:
            males_filt = list(set(males).intersection(uuids_to_test))
            if len(males_filt) > 0:
                male_gts_filt = [l_dict[i] for i in males if i in uuids_to_test]
                num_passing_males = len(male_gts_filt)
                mode_males = int(mode(male_gts_filt).mode[0])
                
            else:
            
                mode_males = int(mode(l_males_all).mode[0])
                num_passing_males = 0
                no_males = True
        else:
            num_passing_males = 0
            no_males = True
        
        female_gts_filt = copy.deepcopy(l_females_all)
        if len(l_females_all) > 0:
            females_filt = list(set(females).intersection(uuids_to_test))
            if len(females_filt) > 0:
                female_gts_filt = [l_dict[i] for i in females if i in uuids_to_test]
                num_passing_females = len(female_gts_filt)
                mode_females = int(mode(female_gts_filt).mode[0])
                
            else:
                no_females = True
                num_passing_females = 0
                mode_females = int(mode(l_females_all).mode[0])
                
        else:
            num_passing_females = 0
            no_females = True
        
        
        nnref = 0
        nnmode = 0
        count = 0
        nref_uuids = []
        alleles_dist = {}
        unique_nref_vars = []
        unique_non_mode_vars = []
        alleles_discovered = []
        
        alleles_males = []
        alleles_females = []
        
        #iterate through gts at each site
#         l = map(int, l) 

        s_to_test = set(uuids_to_test)
        
        all_vars = []
        for s in uuids:
            t = int(l_dict[s])
            
            if not all_lq:
                if s not in s_to_test:
                    continue
                    
            
            all_vars.append(t)
            sex = sex_dict[s]
            
            
            if sex == 'M':
                alleles_males.append(t)
            elif sex == 'F':
                alleles_females.append(t)
            
            # sex chroms
            if chrom == 'Y':
                if sex == 'M':
                    
                    alleles_discovered.append(t)
                    # record alleles only of males for y chrom
                    alleles_dist[t] = alleles_dist.get(t, 0) + 1


                    if t <> 1:
                        nnref += 1
                        unique_nref_vars.append(t)

                    if t <> mode_males:
                        nnmode +=1
                        unique_non_mode_vars.append(t)
                        nref_uuids.append(s)
            
            elif chrom == 'X':
                alleles_discovered.append(t)
                alleles_dist[t] = alleles_dist.get(t, 0) +1
                if sex == 'M':
                    
                    if t <> 1:
                        nnref += 1
                        unique_nref_vars.append(t)  
                    if t <> mode_males:
                        nnmode +=1
                        unique_non_mode_vars.append(t)
                        nref_uuids.append(s)
                else:
                        
                    if t <> 2:
                        nnref +=1
                        unique_nref_vars.append(t)
                    
                    if t <> mode_females:
                        nnmode +=1
                        unique_non_mode_vars.append(t)
                        nref_uuids.append(s)
            # autosomes      
            else:
                alleles_discovered.append(t)
                alleles_dist[t] = alleles_dist.get(t, 0) +1
                # lets also record vars that deviate from 2
                if t <> 2:
                    nnref += 1
                    unique_nref_vars.append(t)
                        
                # vars that deviate from the mode 
                if t <> m:
                    nnmode +=1
                    unique_non_mode_vars.append(t)
                    nref_uuids.append(s)
            
            
            count +=1
        
        alleles_dist_col.append(alleles_dist)
        alleles_dist_males_col.append(dict(Counter(alleles_males)))
        alleles_dist_females_col.append(dict(Counter(alleles_females)))
        
        
        
        # set of variants that deviates from expected diploid CN of reference (2 copies or 1 copy on male sex chroms)
        vars_at_site = list(set(unique_non_mode_vars))
        variants_discovered_col.append(vars_at_site)
       
    
        # number of different alleles taht deviate from expected diploid cn (not MODE cn)
        variant_allele_count_col.append(len(vars_at_site))
        
        count_alleles_col.append(len(alleles_dist))
        
        nnrefs_col.append(nnref)
        all_lq_col.append(all_lq)
        # calculate percent diff from mode CN - do this separately for x and y chroms
        if chrom == 'Y':
            percent_diff_from_mode = nnmode/(len(male_gts_filt))
            nonref_af = nnref/(len(male_gts_filt))
            
        else:
        
            if all_lq:
                num_pass = len(l_all)                
                
            else:
                num_pass = len(l_filt)
            
            percent_diff_from_mode = nnmode/num_pass
            nonref_af = nnref/num_pass
        
        percent_diff_from_mode_col.append(percent_diff_from_mode)
        non_ref_af_col.append(nonref_af)
        # col of uuids that are diff from mode
        
        nref_uuids_col.append(nref_uuids)
        modes_col.append(m)
        modes_male_col.append(mode_males)
        modes_female_col.append(mode_females)
        
        
        # col of unique genotypes at the site (all copy numbers)
        unique_gts_site = list(set(alleles_discovered))
        unique_gts_males = list(set(alleles_males))
        unique_gts_females = list(set(alleles_females))
        
         
        min_allele = min(unique_gts_site)
        max_allele  = max(unique_gts_site)
        
        
        gt_list = [unique_gts_males, unique_gts_females]
        min_max = []
        
        for g in gt_list:
            try:
                min_a = min(g)
            except:
                min_a = False
            try:
                max_a = max(g)
            except:
                max_a = False
            
            min_max.append([min_a, max_a])
            
        min_allele_males, max_allele_males = min_max[0]
        min_allele_females, max_allele_females = min_max[1]
          
        
#         min_allele_females = min(unique_gts_females)
#         max_allele_females = max(unique_gts_females)
        
        unique_gts_col.append(unique_gts_site)
        unique_gts_males_col.append(unique_gts_males)
        unique_gts_females_col.append(unique_gts_females)
        
      
        min_allele_col.append(min_allele)
        max_allele_col.append(max_allele)
        
        min_allele_males_col.append(min_allele_males)
        max_allele_males_col.append(max_allele_males)
        
        min_allele_females_col.append(min_allele_females)
        max_allele_females_col.append(max_allele_females)
        
        variant_class, cnv_class, min_pred_alleles = get_variant_type_info(min_allele, max_allele,
                                                                           min_allele_males,
                                                                       max_allele_males, min_allele_females, 
                                                                       max_allele_females, chrom)

        min_alleles_pred_col.append(min_pred_alleles)
        is_mixed_col.append(variant_class)
    
        
        cnv_subclass = variant_class + '_' + cnv_class
        cnv_subclass_col.append(cnv_subclass)
        
        # this will mark some sex chrom things as UNKNOWN we will revisit them
        bi_class = check_biallelic_status(unique_gts_site, unique_gts_males, 
                                                   unique_gts_females, cnv_class, chrom)
        bi_class_col.append(bi_class)
        
        maf = percent_diff_from_mode
        minor_allele = 'NA'
        if bi_class in ['DUP', 'DEL']:
            gt_dict = transform_to_biallelic(l_dict, uuids_to_test, chrom, sex_dict, 
                                             svtype = cnv_class)

            out_data, data_dict = compute_maf_lightweight(gt_dict, chrom, sex_dict)
            maf = data_dict['MAF']
            minor_allele = data_dict['Minor_Allele']
            
            out_data = out_data + [c, cnv_class]
            maf_data.append(out_data)
            
        maf_col.append(maf)
        minor_allele_col.append(minor_allele)
        
        cnv_class_col.append(cnv_class)
        
        
        
    out_cols = ['ID','cn_mode', 'diff_from_mode', 'percent_diff_from_mode', 'diploid_alleles',
                'diploid_alleles_males', 'diploid_alleles_females', 
                'count_alleles', 
               'alleles_dist', 'alleles_dist_males', 'alleles_dist_females',
                'variant_alleles', 'variant_allele_count', 
                'min_allele', 'max_allele', 'min_allele_males', 'max_allele_males', 'min_allele_females', 'max_allele_females', 'cnv_type','min_predicted_alleles', 'cnv_class', 'cnv_subclass', 'diff_mode_uuids', 'cn_mode_male', 'cn_mode_female', 'lq_samps', 'num_lq', 'num_pass', 'all_lq', 'pass_samps', 'num_pass_males', 'num_pass_females', 'bi_class', 'MAF', 'NONREF_AF']
    
    
    
    data_array = [tran.columns.tolist(), modes_col, nnrefs_col, percent_diff_from_mode_col, unique_gts_col,
                  unique_gts_males_col, unique_gts_females_col,
                  count_alleles_col,
                  alleles_dist_col, alleles_dist_males_col, alleles_dist_females_col,
                  variants_discovered_col, variant_allele_count_col,
                  min_allele_col,max_allele_col, min_allele_males_col, max_allele_males_col,
                  min_allele_females_col, max_allele_females_col,
                  is_mixed_col, min_alleles_pred_col,
                  cnv_class_col, cnv_subclass_col, nref_uuids_col, modes_male_col, modes_female_col, 
                  lq_samps_col, num_lq_col, num_non_lq_col, all_lq_col,  non_lq_samps_col, num_passing_males,
                  num_passing_females, bi_class_col, maf_col, non_ref_af_col]
    # return the result as a dict of lists (this can be directly dataframe converted or added as cols individually

        
    return  dict(zip(out_cols, data_array)), maf_data


# In[9]:


geno_all_with_stitch = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/i2QTL_combined/cns_all_with_stitch.pkl')

sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')

# subsets we will need for unifying genotype fields

samples_ipscore = sample_info[(sample_info.STUDY == 'iPSCORE')].WGS_ID.tolist()

samples_hipsci_fib = sample_info[(sample_info.STUDY == 'HipSci') & (sample_info.CELL_TYPE == 'Fibroblast')].WGS_ID.tolist()

info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/i2QTL_combined/gs_info_stitch_annot_final.pkl')

lq_dict = info.LQ_samps.to_dict()

sex_dict = sample_info.SEX.to_dict()

geno_all_with_stitch = geno_all_with_stitch.reindex(info.index)

tran = geno_all_with_stitch[samples_ipscore].T

t,t2 = gs_annotations(tran, samples_ipscore, sex_dict, lq_dict=lq_dict)

cols = ['NNREF', 'ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'ID', 'cnv_class']
maf_df = pd.DataFrame(t2, columns=cols)

maf_df = maf_df.set_index('ID')

cols = ['ID','cn_mode', 'diff_from_mode', 'percent_diff_from_mode', 'diploid_alleles',
                'diploid_alleles_males', 'diploid_alleles_females', 
                'count_alleles', 
               'alleles_dist', 'variant_alleles', 'variant_allele_count', 
                'min_allele', 'max_allele', 'min_allele_males', 'max_allele_males', 'min_allele_females', 'max_allele_females', 'cnv_type','min_predicted_alleles', 'cnv_class', 'cnv_subclass', 'diff_mode_uuids', 'cn_mode_male', 'cn_mode_female', 'lq_samps', 'num_lq', 'num_pass', 'all_lq', 'pass_samps', 'num_pass_males', 'num_pass_females', 'bi_class']

tdf = pd.DataFrame.from_dict(t)

tdf = tdf.set_index("ID", drop = False)
tdf = tdf[cols]


# In[10]:

def annotate_gencode_genes(cnv_info):
    
    
    tss_bt = pbt.BedTool('/publicdata/gencode_v19_20151104/tss_merged.bed')                     
    genes = pbt.BedTool('/publicdata/gencode_v19_20151104/genes.bed')
    exons = pbt.BedTool('/publicdata/gencode_v19_20151104/exons.bed')
    
    transcript_to_gene = '/publicdata/gencode_v19_20151104/transcript_to_gene.tsv'
    tg= pd.read_table(transcript_to_gene, index_col=0, header=None, squeeze=True)
    
    # need the chr prefix to intersect gencode 
    
    cnv_info.Chr = cnv_info.Chr.astype(str)
    cnv_info['chrom'] = cnv_info.Chr.apply(lambda x : 'chr' + x)
    
    
    cnv_bt = pbt.BedTool.from_dataframe(cnv_info[['chrom','Start', 'End', 'ID']]).sort()
    
    cnv_info.index = cnv_info.ID
    cnv_info.index.name = 'index'
    
    # Find genes that the CNV overlaps.
    res = cnv_bt.intersect(genes, sorted=True, wo=True)
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
    res = cnv_bt.intersect(exons, sorted=True, wo=True)
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


# In[11]:

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


# In[12]:

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
    


# In[13]:

def annotate_somatic_var(cnv_info, UUID, position, svtype):
    
    cnv_info['diff_from_mode_uuids_str'] = cnv_info.diff_mode_uuids.apply(lambda x: ",".join(x))
    chrom = position.split(':')[0]
    start = int(position.split(':')[1].split('-')[0])
    end = int(position.split(':')[1].split('-')[1])
    
    try:

        gdf = cnv_info.groupby('Chr').get_group(chrom)

        # remove all the ones in region that are singleton dels in that one person, excluding other sites

        gdf = gdf[(gdf.Start > start) & (gdf.End < end) & (gdf['diff_from_mode_uuids_str'] == UUID)]
        print 'somatic sites marked:'
        print gdf.shape

        ind = gdf[gdf.cnv_class == svtype].index.tolist()
        
        if len(ind) > 0:
            cnv_info['somatic'] = False
            cnv_info.loc[ind, 'somatic'] = True
    except:
        cnv_info['somatic'] = False
        
    return cnv_info


# In[16]:

def basic_length_annotations(cnv_info):
    
    # ensure that uuids list order is matching the order of the df before processing
    
    
    cnv_info['Length'] = cnv_info.End.astype(int) - cnv_info.Start.astype(int)
    
    size_bins = np.arange(1.6, 6, 0.2 )
    cnv_info['log_length']= np.log10(cnv_info.Length)
    cnv_info['log_size_bins']=pd.cut(cnv_info.log_length, size_bins, labels=size_bins[1:])
    
    return cnv_info


# In[ ]:

def annotate_info_LQ(df, fcns):
    """annotate the LQ counts for each site from LQ matrix onto info df,
    returns info df with new cols"""
    
    df = df.copy()
    
    g =  fcns.apply(pd.value_counts, axis = 1)
    cols = g.columns.tolist()
    
    # if all pass
    for i in  ['LQ', 'PASS']:
        if i not in cols:
            g[i] = 0
    
    
    g['LQ'] = g.LQ.fillna(0)
    g['PASS'] = g.PASS.fillna(0)
    
    df['LQ_NSAMP'] = g['LQ']
    df['PASS_NSAMP'] =  g['PASS']
    
    inds = fcns.index.tolist()
    tdf = fcns.T
    data = []
    for i in inds:
        t  = (tdf[i] == 'LQ')
        t = t[t]
        
        if t.shape[0] != 0:
            samps = t.index.tolist()
            data.append(samps)
        else:
            data.append(False)
    
    df['LQ_samps']  = data
    return df


# In[ ]:

def get_lq_union(df):
    df = df.copy()
    
    lq_cols = ['lq_samps_hipsci_fib', 'lq_samps_hipsci_ipsc', 'lq_samps_ipscore_fb']
    
    
    lq_cols2 = ['lq_samps_hipsci_fib', 'lq_samps_ipscore_fb']
    
    def union_of_cols(df, lq_cols):

        in_data = df[lq_cols].values
        out_unions = []
        num_lq  = []
        for i in in_data:
            s = set()
            for k in i:
                if type(k) == bool:
                    pass
                else:
                    s.update(set(k))
            out_unions.append(s)


            num_lq.append(len(s))
        return out_unions, num_lq
    

    out_unions, num_lq = df.pipe(union_of_cols, lq_cols)
    df['lq_union_all']  = out_unions
    df['num_union_all'] = num_lq
    
    df['lq_union_all_str']  = [",".join(list(i)) for i in out_unions]
    
    
    out_unions, num_lq = df.pipe(union_of_cols, lq_cols2)
    df['lq_union_bf']  = out_unions
    df['num_union_bf'] = num_lq
    
    df['lq_union_bf_str']  = [",".join(list(i)) for i in out_unions]
    return df


# In[20]:

def annotate_specific_set_of_samples(cnv_info, cns, uuids, sex_dict, suffix = 'FALSE', lq_dict = False):
    
    cns = cns[uuids]
    females = [i for i in uuids if sex_dict[i] == 'F']
    males = [i for i in uuids if sex_dict[i] == 'M']
    
     # transpose the df to increase speed (accessing columns is faster than rows
    cns_trans = cns.T
    data_summary, maf_data = gs_annotations(cns_trans, uuids, sex_dict, lq_dict=lq_dict)
    
    
    out_cols = ['ID','cn_mode', 'cn_mode_male', 'cn_mode_female', 'diff_from_mode', 'percent_diff_from_mode',
                'diploid_alleles', 'count_alleles', 'alleles_dist', 'alleles_dist_males',
                'alleles_dist_females','variant_alleles', 'variant_allele_count',
                'max_allele', 'min_allele', 'cnv_type','min_predicted_alleles', 'cnv_class', 
                'cnv_subclass', 'diff_mode_uuids', 'lq_samps', 'num_lq', 'num_pass', 'all_lq', 'MAF', 'NONREF_AF']
    
    
    if suffix == 'FALSE':
        suffix = ''
    
    
    cols = data_summary.keys()
    df_main = pd.DataFrame.from_dict(data_summary)
    df_main = df_main[out_cols]
    df_main = df_main.set_index('ID')
    
    col_names = [ c + suffix for c in df_main.columns]
    df_main.columns = col_names
    
    cols = ['NNREF', 'ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'ID', 'cnv_class']
    maf_df = pd.DataFrame(maf_data, columns=cols)
    maf_df = maf_df.set_index('ID')
    col_names = [ c + suffix for c in maf_df.columns]
    maf_df.columns = col_names
    
        
    return df_main, maf_df

    


# In[ ]:

def add_sample_set_annotations(cnv_info, cns, fcns, uuid_lists, suffix_lists, sex_dict, lq_adjust = False, lq_union = False, add_to_info = True):
    
    if lq_adjust:
        if 'LQ_samps' in cnv_info.columns.tolist():
            lq_dict = cnv_info['LQ_samps'].to_dict()
        else:
            cnv_info =  cnv_info.pipe(annotate_info_LQ, fcns)
            lq_dict = cnv_info['LQ_samps'].to_dict()
        
    else:
        lq_dict = False
    
    dfs_general = []
    dfs_maf = []
    
    for ul, sl in zip(uuid_lists, suffix_lists):
        print 'annotating of {}'.format(sl)
        gen_info, maf_info = annotate_specific_set_of_samples(cnv_info, cns, ul, sex_dict, sl, lq_dict=lq_dict)
        dfs_general.append(gen_info)
        dfs_maf.append(maf_info)
         
    
    to_add = pd.concat(dfs_general, axis = 1)
    maf_all = pd.concat(dfs_maf, axis = 1)
    
    if add_to_info:
        for c in to_add.columns:
            cnv_info[c] = to_add[c]
            
            
    if lq_union:
        try:
            cnv_info = cnv_info.pipe(get_lq_union)
        except:
            print 'could not annotate lq unions'
       
    
    return cnv_info, maf_all, to_add
    


# In[34]:


# geno_all_with_stitch = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/i2QTL_combined/cns_all_with_stitch.pkl')

# sample_info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/sample_info.pkl')

# # subsets we will need for unifying genotype fields

# samples_ipscore = sample_info[(sample_info.STUDY == 'iPSCORE')].WGS_ID.tolist()

# samples_hipsci_fib = sample_info[(sample_info.STUDY == 'HipSci') & (sample_info.CELL_TYPE == 'Fibroblast')].WGS_ID.tolist()

# info = pd.read_pickle('/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/gs_processing_V3/i2QTL_combined/gs_info_stitch_annot_final.pkl')

# # lq_dict = info.LQ_samps.to_dict()

# sex_dict = sample_info.SEX.to_dict()

# geno_all_with_stitch = geno_all_with_stitch.reindex(info.index)

# tran = geno_all_with_stitch[samples_ipscore].T

# t,t2 = gs_annotations(tran, samples_ipscore, sex_dict, lq_dict=lq_dict)

# cols = ['NNREF', 'ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'ID', 'cnv_class']
# maf_df = pd.DataFrame(t2, columns=cols)

# maf_df = maf_df.set_index('ID')

# cols = ['ID','cn_mode', 'diff_from_mode', 'percent_diff_from_mode', 'diploid_alleles',
#                 'diploid_alleles_males', 'diploid_alleles_females', 
#                 'count_alleles', 
#                'alleles_dist', 'variant_alleles', 'variant_allele_count', 
#                 'min_allele', 'max_allele', 'min_allele_males', 'max_allele_males', 'min_allele_females', 'max_allele_females', 'cnv_type','min_predicted_alleles', 'cnv_class', 'cnv_subclass', 'diff_mode_uuids', 'cn_mode_male', 'cn_mode_female', 'lq_samps', 'num_lq', 'num_pass', 'all_lq', 'pass_samps', 'num_pass_males', 'num_pass_females', 'bi_class']

# tdf = pd.DataFrame.from_dict(t)

# tdf = tdf.set_index("ID", drop = False)
# tdf = tdf[cols]


# In[33]:

# sample_lists = [samples_ipscore, samples_hipsci_fib]
# suff_list = ['_ipscore_fb', '_hipsci_fib']


# info, dfs_maf, dfs_general = add_sample_set_annotations(info, geno_all_with_stitch, False, sample_lists, suff_list, sex_dict, lq_adjust=True, add_to_info = False)


# In[32]:

# dfs_general[dfs_general.MAF_ipscore_fb != dfs_general.percent_diff_from_mode_ipscore_fb][['alleles_dist_ipscore_fb', 'MAF_ipscore_fb', 'percent_diff_from_mode_ipscore_fb']].percent_diff_from_mode_ipscore_fb.hist()


# In[18]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-cns", "--cns", dest="cns", metavar='<cns>', help="genome strip cns pickle (Copy Number States)", required=True)
    
    parser.add_argument("-fcns", "--fcns", dest="fcns", metavar='<fcns>', help="genome strip fcns pickle (Filter Column for each sample/site)", required=True)
    
    parser.add_argument("-info", "--info", dest="info", metavar='<info>', help="genome strip info pickle", required=True)
    
    parser.add_argument("-s", "--samples", dest="samples", metavar='<fn_samples1,fn_samples2>', help="sets of samples to annotate", required=True)
    
    parser.add_argument("-gender_map", "--gender_map", dest="gender_map", metavar='<gender_map>', help="gender file UUID Sex, tab delimited", required=True)   
    
    parser.add_argument("-intersect", "--intersections", dest="intersect", metavar='<True/False>', help="intersect with MHC VDJ Centromeres, Telomeres", required=False, default=True)
    
    parser.add_argument("-maf_only", "--maf_only", dest="maf_only", action = 'store_true', help="only annotate maf from a subset of samples, skip other annotations")
    
    
    parser.add_argument("-LQ_adjust", "--lq_adjust", dest="lq_adjust", action = 'store_true', help="adjust for low qual samples during MAF calculations")
    
    
    parser.add_argument("-LQ_union", "--lq_union", dest="lq_union", action = 'store_true', help="annotate the union of lq samples per site, num LQ")
    
    
    parser.add_argument("-ss_suff", "--suffix_subsets", dest="suffix_subsets", metavar='<suff1,suff2,suff3>', help="suffixes to use for naming columns of annotated subsets of samples ex: unrelated, related", required=True, default=False)
    
    
    parser.add_argument("-somatic", "--somatic", dest="somatic", metavar='<Sample,region,svtype>', help="annotate presence of a known somatic variant, svtype indicates which variant type the somatic variant is.  useful if you know that you have a somatic variant with other evidence such as arrays from different cell types on the same individuals", required=False, default=False)
    
    
    parser.add_argument("-genes", "--genes", dest="genes", action = 'store_true', help="annotate intersection with gencode genes")
    
    
    parser.add_argument("-conv_int", "--conv_int", dest="convert_int", action = 'store_true', help="force data in cns columns to be integer dtype")
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    parser.add_argument("-pre", "--prefix", dest="prefix", metavar='<prefix>', help="prefix to name files", default = False)
    
    parser.add_argument("-file_suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to annotate extracted info and genotypes from genome_strip_extract utility with various things such as gene intersections, minor allele frequency etc.')
    add_arguments_to_parser(parser)
    return parser


# In[23]:

def run_from_args(args):
    gs_cns = pd.read_pickle(args.cns)
    gs_info = pd.read_pickle(args.info)
    gs_fcns = pd.read_pickle(args.fcns)
    
    # ensure ordering is the same from the start
    gs_cns = gs_cns.reindex(gs_info.index)
    gs_fcns = gs_fcns.reindex(gs_info.index)
    
    
    sample_files = args.samples
    sample_files = sample_files.split(',')
    
    sample_lists = []
    for fn in sample_files:
        samples = [line.rstrip() for line in open(fn)]
        sample_lists.append(samples)
        
    suffixes = args.suffix_subsets
    suffixes = suffixes.split(',')
    
    
    
    
    gender_file = args.gender_map
    
    sex_dict = {line.rstrip().split()[0]:line.rstrip().split()[1] for line in open(gender_file)}
   
    
    print 'starting annotation'
    print CM.datestring(hour=True, minute=True)
    
    if args.convert_int:
        for s in samples:
            gs_cns[s] = gs_cns[s].astype(int)



        
    output_location = args.output_dir
    
    if args.maf_only:
        print "annotating subset(s)"
        
        gs_info, maf_df, to_add = add_sample_set_annotations(gs_info, gs_cns, gs_fcns, sample_lists, 
                                                             suffixes, sex_dict,lq_adjust=args.lq_adjust, 
                                                             lq_union = args.lq_union)
    
    else:
        print 'annotating data subsets'
        
        gs_info, maf_df, to_add = add_sample_set_annotations(gs_info, gs_cns, gs_fcns, sample_lists, 
                                                             suffixes, sex_dict,
                                                             lq_adjust=args.lq_adjust, lq_union = args.lq_union)
        
        gs_info = basic_length_annotations(gs_info)
        
        if args.intersect:
            print "annotating centromere/telomere distances, MHC, VDJ Regions"
            gs_info = annotate_filters(gs_info)
    #         print gs_info

        if args.somatic:
            print "annotating somatic variants"
            spl = args.somatic.split(',')
            uuid = spl[0]
            region = spl[1]
            svtype = spl[2]
            gs_info = annotate_somatic_var(gs_info, uuid, region, svtype)

        if args.genes == True:
            print "annotating gencode genes"
            gs_info = annotate_gencode_genes(gs_info)




    
    if args.suffix:
        fn_info = os.path.join(output_location, 'gs_info' + args.suffix)
        fn_cns = os.path.join(output_location, 'gs_cns' + args.suffix)
        
        var_name_info = 'gs_info' + args.suffix
        var_name_cns  = 'gs_cns' + args.suffix
        var_name_maf = 'gs_maf_bi' + args.suffix
    
    else:
        fn_info = os.path.join(output_location, 'gs_info')
        fn_cns = os.path.join(output_location, 'gs_cns')
        var_name_info = 'gs_info'
        var_name_cns  = 'gs_cns'
        var_name_maf = 'gs_maf_bi'
        
    
    print 'data annotated'
    # is this necessary- probably don't need to save the cns frame again
#     if args.cns_reset_ind:
#         CM.save_dataframe(var_name_cns, gs_cns, output_location, print_vars_recorded_loc=False, reset_index = True, index = False)
    
#     else:
#         CM.save_dataframe(var_name_cns, gs_cns, output_location, print_vars_recorded_loc=False)
    
    CM.save_dataframe(var_name_info, gs_info, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_maf, maf_df, output_location, print_vars_recorded_loc=False)


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

