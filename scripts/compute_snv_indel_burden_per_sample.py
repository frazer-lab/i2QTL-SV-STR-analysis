#!/usr/bin/env python
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
import scipy.stats as stats

import argparse
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import  make_axes_locatable
import datetime
import time


# In[2]:


def find_header_end(fn):
    """ find header end line number of vcf or gzipped vcf"""

    if fn.split('.').pop() == 'gz':
        import gzip
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')


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


# In[3]:



def parse_info_col(t, lab, type_out = str, alternative = 'None'):
    t = t.split(';')
    # filter out the tags with no keyword, if any
    t = [i for i in t if i.find('=') != -1]
    cols = [i.split('=')[0] for i in t]
    try:
        vals = [i.split('=')[1] for i in t]    
    except:
        return "PARSE_ERROR"
    try:
        ind = cols.index(lab)
        v = vals[ind]
        
        try:
            return type_out(v)
        except:
            return alternative
    except:
        return 'Column_Not_Present'


# In[4]:


def parse_format_col(t, col, lab, type_out = str, alternative = 'None', change_dict = False):
    """ use the format column to pull out info of a specific VCF record
    in a genotyping matrix (unparsed)"""
    
    f= t['FORMAT']
    c = t[col]
    c  = c.split(':')
    f = f.split(':')
    
    try:
        ind = f.index(lab)
        v = c[ind]
        try:
            if change_dict:
                try:
                    v = change_dict.get(v, v)
                except:
                    return 'ERROR'
            
            return type_out(v)
        except:
            return alternative
            
    except:
        return 'Column_Not_Present'
    


# In[5]:


def parse_format_mod(t, col, lab, except_out = 'None'):
    """ use the format column to pull out info of a specific VCF record
    in a genotyping matrix (unparsed)"""
    
    f= t['FORMAT']
    c = t[col]
    c  = c.split(':')
    f = f.split(':')
    
    try:
        ind = f.index(lab)
        v = c[ind]
        return v
    except:
        return except_out


# In[6]:


def coord_extract(df,chrom, start, end, contained = True):
    
    if contained:
        return df[(df.Chr==chrom) & (df.POS >= start) & (df.END <= end)]


# In[169]:


def get_geno_fields(labs, samples, lin_spl, header_dict, format_dict):
    gts_out = [{} for l in labs]
    missing_samps = []
    for s in samples:
        d = lin_spl[header_dict[s]].split(':')
        # cases where we don't have format info for all categories
        if len(d) == 1:
            if d[0] == '.':
                missing_samps.append(s)
        
        for i,l in enumerate(labs):
            if len(d) == 1:
                if d[0] == '.':
                    gt = './.'
                
            else:
                gt = d[format_dict[l]]
                if gt in ['.', './.']:
                    gt = './.'
                    missing_samps.append(s)
                
                
            
            gts_out[i][s] = gt
    missing_samps = list(set(missing_samps))
    return gts_out, missing_samps


# In[9]:


def safe_div(x, y, alt=0):
    try:
        return x/y
    except:
        return alt


# In[10]:


def get_max_min_dict(d, samples):
    v = [d[s] for s in samples]
    v = [float(i) for i in v if i != './.']
    max_v = round(max(v),4)
    min_v = round(min(v), 4)
    
    return max_v, min_v


# In[4]:


def add_counts_to_per_sample(d, samples, cumulative_dict):
    
    for s in samples:
        c = str(d[s])
        
        if c == ['0/0'] :
            cat = 'REF'
        
        elif c in ['./.', '.']:
            cat = 'MISSING'
        
        else:
            cat = 'NREF'
        
        cumulative_dict[s][cat] = cumulative_dict[s].get(cat, 0) + 1


# In[10]:


def add_counts_to_per_burden(d, samples, cumulative_dict, length, cat):  
    copies = {'0/0':0, '0/1': 1, '1/1':2, './.':0, '.':0, '1/0': 1}
    for s in samples:
        c = str(d[s])
        try:
            bp = length * copies[c]
        except:
            print set(d.values())
            return 'FAIL'
            break
        cumulative_dict[cat][s]['Burden'] = cumulative_dict[cat][s].get('Burden', 0) + bp


# In[127]:


def calculate_pair_concordance(pair_samples, data_dict):
         
    number_with_var = 0
    number_concordant_with_var = 0
    number_discordant_with_var = 0
    number_missing = 0
    per_pair_data = []
    for tp in pair_samples:

        vals =  map(str, [data_dict[tp[0]], data_dict[tp[1]]])

        # check one way- first twin- for variant
        if vals[0] == './.':
            number_missing +=1
            per_pair_data.append('MISSING')
        
        if vals[0] not in ['0/0', './.', '0']:

            number_with_var +=1
            if vals[0] == vals[1]:
#                 tps_concordant.append("_".join(tp))
                number_concordant_with_var +=1
                per_pair_data.append('CONCORDANT')
            else:
#                 tps_discordant.append("_".join(tp))
                number_discordant_with_var +=1 
                per_pair_data.append('DISCORDANT')
        elif vals[0] in ['0/0', '0']:
            per_pair_data.append('REF')

    if number_with_var > 0:
        # replication rate- ranges from 0-100% 100% means every variant is a match
        replication_rate = 1 - (number_discordant_with_var/number_with_var)
    else:
        replication_rate = 'None'
    
    if number_with_var > 0:
        return [per_pair_data, number_with_var, number_concordant_with_var, number_discordant_with_var,
                number_missing, replication_rate]
    else:
        return False


# In[116]:


def calculate_alt_allele_freq(gts_dict, samples, chrom, num_samples):
    """ calculate alt allele freq """
    
 
    #uuids containing the 1 allele
    non_ref = 0
    non_ref_uuids = []

    #uuids containing the 0 allele 
    ref = 0
    ref_uuids = []


    # for chroms that are diploid
    # correct for allelic probs on sex chroms

    missing = 0
    alleles_dist = {}
    
    for s in samples:
        c = str(gts_dict[s])
        
        alleles_dist[c] = alleles_dist.get(c, 0) + 1
        
        
        if c == './.':
            missing +=1
        
        if c not in ['0|0','./.', '0']:
            non_ref +=1
            non_ref_uuids.append(s)

        if c in ['0|0', '0']:
            ref +=1
            ref_uuids.append(s)


    percent_missing = round((missing/num_samples), 3)
    
    if (ref == 0) & (non_ref > 0): 
        NNREF_Freq = 1.0
          
    else:
        NNREF_Freq = round(safe_div(non_ref, ref, alt= 0), 4)
    

    
    out_names = ['NREF', 'REF', 'NMissing', 'NNREF_UUIDs','NNREF_AF', 'ALLELES_DIST', 'PERC_MISSING']
    
    
    alleles_dist = ",".join(["{}:{}".format(i,k) for i,k in alleles_dist.iteritems()])
    out_data = [non_ref, ref, missing, ",".join(non_ref_uuids), NNREF_Freq, alleles_dist, percent_missing]
    data_dict = dict(zip(out_names, out_data))
    return data_dict, out_names, out_data


# In[168]:


def prep_per_sample(df):
    df = df.copy()
    df = df.T
    return df


# In[3]:


def classify_snv_indel(REF, ALT):
    
    len_ref = len(REF)
    len_alt = len(ALT)
    
    variant_type = 'SNV'
    
    if ALT == '*':
        variant_type = 'DEL'
        bp_change = len_ref
        
    elif len_ref == len_alt:
        bp_change = 0
    
    elif len_alt > len_ref:
        variant_type = 'INS'
        bp_change = len_alt - len_ref
    
    elif len_ref > len_alt:
        
        variant_type = 'DEL'
        bp_change = len_ref - len_alt
        
    else:
        pass
    
    return variant_type, bp_change


# In[4]:


def calculate_alt_allele_freq(gts_dict, samples, chrom):
    """ calculate alt allele freq """
    
 
    #uuids containing the 1 allele
    non_ref = 0
    non_ref_uuids = []

    #uuids containing the 0 allele 
    ref = 0
    ref_uuids = []
    
    ref_allele = 0
    alt_allele = 0


    # for chroms that are diploid
    # correct for allelic probs on sex chroms

    missing = 0
    alleles_dist = {}
    
    allele_dict = {'0/0':0, '0/1':1, '1/1':2, './.':0, '1/0':1}
    allele_dict_sex = {'0/0':0, '0/1':1, '1/1':1, './.':0}
    
    # most common
    s = set(samples)
    gts_filt = {k:v for k,v in gts_dict.iteritems() if k in s}
    vals = gts_filt.values()
    C = Counter(vals)
    mc = C.most_common()
    try:
        mc_allele, count_mc_allele = mc[0]
    except:
        pass
    
    
    for s in samples:
        c = str(gts_dict[s])
        
        
        alleles_dist[c] = alleles_dist.get(c, 0) + 1
        aa = allele_dict[c]
        ra = 2 - aa
        
        if c == './.':
            missing +=1
        
        if c not in ['0|0', '0/0', './.', '0']:
            non_ref +=1
            non_ref_uuids.append(s)

        if c in ['0|0', '0', '0/0']:
            ref +=1
            ref_uuids.append(s)
        
        
        ref_allele += ra
        alt_allele += aa
        
    
    percent_missing = round((missing/len(samples)), 3)
    
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
    
    if (ref == 0) & (non_ref > 0): 
        NNREF_Freq = 1.0
          
    else:
        NNREF_Freq = round(safe_div(non_ref, ref + non_ref, alt= 0), 4)
        
    
    
    major_allele = {'ALT':'REF', 'REF':'ALT'}[minor_allele]
    out_names = ['NREF', 'REF', 'NMissing', 'NNREF_UUIDs','NNREF_AF', 'ALLELES_DIST', 'PERC_MISSING', 'MAF', 
                 'MINOR_ALLELE', 'MAJOR_ALLELE', 'MODE_ALLELE']
    alleles_dist = ",".join(["{}:{}".format(i,k) for i,k in alleles_dist.iteritems()])
    out_data = [non_ref, ref, missing, ",".join(non_ref_uuids), NNREF_Freq, alleles_dist, percent_missing, maf,
                minor_allele, major_allele, mc_allele]
    data_dict = dict(zip(out_names, out_data))
    return data_dict


# In[159]:


def prep_per_pair(df):
    df = df.copy()
    df = df.T
    df['RR'] = df['CONCORDANT'] / (df['CONCORDANT'] + df['DISCORDANT'])
    return df


# In[1]:


from collections import Counter


# In[5]:


def process_vcf_and_generate_qc_info(fn, out_dir, ipscore_samples, pair_samples, suff = 'auto', chroms = False,
                                    singleton_only = False, maf_thresh = False):
    
    
    count = 0
    progress_level = 0
    progress_increment = 100000
    
    flattened_pairs = list(set([i for sublist in pair_samples for i in sublist]))
    
    num_samples = len(ipscore_samples)
    if fn.split('.').pop() == 'gz':
        if not chroms:
            chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
        
        command = "bcftools view -r {} {}".format(",".join(chroms), fn)
        F = subprocess.Popen(command, shell=True, stdout= subprocess.PIPE).stdout
        
        
    else:
        F = open(fn, 'rU')
    

    replication_snv_file = open(os.path.join(out_dir, 'snv_replication_gt_info.{}.tsv'.format(suff)), 'w')
    replication_indel_file = open(os.path.join(out_dir, 'indel_replication_length_info.{}.tsv'.format(suff)), 'w')
    indel_info_file = open(os.path.join(out_dir, 'indel_length_info.{}.tsv'.format(suff)), 'w')
   
    
  
    per_pair_data_cumulative_snv_gt = {"_".join(tp):{'CONCORDANT':0, 'DISCORDANT':0, 'REF':0, 'MISSING':0} for tp in pair_samples}
    per_pair_data_cumulative_ins_gt = {"_".join(tp):{'CONCORDANT':0, 'DISCORDANT':0, 'REF':0, 'MISSING':0} for tp in pair_samples}
    
    per_pair_data_cumulative_del_gt = {"_".join(tp):{'CONCORDANT':0, 'DISCORDANT':0, 'REF':0, 'MISSING':0} for tp in pair_samples}
    
    
    per_pair_dicts = {'SNV': per_pair_data_cumulative_snv_gt, 'DEL': per_pair_data_cumulative_del_gt,
                      'INS': per_pair_data_cumulative_ins_gt}
    
    

    overall_count_dict = {}
    
    count_missing = 0
    num_variants = 0
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Starting Variant Processing: {}".format(ts)
    
    h = True
    for line in F:

        line = line.rstrip()
        lin_spl = line.split()
        if h:
            count +=1       
            if line.find('#CHROM') == 0:
                h = False
                print 'Encountered Header End'
                header_count = copy.deepcopy(count)
                header = copy.deepcopy(lin_spl)
                header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
    #             print header_dict

                samples = header[9:]
                info_cols = header[:9]

                per_sample_data_snv_cumulative = {s:{} for s in samples}
                per_sample_data_del_cumulative = {s:{} for s in samples}
                per_sample_data_ins_cumulative = {s:{} for s in samples}
                per_sample_dicts = {'SNV': per_sample_data_snv_cumulative, 'DEL': per_sample_data_del_cumulative,
                          'INS': per_sample_data_ins_cumulative}
                
                per_sample_dicts_noxy = copy.deepcopy(per_sample_dicts)
                per_sample_burden_dicts_noxy = copy.deepcopy(per_sample_dicts)


                # may not use these
                cols_af =  ['NREF', 'REF', 'NMissing', 'NNREF_UUIDs','NNREF_AF', 'ALLELES_DIST', 'FRAC_MISSING']
                cols_af_i2QTL_unrel = ["{}_i2QTL_unrel".format(i) for i in cols_af]    

                cols_replication = ['NP_VAR', 'NP_CONC', 'NP_DISC', 'NP_MISSING', 'RR']
                indel_info_header = ['CHROM', 'POS', 'ID', 'VARIANT_TYPE', 'BP_CHANGE']
                


                replication_file_header = (['ID'] + cols_replication + ['variant_type', 'bp_change'])
                replication_indel_file.write("\t".join(replication_file_header) + '\n')
                replication_snv_file.write("\t".join(replication_file_header) + '\n')
                indel_info_file.write("\t".join(indel_info_header) + '\n')
                
            else:
                continue


        else:
      
            if progress_level == progress_increment:
                
                d = datetime.datetime.now()
                ts = d.strftime('%D- %H:%M')
                print "processed {} variants {}".format(num_variants, ts)
                progress_level = 0

            progress_level +=1
            
            count +=1

            format_fields = lin_spl[header_dict['FORMAT']].split(':')
            format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}

            num_variants +=1

            info_col = header_dict['INFO']

            POS = lin_spl[header_dict['POS']]
            ID = lin_spl[header_dict['ID']]
            info = lin_spl[header_dict['INFO']]
            chrom = str(lin_spl[header_dict['#CHROM']])
            REF = lin_spl[header_dict['REF']]
            ALT = lin_spl[header_dict['ALT']]
            variant_type, bp_change = classify_snv_indel(REF, ALT)
            
            if variant_type != 'SNV':
                out = [chrom, POS, ID, variant_type, bp_change]
                indel_info_file.write("\t".join(map(str, out)) + '\n')
            
            overall_count_dict[variant_type] = overall_count_dict.get(variant_type, 0) + 1
            
            
            fields_desired = ['GT']
            geno_fields_dict = {fields_desired[i]:i for i in range(0, len(fields_desired))}

            geno_extracted, missing_samps = get_geno_fields(fields_desired, samples, lin_spl, 
                                                            header_dict, format_dict)
            gt_dict = geno_extracted[geno_fields_dict['GT']]
            
            singleton = False
            samps_to_use = set(ipscore_samples).difference(missing_samps)
            alleles_count = Counter([gt_dict[i] for i in samps_to_use])
            mc = alleles_count.most_common()
            if len(mc) == 2:
                singleton = (mc[-1][1] == 1)
            
            
            subset_missing = [i for i in missing_samps if i in ipscore_samples]
            num_missing = len(subset_missing)
            num_passing = num_samples - num_missing
            
            if singleton_only:
                if not singleton:
                    continue

            # ======= Add info to cumulative dict ==========
            bp = 1
            if variant_type != 'SNV':
                bp = bp_change
                
            add_counts_to_per_sample(gt_dict, samples,  per_sample_dicts[variant_type])


                #================TWIN/PAIR Replication Calculations ==========
                
                # check for non_missing_gt
            if chrom not in ['X', 'Y']:
                if not singleton_only:
#                     data_dict = calculate_alt_allele_freq(gt_dict, ipscore_samples, chrom)
#                     out_line = [ID, chrom, POS, caller, data_dict['NNREF_AF'], 
#                                 data_dict['ALLELES_DIST'], data_dict['NREF'],
#                                 svtype, FILTER, data_dict['MAF'], data_dict['MODE_ALLELE']]
                    
                    
                    if maf_thresh:
                        if singleton:
                            continue
                        
                        cols_info_extracted = ['MAF_UNREL']
                        cols_dict = dict(zip(cols_info_extracted, range(0, len(cols_info_extracted))))
                        info_out = []
                        info_cols_dict = {}
                    
                        for l in cols_info_extracted:
                            c = parse_info_col(info, l)
                            info_out.append(c)
                            info_cols_dict[l] = c
                        try:
                            MAF = float(info_cols_dict['MAF_UNREL'])
                        except:
                            continue

                        if ((singleton == False) & (MAF > maf_thresh[0]) & (MAF <= maf_thresh[1])):
                            pass
                        else:
                            continue
               
                TEST = add_counts_to_per_burden(gt_dict, samples,  per_sample_burden_dicts_noxy, bp, variant_type)
                if TEST == 'FAIL':
                    print REF, ALT, ID, chrom, POS
                    print lin_spl
                    break
                    return
        
                
                
                add_counts_to_per_sample(gt_dict, samples, per_sample_dicts_noxy[variant_type])
                
                if singleton:
                    # don't collect any information if its just singletons
                    continue
                t = [gt_dict[s] for s in flattened_pairs if gt_dict[s] != './.']
                if len(t) > 0:
                    pair_data_gt = calculate_pair_concordance(pair_samples, gt_dict)
                    if pair_data_gt != False:

                        (per_pair_data_gt, number_with_var, 
                         number_concordant_with_var, number_discordant_with_var,
                         number_missing, replication_rate) = pair_data_gt

                        OUT_PAIR_DATA_GT = copy.deepcopy([number_with_var, number_concordant_with_var,
                                         number_discordant_with_var, number_missing, replication_rate])
                        
                        
                        
                        for tp, d_gt in zip(pair_samples, per_pair_data_gt):
                            p = "_".join(tp)
                            per_pair_dicts[variant_type][p][d_gt] = (per_pair_dicts[variant_type][p].get(d_gt, 0) + 1)
                        
                        OUT_LINE_REPLICATION_GT = ([ID] + OUT_PAIR_DATA_GT  + [variant_type, bp_change])
                        
                        name = "{}_seg".format(variant_type)
                        overall_count_dict[name] = overall_count_dict.get(name, 0) + 1 
                        
                        if variant_type == 'SNV':                    
                            replication_snv_file.write("\t".join(map(str,OUT_LINE_REPLICATION_GT)) + '\n') 
                        else:
                            replication_indel_file.write("\t".join(map(str,OUT_LINE_REPLICATION_GT)) + '\n')

        
    
    replication_snv_file.close()
    replication_indel_file.close()
    print "saving files"
    cats = ['SNV', 'INS', 'DEL']
    for c in cats:
        if not singleton_only:
            
            ppd = per_pair_dicts[c]
            fn_per_pair_gt = os.path.join(out_dir, '{}_per_pair_gt.{}.tsv'.format(c, suff))
            per_pair_gt  = pd.DataFrame(ppd).pipe(prep_per_pair)
            per_pair_gt.to_csv(fn_per_pair_gt, sep = '\t')
            print fn_per_pair_gt

        fn_per_sample = os.path.join(out_dir, '{}_per_sample_gt.{}.tsv'.format(c, suff))
        psd = per_sample_dicts[c] 
        per_sample_gt = pd.DataFrame(psd).pipe(prep_per_sample)
        per_sample_gt.to_csv(fn_per_sample, sep = '\t')
        print fn_per_sample
        
        fn_per_sample = os.path.join(out_dir, '{}_per_sample_burden_auto.{}.tsv'.format(c,suff))
        psd = per_sample_burden_dicts_noxy[c] 
        per_sample_gt = pd.DataFrame(psd).pipe(prep_per_sample)
        per_sample_gt.to_csv(fn_per_sample, sep = '\t')
        print fn_per_sample

        fn_per_sample_noxy = os.path.join(out_dir, '{}_per_sample_auto_gt.{}.tsv'.format(c, suff))
        psd = per_sample_dicts_noxy[c] 
        per_sample_gt = pd.DataFrame(psd).pipe(prep_per_sample)
        per_sample_gt.to_csv(fn_per_sample_noxy, sep = '\t')
        print fn_per_sample

    fn_overall_count =  os.path.join(out_dir, 'variant_counts.{}.tsv'.format(suff))
    tdf = pd.Series(overall_count_dict).to_frame('Num_Variants')
    tdf['VARIANT_TYPE'] = tdf.index
    tdf.to_csv(fn_overall_count, sep = '\t', index =False, header = False)

    print "number of variants: {}".format(num_variants)
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Completed Variant Processing: {}".format(ts)
    print 'COMPLETE'
    
    return
                                


# In[111]:


def add_arguments_to_parser(parser):
    

    parser.add_argument("-vcf", "--vcf", dest="vcf_file", metavar='<vcf_file>', help="vcf file from lumpy/speedseq pipeline, may be gzipped or not", required=True)
    
    parser.add_argument("-pairs", "--pairs", dest="pairs_fn", metavar='<pair_fn>', help="file of pairs of samples, likely genetically duplicate samples (twins, replicates, biological replicates) (tsv with no header)", required=True)
    
    parser.add_argument("-s", "--samples", dest="samples_fn", metavar='<samples_fn>', help="list of samples on which to compute MAF and other general stats", required=True)
    
    
    parser.add_argument("-suff", "--suff", dest="suffix", metavar='<suffix>', help="suffix of info columns", required=False, default = False)
        
    parser.add_argument("-chroms", "--chroms", dest="chroms", metavar='<chroms>', help="chroms to process", required=False, default = False)
    
    parser.add_argument("-maf_min_max", "--maf_min_max", dest="maf_min_max", metavar='<maf_min_max>', help="minor allele frequency min and max, comma separated, to use", required=False, default = False)
  
    parser.add_argument("-singleton", "--singleton", dest="singleton", action = 'store_true', help="compute burdens/etc for singletons only")
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:


def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to extract info and calculate replication rate for HipSTR vcf files')
    add_arguments_to_parser(parser)
    return parser



# In[1]:


def run_from_args(args):
    vcf_fn = args.vcf_file
    samples_fn = args.samples_fn
    suff = args.suffix
    pair_fn = args.pairs_fn
    out_dir = args.output_dir
    pairs = [line.rstrip().split() for line in open(pair_fn)]
    samples = [line.rstrip() for line in open(samples_fn)]   
    chroms = args.chroms
    maf_min_max = args.maf_min_max
    if chroms:
        chroms = chroms.split(",")
    if maf_min_max:
        maf_min_max = map(float, maf_min_max.split(','))
        
    process_vcf_and_generate_qc_info(vcf_fn, out_dir, samples, pairs, suff = suff, chroms = chroms, 
                                     maf_thresh = maf_min_max, singleton_only = args.singleton)


# In[87]:


if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))


# In[165]:


# fn = '/frazer01/projects/hipsci/pipeline/WGS/HipSTR/combined_results/hipstr_ipscore.vcf'

# out_dir = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/hipstr_qc_analysis/test'



# twin_fn = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/ipscore_sample_info/twins.tsv'
# pairs = [line.rstrip().split() for line in open(twin_fn)]    

# fn_ipscore = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/sample_info_combined/samples_ipscore.txt'

# samples = [line.rstrip() for line in open(fn_ipscore)]

# p,d = process_vcf_and_generate_qc_info(fn, out_dir, samples, pairs)

# fn = '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/hipstr_qc_analysis/ipscore_geno/hipstr_per_sample_length.tsv'

# fn =  '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/hipstr_qc_analysis/test/hipstr_per_sample_gt.tsv'

# fn =  '/frazer01/projects/hipsci/analysis/i2QTL-sv-analysis/private_output/hipstr_qc_analysis/test/hipstr_per_sample_length.tsv'

# t = pd.read_table(fn)

