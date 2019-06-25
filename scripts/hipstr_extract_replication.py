
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


# In[170]:

def add_counts_to_per_sample(d, samples, cumulative_dict):
    
    for s in samples:
        c = str(d[s])
        
        if c== '0':
            cat = 'REF'
        
        elif c in ['./.', '.']:
            cat = 'MISSING'
        
        else:
            cat = 'NREF'
        
        cumulative_dict[s][cat] = cumulative_dict[s].get(cat, 0) + 1


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
        
        if vals[0] not in ['0|0', './.', '0']:

            number_with_var +=1
            if vals[0] == vals[1]:
#                 tps_concordant.append("_".join(tp))
                number_concordant_with_var +=1
                per_pair_data.append('CONCORDANT')
            else:
#                 tps_discordant.append("_".join(tp))
                number_discordant_with_var +=1 
                per_pair_data.append('DISCORDANT')
        elif vals[0] in ['0|0', '0']:
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
        
#         sex = sex_dict[s]
        
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


# In[159]:

def prep_per_pair(df):
    df = df.copy()
    df = df.T
    df['RR'] = df['CONCORDANT'] / (df['CONCORDANT'] + df['DISCORDANT'])
    return df


# In[168]:

def prep_per_sample(df):
    df = df.copy()
    df = df.T
    return df


# In[ ]:

def process_vcf_and_generate_qc_info(fn, out_dir, ipscore_samples, unrelated, pair_samples):
    
    header_end = find_header_end(fn)
    
    count = 0
    progress_level = 0
    progress_increment = 50000
    
    flattened_pairs = list(set([i for sublist in pair_samples for i in sublist]))
    
    num_samples = len(ipscore_samples)
    num_samples_unrel = len(unrelated)
    
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
#     fn_annotated_vcf = os.path.join(out_dir, 'hipstr_output_annot.vcf')
    
#     gt_file = open(os.path.join(out_dir, 'hipstr_pair_gt.tsv'), 'w')
    info_file_gt = open(os.path.join(out_dir, 'hipstr_info_gt.tsv'), 'w')
    info_file_dosage = open(os.path.join(out_dir, 'hipstr_info_dosage.tsv'), 'w')
    
#     info_file_gt_unrel = open(os.path.join(out_dir, 'hipstr_info_gt_unrel.tsv'), 'w')
#     info_file_dosage_unrel = open(os.path.join(out_dir, 'hipstr_info_dosage_unrel.tsv'), 'w')
    
    
    replication_gt_file = open(os.path.join(out_dir, 'hipstr_replication_gt_info.tsv'), 'w')
    replication_length_file = open(os.path.join(out_dir, 'hipstr_replication_length_info.tsv'), 'w')
    
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
    
    per_pair_data_cumulative_gt = {"_".join(tp):{'CONCORDANT':0, 'DISCORDANT':0, 'REF':0, 'MISSING':0} for tp in pair_samples}
    
    per_pair_data_cumulative_length = {"_".join(tp):{'CONCORDANT':0, 'DISCORDANT':0, 'REF':0, 'MISSING':0} for tp in pair_samples}
    

  
    count_missing = 0
    num_variants = 0
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Starting Variant Processing: {}".format(ts)
    
    
    count_num_length = 0
    for line in F:
        
        if progress_level == progress_increment:
            d = datetime.datetime.now()
            ts = d.strftime('%D- %H:%M')
            print "processed {} variants {}".format(count, ts)
            progress_level = 0
            
        progress_level +=1

        line = line.rstrip()
        lin_spl = line.split('\t')

        if count < header_end-1:
            count +=1 
            continue

        if count == header_end-1:



            header = copy.deepcopy(lin_spl)
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
#             print header_dict

            samples = header[9:]
            info_cols = header[:9]
            per_sample_data_cumulative_noxy = {s:{} for s in samples}
            per_sample_data_cumulative = {s:{} for s in samples}
    
    



            # may not use these
            cols_af =  ['NREF', 'REF', 'NMissing', 'NNREF_UUIDs','NNREF_AF', 'ALLELES_DIST', 'FRAC_MISSING']
            cols_af_i2QTL_unrel = ["{}_i2QTL_unrel".format(i) for i in cols_af]
            
            
            
            cols_per_samp_max_mins = ['MIN_Q', 'MAX_FLANKINDEL', 'MAX_DSTUTTER', 'MIN_AB', 'MIN_FS',
                                     'MAX_DOSAGE', 'MIN_DOSAGE']
            cols_info_extracted = ['DP','DSTUTTER', 'DFLANKINDEL', 'PERIOD', 'NFILT', 'INFRAME_UP',
                                   'INFRAME_DOWN']
            
            cols_info_generated = ['FRAC_DP', 'FRAC_DSTUTTER', 'FRAC_DFLANKINDEL']

            cols_replication = ['NP_VAR', 'NP_CONC', 'NP_DISC', 'NP_MISSING', 'RR']


            genotypes_header = ['ID'] + flattened_pairs
            dosages_header = ['ID'] + flattened_pairs
            info_file_header = (['ID'] + cols_af +  cols_af_i2QTL_unrel + cols_per_samp_max_mins + 
                                cols_info_extracted + cols_info_generated)
            
            
            
            
            replication_file_header = (['ID'] + cols_replication + cols_per_samp_max_mins + 
                                       cols_info_extracted + cols_info_generated + ['PERC_MISSING'])
            


            #====== Write Headers For Each File ========
            info_file_gt.write("\t".join(info_file_header) + '\n')
            info_file_dosage.write("\t".join(info_file_header) + '\n')
            
            replication_gt_file.write("\t".join(replication_file_header) + '\n')
            replication_length_file.write("\t".join(replication_file_header) + '\n')


        elif count > header_end:

            if count == header_end+1:

                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}
#                 print format_dict

            num_variants +=1

            info_col = header_dict['INFO']

            POS = lin_spl[header_dict['POS']]
            ID = lin_spl[header_dict['ID']]
            info = lin_spl[header_dict['INFO']]
            chrom = str(lin_spl[header_dict['#CHROM']])

            # for all samples

            fields_desired = ['GT', 'GB', 'Q', 'DP', 'DSTUTTER', 'DFLANKINDEL', 'AB', 'FS']
            geno_fields_dict = {fields_desired[i]:i for i in range(0, len(fields_desired))}

            geno_extracted, missing_samps = get_geno_fields(fields_desired, samples, lin_spl, 
                                                            header_dict, format_dict)
            gb_dict = geno_extracted[geno_fields_dict['GB']]
            gt_dict = geno_extracted[geno_fields_dict['GT']]
            
            
         
            

            # compute summed dosages from haplotypes
            try:
                
                dosage_dict = {(s):(sum(map(int, x.split("|"))) if x != './.' else "./.") for s,x 
                            in gb_dict.iteritems()}
            except:
                print gb_dict
                return
        
            subset_missing = [i for i in missing_samps if i in ipscore_samples]
            num_missing = len(subset_missing)

#             num_missing = len(missing_samps)

            num_passing = num_samples - num_missing
            
          
            
            percent_missing = round((num_missing/num_samples), 3)
            # lets omit sites with all missings 
            if percent_missing == 1:
                count_missing +=1
            if percent_missing < 0.2:
                

                #===== Allele Freq Data Extraction =============
                (data_dict, out_names, 
                 OUT_AF) = calculate_alt_allele_freq(gt_dict, ipscore_samples, chrom, num_samples)
                
                (data_dict_unrel, out_names_unrel, 
                 OUT_AF_unrel) = calculate_alt_allele_freq(gt_dict, unrelated, chrom, num_samples_unrel)
                
                (z, n, 
                 OUT_AF_DOSAGE) = calculate_alt_allele_freq(dosage_dict, ipscore_samples, chrom, num_samples)
                
                
                (z_unrel, n_unrel, 
                 OUT_AF_DOSAGE_unrel) = calculate_alt_allele_freq(dosage_dict, unrelated, chrom, num_samples_unrel)


                # per sample stats
                try:
                    max_Q, min_Q = get_max_min_dict(geno_extracted[geno_fields_dict['Q']], ipscore_samples)
                except:
                    print num_missing, num_samples, num_passing, percent_missing, geno_extracted[geno_fields_dict['Q']]
                    
                max_DS, min_DS = get_max_min_dict(geno_extracted[geno_fields_dict['DSTUTTER']], ipscore_samples)
                max_FI, min_FI = get_max_min_dict(geno_extracted[geno_fields_dict['DFLANKINDEL']], ipscore_samples)
                max_AB, min_AB = get_max_min_dict(geno_extracted[geno_fields_dict['AB']], ipscore_samples)
                max_FS, min_FS = get_max_min_dict(geno_extracted[geno_fields_dict['FS']], ipscore_samples)
                max_dosage, min_dosage = get_max_min_dict(dosage_dict, ipscore_samples)
                # output line
                OUT_PS = [min_Q, max_FI, max_DS, min_AB, min_FS, max_dosage, min_dosage]
                   
                # ======= Add info to cumulative dict ==========
                add_counts_to_per_sample(dosage_dict, samples, per_sample_data_cumulative)

                #====== INFO Col extractions ========================
            
                
                cols_dict = dict(zip(cols_info_extracted, range(0, len(cols_info_extracted))))
                info_out = []
                info_cols_dict = {}
                for l in cols_info_extracted:
                    c = parse_info_col(info, l)
                    info_out.append(c)
                    info_cols_dict[l] = c
                #====================================================
                # aggregate locus stats
                DP = int(info_cols_dict['DP'])
                DSTUTTER = int(info_cols_dict['DSTUTTER'])
                DFLANKINDEL = int(info_cols_dict['DFLANKINDEL'])
                PERIOD = info_cols_dict['PERIOD']
                NFILT = info_cols_dict['NFILT']
                INFRAME_UP = info_cols_dict['INFRAME_UP']
                INFRAME_DOWN  = info_cols_dict['INFRAME_DOWN']
                
                perc_dp = round((DP/num_passing), 4)
                perc_dstutter = round((DSTUTTER/DP),4)
                perc_flankindel = round((DFLANKINDEL/DP), 4)

                OUT_AGG = [DP, DSTUTTER, DFLANKINDEL, PERIOD, NFILT, INFRAME_UP, INFRAME_DOWN, 
                           perc_dp, perc_dstutter, perc_flankindel]


                #====== Write to info file =========

                OUT_LINE_INFO_GT = [ID] + OUT_AF + OUT_AF_unrel + OUT_PS + OUT_AGG 
#                 OUT_LINE_INFO_GT_UNREL = [ID] + OUT_AF_unrel
                
                
                OUT_LINE_INFO_DOSAGE = [ID] + OUT_AF_DOSAGE +  OUT_AF_DOSAGE_unrel + OUT_PS + OUT_AGG 
#                 OUT_LINE_INFO_DOSAGE_UNREL = [ID] + OUT_AF_DOSAGE_unrel
                
                info_file_gt.write("\t".join(map(str, OUT_LINE_INFO_GT)) + '\n') 
                info_file_dosage.write("\t".join(map(str, OUT_LINE_INFO_DOSAGE)) + '\n') 
                
                
                info_file_gt.write("\t".join(map(str, OUT_LINE_INFO_GT)) + '\n') 
                info_file_dosage.write("\t".join(map(str, OUT_LINE_INFO_DOSAGE)) + '\n') 
                

                #================TWIN/PAIR Concordance Calculations ==========
                
                # check for non_missing_gt
                if chrom not in ['X', 'Y']:
                    add_counts_to_per_sample(dosage_dict, samples, per_sample_data_cumulative_noxy)

                    t = [gt_dict[s] for s in flattened_pairs if gt_dict[s] != './.']
                    if len(t) > 0:
                        pair_data_gt = calculate_pair_concordance(pair_samples, gt_dict)
                        pair_data_length = calculate_pair_concordance(pair_samples, dosage_dict)
                        


                        if pair_data_gt != False:
                            
                            (per_pair_data_gt, number_with_var, 
                             number_concordant_with_var, number_discordant_with_var,
                             number_missing, replication_rate) = pair_data_gt
                            
                            OUT_PAIR_DATA_GT = copy.deepcopy([number_with_var, number_concordant_with_var,
                                             number_discordant_with_var, number_missing, replication_rate])
                            
                            
                            for tp, d_gt in zip(pair_samples, per_pair_data_gt):
                                p = "_".join(tp)
                                per_pair_data_cumulative_gt[p][d_gt] = (per_pair_data_cumulative_gt[p].get(d_gt, 0) + 1)
                                
                            OUT_LINE_REPLICATION_GT = ([ID] + OUT_PAIR_DATA_GT +  OUT_PS + 
                                                       OUT_AGG + [percent_missing])
                            
                            replication_gt_file.write("\t".join(map(str,OUT_LINE_REPLICATION_GT)) + '\n')
                            
                            

                        if pair_data_length != False:
                            
                            
                            (per_pair_data_length, number_with_var, 
                             number_concordant_with_var, number_discordant_with_var,
                             number_missing, replication_rate) = pair_data_length
                         
                            
                            OUT_PAIR_DATA_LENGTH = copy.deepcopy([number_with_var, number_concordant_with_var,
                                             number_discordant_with_var, number_missing, replication_rate])
                           
                             
                            for tp, d_len in zip(pair_samples, per_pair_data_length):
                                p = "_".join(tp)
                                per_pair_data_cumulative_length[p][d_len]=(per_pair_data_cumulative_length[p].get(d_len, 0) + 1)
                            
                            OUT_LINE_REPLICATION_LENGTH = ([ID] + OUT_PAIR_DATA_LENGTH +  OUT_PS + 
                                                       OUT_AGG + [percent_missing])
                            
                            replication_length_file.write("\t".join(map(str,OUT_LINE_REPLICATION_LENGTH)) + '\n')
                                
                            
                                                        
            
            else:
                count_missing +=1
                      


                            
                            
#         if count == 10000:
#             replication_gt_file.close()
#             replication_length_file.close()
#             info_file_gt.close()
#             info_file_dosage.close()


#             per_sample_gt = pd.DataFrame(per_pair_data_cumulative_gt).pipe(prep_per_sample)
#             per_sample_length = pd.DataFrame(per_pair_data_cumulative_length).pipe(prep_per_sample)

#             fn_per_sample_gt = os.path.join(out_dir, 'hipstr_per_sample_gt.tsv')
#             fn_per_sample_length = os.path.join(out_dir, 'hipstr_per_sample_length.tsv')

#             per_sample_gt.to_csv(fn_per_sample_gt, sep = '\t')
#             per_sample_length.to_csv(fn_per_sample_length, sep = '\t')

#             print "number of variants: {}".format(num_variants), "number_all_missing: {}".format(count_missing)

#             d = datetime.datetime.now()
#             ts = d.strftime('%D- %H:%M')
#             print "Completed Variant Processing: {}".format(ts)
#             print 'TEST COMPLETE'
#             return per_sample_gt, per_sample_length
            
        
            
        count +=1
    
    replication_gt_file.close()
    replication_length_file.close()
    info_file_gt.close()
    info_file_dosage.close()


    per_pair_gt = pd.DataFrame(per_pair_data_cumulative_gt).pipe(prep_per_pair)
    per_pair_length = pd.DataFrame(per_pair_data_cumulative_length).pipe(prep_per_pair)
    
    per_sample_length  = pd.DataFrame(per_sample_data_cumulative).pipe(prep_per_sample)
    per_sample_noxy_length  = pd.DataFrame(per_sample_data_cumulative_noxy).pipe(prep_per_sample)
    
    fn_per_pair_gt = os.path.join(out_dir, 'hipstr_per_pair_gt.tsv')
    fn_per_pair_length = os.path.join(out_dir, 'hipstr_per_pair_length.tsv')
    fn_per_sample_length = os.path.join(out_dir, 'hipstr_per_sample_length.tsv')
    fn_per_sample_noxy_length = os.path.join(out_dir, 'hipstr_per_sample_noxy_length.tsv')

    
    per_pair_gt.to_csv(fn_per_pair_gt, sep = '\t')
    per_pair_length.to_csv(fn_per_pair_length, sep = '\t')
    per_sample_length.to_csv(fn_per_sample_length, sep = '\t')
    per_sample_noxy_length.to_csv(fn_per_sample_noxy_length, sep = '\t')
    
    
    
    print "number of variants: {}".format(num_variants), "number_all_missing: {}".format(count_missing)

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
    
    
#     parser.add_argument("-gender", "--gender_map", dest="gender_map", metavar='<gender.txt>', help="gender_map_file", required=True)
        
    
    parser.add_argument("-unr", "--unrelated_samples", dest="unr_samples", metavar='<unrelated_samples.txt>', help="list of sample names that indicate unrelated samples", required=True)
        
       
#     parser.add_argument("-non_ipsc", "--non_ipsc_samples", dest="non_ipsc_samples", metavar='<non_ipsc_samples.txt>', help="list of sample names that indicate samples that are not iPSC", required=True)
  
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to extract info and calculate replication rate for HipSTR vcf files')
    add_arguments_to_parser(parser)
    return parser



# In[ ]:

def run_from_args(args):
    vcf_fn = args.vcf_file
    samples_fn = args.samples_fn
    unrelated_samples_fn = args.unr_samples
    
    pair_fn = args.pairs_fn
    out_dir = args.output_dir
    pairs = [line.rstrip().split() for line in open(pair_fn)]
    samples = [line.rstrip() for line in open(samples_fn)]
    unr_samples =  [line.rstrip() for line in open(unrelated_samples_fn)]
    
    process_vcf_and_generate_qc_info(vcf_fn, out_dir, samples, unr_samples, pairs)
    


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

