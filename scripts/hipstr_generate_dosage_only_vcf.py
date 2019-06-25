
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
from collections import Counter


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


# In[1]:

def MinMaxScaler(features, minimum = 0, maximum = 2):
    min_features = min(features)
    max_features = max(features)
    
    if max_features - min_features == 0:
        return
    
    def transform(x, min_features, max_features, min_target, max_target):
                
        xstd = (x - min_features) / (max_features - min_features)
        xscaled = xstd * (max_target - min_target) + min_target
        
        return xscaled
        
    features_out = [transform(i, min_features, max_features, minimum, maximum) for i in features]
    return features_out

def get_rank_dict(ds_dict, sample_subset, type_out = int):
    """adjusted dosages - rank the dosages and put them on a 0-2 scale"""
    
    values = []
    for k,v in ds_dict.iteritems():
        
        if k in sample_subset:
            if v in ['./.', '.']:
                pass
            else:
                values.append(type_out(v))

    rank_vals = stats.rankdata(values)
    
    scale_ranks = MinMaxScaler(rank_vals)
    scale_ranks = [round(i, 4) for i in scale_ranks]
    rank_dict = dict(zip(map(str, values), scale_ranks))
    return rank_dict

def get_adjusted_ds_dict(ds_dict, rank_dict, type_in = int):
    
    adj_ds_dict = {}
    
    for k,v in ds_dict.iteritems():
        try:
            key_out = str(type_in(v))
        except:
            key_out = v
            
        rank_convert = rank_dict.get(key_out, -1)
        adj_ds_dict[k] = rank_convert
        
    return adj_ds_dict
  


# In[36]:

def compute_maf_dosage(ds_dict, samples):
    s = set(samples)
    ds_dict_sub = {k:v for k,v in ds_dict.iteritems() if k in s}
    vals = ds_dict_sub.values()
    ignore = {'./.', '.'}
    vals_corrected = [i for i in vals if i not in ignore]
    C = Counter(vals_corrected)
    mc = C.most_common()
    mc_allele, count_mc_allele = mc[0]
    
    num_passing = len(vals_corrected)
    num_alt = count_mc_allele -  num_passing
    if num_alt > 0:
        maf = num_alt/num_passing
    
    else:
        maf = 0
    
    return maf, count_mc_allele, mc_allele, num_passing, num_alt


# In[3]:

def process_vcf_and_annotate(fn, out_dir, sample_subset, fn_out):
    
    header_end = find_header_end(fn)
    
    count = 0
    progress_level = 0
    progress_increment = 50000
    
    
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
    fn_annotated_vcf = os.path.join(out_dir, fn_out)
#     fn_maf = os.path.join(out_dir, 'hipstr_maf_rna.tsv')
    
    annotated_vcf = open(fn_annotated_vcf, 'w')
#     maf_file = open(fn_maf, 'w')
    
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
    

#     IDs_set = set(IDs)
    
    count_missing = 0
    num_variants = 0
    num_variants_out = 0
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Starting Variant Processing: {}".format(ts)
    
#     info_fields = 'INFRAME_PGEOM INFRAME_UP INFRAME_DOWN OUTFRAME_PGEOM OUTFRAME_UP OUTFRAME_DOWN BPDIFFS START END AN REFAC AC NSKIP NFILT DP DSNP DSTUTTER DFLANKINDEL'.split()
            
    add_id = False
    add_info = False
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
            if not add_id:
                if line.find('##FORMAT=<ID=DS') == 0:
                    add_id= True
                    
            else:
                # remove the format column
                if line.find('##FORMAT') == 0:
                    count += 1
                    continue
            
            if not add_info:
                if line.find('##INFO') == 0:
                    add_info= True
         
                    
            annotated_vcf.write(line + '\n')
            count +=1 
            continue

        if count == header_end-1:

            
            header = copy.deepcopy(lin_spl)
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            info_head = header[:9]
            samples = header[9:]
            info_cols = header[:9]
            sample_order_dict = {s:header.index(s) for s in samples}
            
            try:
                reodered_samples = [header[sample_order_dict[i]] for i in sample_subset]
            except:
                print sample_order_dict
                print desired_samples
                print samples
                break
            
            header_mod = info_cols + reodered_samples
            annotated_vcf.write("\t".join(header_mod) + '\n')

            

        elif count > header_end:

            if count == header_end+1:

                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}
#                 print format_dict


            POS = lin_spl[header_dict['POS']]
            ID = lin_spl[header_dict['ID']]
            INFO = lin_spl[header_dict['INFO']]
            chrom = str(lin_spl[header_dict['#CHROM']])
            
            cols_info_extracted = ['MAF_RNA']
            cols_dict = dict(zip(cols_info_extracted, range(0, len(cols_info_extracted))))
            info_out = []
            info_cols_dict = {}
            for l in cols_info_extracted:
                c = parse_info_col(INFO, l)
                info_out.append(c)
                info_cols_dict[l] = c
            
            maf = float(info_cols_dict['MAF_RNA'])
            
            if maf >= 0.01:
                
            # for all samples

                fields_desired = ['DS']
                geno_fields_dict = {fields_desired[i]:i for i in range(0, len(fields_desired))}

                geno_extracted, missing_samps = get_geno_fields(fields_desired, samples, lin_spl, 
                                                                header_dict, format_dict)

                ds_dict = geno_extracted[geno_fields_dict['DS']]
                # set the format col to just DS
                lin_spl[header_dict['FORMAT']] = 'DS'
                info_cols = lin_spl[:9]
                reodered_samples = [ds_dict[s] for s in sample_subset]
                out_line_mod = info_cols + reodered_samples
                annotated_vcf.write("\t".join(map(str, out_line_mod)) + '\n')
            
                num_variants_out +=1
        
            num_variants +=1
                      
        count +=1
    annotated_vcf.close()
    print "number of variants processed: {}".format(num_variants)
    print "number of variants outputted: {}".format(num_variants_out)

    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Completed Variant Processing - subsetted to dosage only: {}".format(ts)
    print 'COMPLETE'
    print 'files generated: {}'.format(fn_annotated_vcf)
    return
                                


# In[38]:

def add_arguments_to_parser(parser):
    

    parser.add_argument("-vcf", "--vcf", dest="vcf_file", metavar='<vcf_file>', help="vcf file from lumpy/speedseq pipeline, may be gzipped or not", required=True)
    
    parser.add_argument("-s", "--samples", dest="samples_fn", metavar='<samples_fn>', help="file with samples to use to compute the maf and converted dosages", required=False)
    
    parser.add_argument("-id_fn", "--id_fn", dest="id_fn", metavar='<IDs to keep in the final output>', help="file of pairs of IDs one per line to keep in the output", required=False)
    
    parser.add_argument("-fn_out", "--fn_out", dest="fn_out", metavar='<file name of vcf output (NOT dir)>', help="file name of the vcf file- excluding the directory, just file name", default = False, required=False)
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    parser.set_defaults(entry_point=run_from_args)


# In[48]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to extract info and calculate replication rate for HipSTR vcf files')
    add_arguments_to_parser(parser)
    return parser



# In[47]:

def run_from_args(args):
    vcf_fn = args.vcf_file
    id_fn = args.id_fn
    samples_fn = args.samples_fn
    fn_out = args.fn_out
    out_dir = args.output_dir
    samples = [line.rstrip() for line in open(samples_fn)]
    process_vcf_and_annotate(vcf_fn, out_dir, samples, fn_out)    


# In[87]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))

