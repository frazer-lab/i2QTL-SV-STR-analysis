
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

def add_dosage_to_per_sample(d, samples, cumulative_dict):
    
    for s in samples:
        c = str(d[s])
        if c in ['./.', '.']:
            continue
        else:
            c = int(c)
            
            if c > 0:
                cumulative_dict[s]['INS'] = cumulative_dict[s]['INS'] + c
                
            elif c < 0:
                cumulative_dict[s]['DEL'] = cumulative_dict[s]['DEL'] + abs(c)
                
            else:
                pass



def process_vcf_and_generate_qc_info(fn, out_dir):
    
    header_end = find_header_end(fn)
    
    count = 0
    progress_level = 0
    progress_increment = 50000
    

    
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    

    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
  
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

            samples = header[9:]
            num_samples = len(samples)
            
            info_cols = header[:9]
            per_sample_data_cumulative_noxy = {s:{'INS':0, 'DEL':0} for s in samples}
            per_sample_data_cumulative = {s:{'INS':0, 'DEL':0} for s in samples}
        


        elif count > header_end:

            if count == header_end+1:

                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}


            num_variants +=1
            info_col = header_dict['INFO']

            POS = lin_spl[header_dict['POS']]
            ID = lin_spl[header_dict['ID']]
            info = lin_spl[header_dict['INFO']]
            chrom = str(lin_spl[header_dict['#CHROM']])

            # for all samples

            fields_desired = ['GB']
            geno_fields_dict = {fields_desired[i]:i for i in range(0, len(fields_desired))}

            geno_extracted, missing_samps = get_geno_fields(fields_desired, samples, lin_spl, 
                                                            header_dict, format_dict)
            gb_dict = geno_extracted[geno_fields_dict['GB']]
            
            try:
                dosage_dict = {(s):(sum(map(int, x.split("|"))) if x != './.' else "./.") for s,x 
                            in gb_dict.iteritems()}
            except:
                print gb_dict
            
            
            if len(missing_samps) < num_samples:
                add_dosage_to_per_sample(dosage_dict, samples, per_sample_data_cumulative)
                
                # check for non_missing_gt
                if chrom not in ['X', 'Y']:
                    add_dosage_to_per_sample(dosage_dict, samples, per_sample_data_cumulative_noxy)

            
            else:
                count_missing +=1
                      

            
        count +=1
    

    per_sample_length  = pd.DataFrame(per_sample_data_cumulative).pipe(prep_per_sample)
    per_sample_noxy_length  = pd.DataFrame(per_sample_data_cumulative_noxy).pipe(prep_per_sample)
    
    fn_per_sample_length = os.path.join(out_dir, 'hipstr_per_sample_dosage_total.tsv')
    fn_per_sample_noxy_length = os.path.join(out_dir, 'hipstr_per_sample_noxy_dosage_total.tsv')
    
    per_sample_length.to_csv(fn_per_sample_length, sep = '\t')
    per_sample_noxy_length.to_csv(fn_per_sample_noxy_length, sep = '\t')
    
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Completed Variant Processing: {}".format(ts)
    print 'COMPLETE'
    return
                                

def add_arguments_to_parser(parser):
    

    parser.add_argument("-vcf", "--vcf", dest="vcf_file", metavar='<vcf_file>', help="vcf file from hipstr, may be gzipped or not", required=True)
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    parser.set_defaults(entry_point=run_from_args)

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to extract info and calculate replication rate for HipSTR vcf files')
    add_arguments_to_parser(parser)
    return parser

def run_from_args(args):
    vcf_fn = args.vcf_file
    out_dir = args.output_dir
    
    process_vcf_and_generate_qc_info(vcf_fn, out_dir)

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))

