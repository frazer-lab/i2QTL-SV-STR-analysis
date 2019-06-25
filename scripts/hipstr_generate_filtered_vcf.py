
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
                number_concordant_with_var +=1
                per_pair_data.append('CONCORDANT')
            else:
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


# In[4]:

def process_vcf_and_annotate(fn, out_dir, IDs, suffix = False, tag_hipsci = False, tag_ipscore = False):
    
    header_end = find_header_end(fn)
    
    count = 0
    progress_level = 0
    progress_increment = 50000
    
    
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
    fn_annotated_vcf = os.path.join(out_dir, 'hipstr_output_annot.vcf')
    annotated_vcf = open(fn_annotated_vcf, 'w')   
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']

    
    if tag_hipsci:
        tag = 'HIPSCI_REF'
    if tag_ipscore:
        tag = 'iPSCORE_REF'
        
    add_tag = False
    for i in [tag_hipsci, tag_ipscore]:
        if i == True:
            add_tag = True
    
    
    IDs_set = set(IDs)
    
    count_missing = 0
    num_variants = 0
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Starting Variant Processing: {}".format(ts)
    
    info_fields = 'INFRAME_PGEOM INFRAME_UP INFRAME_DOWN OUTFRAME_PGEOM OUTFRAME_UP OUTFRAME_DOWN BPDIFFS START END AN REFAC AC NSKIP NFILT DP DSNP DSTUTTER DFLANKINDEL'.split()
            
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
            if suffix:
                if line.find('##INFO') == 0:
                    for p in info_fields:
                        line = line.replace("ID={}".format(p), 'ID={}_{}'.format(p, suffix))
        
            
            if not add_id:
                if line.find('##FORMAT') == 0:
                    add_id= True
                    l = '##FORMAT=<ID=BPS,Number=1,Type=String,Description="Change in base pairs from reference (sum change of both alleles) sum of BPs">'
                    annotated_vcf.write(l + '\n')
                    
                    
                    
            if not add_info:
                if line.find('##INFO') == 0:
                    add_info= True
                    if tag_hipsci == True:

                        l= '##INFO=<ID=HIPSCI_REF,Number=1,Type=String,Description="Tag indicating site found previously in HipSci genomes, genotyped in iPSCORE">'
                    
                    if tag_ipscore == True:  
                        l= '##INFO=<ID=iPSCORE_REF,Number=1,Type=String,Description="Tag indicating site found in iPSCORE genomes, is not genotyped in hipsci">'
                    annotated_vcf.write(l + '\n')
                     
                        

            
            
            
            annotated_vcf.write(line + '\n')
            count +=1 
            continue

        if count == header_end-1:
            annotated_vcf.write(line + '\n')



            header = copy.deepcopy(lin_spl)
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            samples = header[9:]

            

        elif count > header_end:

            if count == header_end+1:

                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}
#                 print format_dict

           
    


            info_col = header_dict['INFO']

            POS = lin_spl[header_dict['POS']]
            ID = lin_spl[header_dict['ID']]
            INFO = lin_spl[header_dict['INFO']]
            chrom = str(lin_spl[header_dict['#CHROM']])
            
            
            
            if ID in IDs_set:
                

                # for all samples

                fields_desired = ['GT', 'GB']
                geno_fields_dict = {fields_desired[i]:i for i in range(0, len(fields_desired))}

                geno_extracted, missing_samps = get_geno_fields(fields_desired, samples, lin_spl, 
                                                                header_dict, format_dict)
                gb_dict = geno_extracted[geno_fields_dict['GB']]
                gt_dict = geno_extracted[geno_fields_dict['GT']]


                # compute summed dosages from haplotypes
                try:

                    dosage_dict = {(s):(sum(map(int, x.split("|"))) if x != './.' else ".") for s,x 
                                in gb_dict.iteritems()}
                except:
                    print gb_dict
                    return


                for s in samples:
                    lin_spl[header_dict[s]] = lin_spl[header_dict[s]] + ':{}'.format(dosage_dict[s])

                FORMAT = lin_spl[header_dict['FORMAT']]
                lin_spl[header_dict['FORMAT']] =  FORMAT + ':BPS'

                if add_tag:
                    lin_spl[header_dict['INFO']] = INFO + ';{}'.format(tag)
                    
                
                # write to file
                
                INFO = lin_spl[header_dict['INFO']]
            
                if suffix:
                    for p in info_fields:
                        INFO = INFO.replace("{}=".format(p), '{}_{}='.format(p, suffix))
                        
                
                    lin_spl[header_dict['INFO']] = INFO
                
                
                line = "\t".join(lin_spl) 
                annotated_vcf.write(line + '\n')
                num_variants +=1
            
            else:
                count_missing +=1
                
            
            
        
                      
        count +=1
    annotated_vcf.close()
    print "number of variants: {}".format(num_variants), "number_all_missing: {}".format(count_missing)

    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    print "Completed Variant Processing: {}".format(ts)
    print 'COMPLETE'
    return
                                


# In[6]:

def add_arguments_to_parser(parser):
    

    parser.add_argument("-vcf", "--vcf", dest="vcf_file", metavar='<vcf_file>', help="vcf file from lumpy/speedseq pipeline, may be gzipped or not", required=True)
    
    parser.add_argument("-id_fn", "--id_fn", dest="id_fn", metavar='<IDs to keep in the final output>', help="file of pairs of IDs one per line to keep in the output", required=True)
    
    parser.add_argument("-suff", "--suffix", dest="suffix", metavar='<Suffix for INFO COL TAGS>', help="file of pairs of IDs one per line to keep in the output", default = False, required=False)
    
    parser.add_argument("-tag_ipscore", dest="tag_ipscore", action =  'store_true',
                        help="tag ipscore variants", default = False, required = False)
    parser.add_argument("-tag_hipsci", dest="tag_hipsci", action =  'store_true',
                        help="tag ipscore variants", default = False, required = False)
    
    
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to extract info and calculate replication rate for HipSTR vcf files')
    add_arguments_to_parser(parser)
    return parser



# In[7]:

def run_from_args(args):
    vcf_fn = args.vcf_file
    id_fn = args.id_fn
    t_i = args.tag_ipscore
    t_h = args.tag_hipsci
    suff = args.suffix
    
    out_dir = args.output_dir
    IDs = [line.rstrip() for line in open(id_fn)]
    
    process_vcf_and_annotate(vcf_fn, out_dir, IDs, suffix= suff, tag_hipsci= t_h, tag_ipscore = t_i)
    


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

