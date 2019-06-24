
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

import pandas as pd
import csv

import copy 


import itertools
import tempfile
import six
import networkx as nx
import scipy.stats as stats
import argparse

import datetime
import gzip


# In[2]:

def find_header_end(fn):
    """ find header end of vcf or gzipped vcf"""
    
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
                return 'incomplete or missing header, or very long header'
        except:
            break
    
    F.close()

def parse_info_col(t, lab, type_out = str, alternative = 'None'):
    """ parse the info column"""
    
    t = t.split(';')
    
    # remove any of the tags that have no label
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

def bp_precision(x):
    beg, end =[int(i) for i in x.split(',')]
    interval = abs(end - beg)
    return interval


def SVLEN_convert(x):
    try:
        return np.abs(int(x))
    except:
        return 0

def parse_format_col(t, col, lab, type_out = str, alternative = 'None', change_dict = False):
    """ use the format column to pull out info of a specific VCF record
    in a genotyping matrix (unparsed)
    change_dict- dict to return something else for any given value in the format col
    """
    
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


def coord_extract(df,chrom, start, end, contained = True):
    
    if contained:
        return df[(df.Chr==chrom) & (df.POS >= start) & (df.END <= end)]


# In[3]:

def add_annotations_extract_gts(fn, out_dir, unrelated_samples, non_ipsc_samples, sex_dict):
    
    header_end = find_header_end(fn)
    
    count = 0
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
    
    gt_file = open(os.path.join(out_dir, 'melt_gt.tsv'), 'w')
    info_file = open(os.path.join(out_dir, 'melt_info.tsv'), 'w')
    
    fixed_header_vcf = open(os.path.join(out_dir, 'melt_final.vcf'), 'w')
    filtered_vcf = open(os.path.join(out_dir, 'melt_final_filt.vcf'), 'w')
    
    gt_filt_file = open(os.path.join(out_dir, 'melt_gt_filt.tsv'), 'w')
#     ab_filt_file = open(os.path.join(out_dir, 'lumpy_ab_filt.tsv'), 'w')
    
    info_filt_file = open(os.path.join(out_dir, 'melt_info_filt.tsv'), 'w')
    
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
    
    add_filters = False
    add_info = False
    for line in F:
        line = line.rstrip()
                         
        lin_spl = line.split('\t')
                         
        if count < header_end-1:
            if not add_info:
                if line.find('##INFO') == 0:
                    add_info = True
                    l= '##INFO=<ID=1KGP_ID,Number=1,Type=String,Description="ID indicating site found previously in 1KGP (site included in priors)">'
                    fixed_header_vcf.write(l + '\n')
                    filtered_vcf.write(l + '\n')
            
            if not add_filters:
                if line.find('##FILTER') == 0:
                    add_filters = True

                    l2= '##FILTER=<ID=iPSC,Description="Filter indicating site present only in iPSC">'
                    l1 = '##FILTER=<ID=ASSESS,Description="Filter indicating site has ASSESS tranche less than 5">'
                    l3 = '##FILTER=<ID=NMissing,Description="at least 10% of samples have a missing GT (iPSCORE fib/blood + HipSci fib)>'
                    l4 =  '##FILTER=<ID=OTHER_CONTIG,Description="not on autosomes or sex chroms">'
                    
                    l = "\n".join([l1, l2, l3, l4])
                    
                    fixed_header_vcf.write(l + '\n') 
                    filtered_vcf.write(l + '\n')
                    
            fixed_header_vcf.write(line + '\n')
            filtered_vcf.write(line + '\n')
        
                  
        if count == header_end-1:
 
            # fix the .mdup issue in the sample naming from processing
            
            lin_spl = [i.replace('.mdup','') for i in lin_spl]
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            
            samples = lin_spl[9:]
            info_cols = lin_spl[:9]
            
            cols_maf = ['NNREF', 'NNREF_UUIDs','ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'UUIDS_REF',
                        'NMissing']
            cols_extracted = ['SVTYPE', 'INTERNAL', 'RP', 'LP', 'SR', 'RA', 'SVLEN', 'ASSESS', 'TSD','MEINFO'] 
            cols_maf_unr = [i+ '_unr' for i in cols_maf]
            cols_maf_non_ipsc = [i+ '_non_ipsc' for i in cols_maf]
            
            #remove the info column for space savings
        
            
            info_header = ['NAME'] +info_cols[:7] + cols_extracted + ['END','Missing_GTs'] + cols_maf + cols_maf_unr + cols_maf_non_ipsc + ['FILTER_ORIGINAL']
            
            genotypes_header = ['NAME'] + samples
            
            gt_file.write("\t".join(genotypes_header) + '\n')
            gt_filt_file.write("\t".join(genotypes_header) + '\n')
            
            info_file.write("\t".join(info_header) + '\n')
            info_filt_file.write("\t".join(info_header) + '\n')
            
            fixed_header_vcf.write( "\t".join(lin_spl) + '\n')
            filtered_vcf.write( "\t".join(lin_spl) + '\n')
            
            
        elif count > header_end:
            
            if count == header_end+1:
            
                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}
                  
            
            sample_gts = []
            # original first 9 columns
            
            info_col = header_dict['INFO']
            
            
            first_cols = lin_spl[:7]
            
            POS = lin_spl[header_dict['POS']]
            FILTER = lin_spl[header_dict['FILTER']]
            FILTER_ORIGINAL = copy.deepcopy(FILTER)
            
            ID = lin_spl[header_dict['ID']]
            ALT = lin_spl[header_dict['ALT']]
            vtype = ALT.split(':')[-1][:-1]
            info = lin_spl[header_dict['INFO']]
            chrom = lin_spl[header_dict['#CHROM']]
            
            # new ID
#             ID_MOD = "{}_{}_{}".format(vtype, chrom, POS)
            
            
    
            cols_to_parse =  ['SVTYPE', 'INTERNAL', 'RP', 'LP', 'SR', 'RA', 'SVLEN', 'ASSESS', 'TSD', 'MEINFO'] 
            cols_dict = dict(zip(cols_to_parse, range(0, len(cols_to_parse))))
            
            info_out = []
            
            info_cols_dict = {}
            for l in cols_to_parse:
                c = parse_info_col(info, l)
                info_out.append(c)
                info_cols_dict[l] = c

            ASSESS = int(info_cols_dict['ASSESS'])
            
            # give the variants more useful IDs 
            END = int(POS) + 1
            ID_MOD = "{}_{}_{}_{}".format(vtype, chrom, POS, END)
            
            if ID != '.':
                lin_spl[header_dict['INFO']] += ';1KGP_ID={}'.format(ID)
                # swap in my ID because its more useful- but add in the 1KGP_ID to the info column
                lin_spl[header_dict['ID']] = ID_MOD
                
            else:
                lin_spl[header_dict['ID']] = ID_MOD
                
            
            
            # add minor allele frequency:
            
            gts_out = []
            
            # all samples in this case
            
            for s in samples:
                d = lin_spl[header_dict[s]].split(':')
                gt = d[format_dict['GT']]
                gts_out.append(gt)
            
            
            missing_gts = ('./.' in gts_out)
          
            
            # minor allele frequencies
            maf_data_all = calculate_minor_allele_freq(gts_out, samples, chrom, sex_dict)
        
            s = [gts_out[samples.index(i)] for i in unrelated_samples]
            maf_data_unr = calculate_minor_allele_freq(s, unrelated_samples, chrom, sex_dict)
 
            s = [gts_out[samples.index(i)] for i in non_ipsc_samples] 
            maf_data_non_ipsc = calculate_minor_allele_freq(s, non_ipsc_samples, chrom, sex_dict)
            
            # write out data
            num_missing_non_ipsc = int(maf_data_non_ipsc[8])
#             if ID_MOD == 'ALU_1_24328785_24328786':
#                 print num_missing_non_ipsc
            perc_missing = (num_missing_non_ipsc/478)
            chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
        
        
            if str(chrom) not in chroms:
                if FILTER == 'PASS':
                    FILTER_MOD = 'OTHER_CONTIG'
                else:
                    FILTER_MOD = FILTER + ';OTHER_CONTIG'

                lin_spl[header_dict['FILTER']] = FILTER_MOD
                FILTER = FILTER_MOD
                
            
            if (perc_missing > 0.1):
                if FILTER == 'PASS':
                    FILTER_MOD = 'NMissing'
                else:
                    FILTER_MOD = FILTER + ';NMissing'
                
                lin_spl[header_dict['FILTER']] = FILTER_MOD
                FILTER = FILTER_MOD
                
            
            
            if ASSESS < 5:
                if FILTER == 'PASS':
                    FILTER_MOD = 'ASSESS'
                else:
                    FILTER_MOD = FILTER + ';ASSESS'
                
                lin_spl[header_dict['FILTER']] = FILTER_MOD
                FILTER = FILTER_MOD
                
            
            
            if int(maf_data_non_ipsc[0]) == 0:
                if FILTER == 'PASS':
                    FILTER_MOD = 'iPSC'
                else:
                    FILTER_MOD = FILTER + ';iPSC'
            
                lin_spl[header_dict['FILTER']] = FILTER_MOD
                
                FILTER = FILTER_MOD
            
            
            first_cols = lin_spl[:7]
            out_line_info =  [ID_MOD] + first_cols + info_out + [END] + [missing_gts] + maf_data_all + maf_data_unr + maf_data_non_ipsc + [FILTER_ORIGINAL]
            out_line_info = map(str, out_line_info)
            out_line_gt = [ID_MOD] + gts_out
            
            # filtered_data
            if FILTER == 'PASS':
                info_filt_file.write("\t".join(out_line_info) + '\n')
                gt_filt_file.write("\t".join(out_line_gt) + '\n')
                filtered_vcf.write("\t".join(lin_spl) + '\n')
                
            
            info_file.write("\t".join(out_line_info) + '\n')
            
            gt_file.write("\t".join(out_line_gt) + '\n')
            
            fixed_header_vcf.write("\t".join(lin_spl) + '\n')
           
        
        
        count +=1
    
    # close out files
    F.close()
    gt_file.close() 
    fixed_header_vcf.close()
    info_file.close()
    gt_filt_file.close()
    info_filt_file.close()
    filtered_vcf.close()


# In[1]:

def safe_div(x, y, alt=0):
    try:
        return x/y
    except:
        return alt


# In[4]:

def calculate_minor_allele_freq(gts, uuids, chrom, sex_dict):
    """ calculate minor allele frequency, assume missing genotypes are reference"""
    
    


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
    gts_dict = {'0/0':0, '0/1':1, '1/1':2, './.':0}
    gts_dict_sex = {'0/0':0, '0/1':1, '1/1':1, './.':0}
    missing = 0
    
    count = 0
    for c in gts:
        uuid = uuids[count]
        sex = sex_dict[uuid]
        
        if c == './.':
            missing +=1
            # skip onward- we are ignoring missing
            continue
        

        if chrom=='Y':
            if sex == 'M':
                aa = gts_dict_sex[c]
                ra = 1 - aa
                
            else:
                count +=1
                continue



        elif chrom=='X':
            if sex == 'M':
                aa = gts_dict_sex[c]
                ra = 1 - aa
            else:
                aa = gts_dict[c]
                ra = 2 - aa


        else:
            aa = gts_dict[c]
            ra = 2 - aa
#                 ref_allele += ra
#                 alt_allele += aa

        ref_allele += ra
        alt_allele += aa

        if c not in ['0/0','./.']:
            non_ref +=1
            non_ref_uuids.append(uuids[count])

        if c not in ['1/1','0/1']:
            ref +=1
            ref_uuids.append(uuids[count])




        count +=1



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
            
    
#     ref_af = ref_allele/tot_alleles
#     alt_af = alt_allele/tot_alleles

#     afs = [ref_af, alt_af]
#     maf = min(afs)
    

#     if afs.index(maf) == 0:
#         minor_allele = 'REF'
#     else:
#         minor_allele = 'ALT'
    

    out_names = ['NNREF', 'NNREF_UUIDs','ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'UUIDS_REF', 'Missing']
    out_data = [non_ref, ",".join(non_ref_uuids), alt_af, ref_af, maf, minor_allele, ref, ",".join(ref_uuids),
               missing]
    
    data_dict = dict(zip(out_names, out_data))
    
    
    return out_data


# In[217]:

def add_arguments_to_parser(parser):
    

    parser.add_argument("-vcf", "--vcf", dest="vcf_file", metavar='<vcf_file>', help="vcf file from lumpy/speedseq pipeline, may be gzipped or not", required=True)

    
    parser.add_argument("-gender", "--gender_map", dest="gender_map", metavar='<gender.txt>', help="gender_map_file", required=True)
        
    
    parser.add_argument("-unr", "--unrelated_samples", dest="unr_samples", metavar='<unrelated_samples.txt>', help="list of sample names that indicate unrelated samples", required=True)
    
      
    parser.add_argument("-non_ipsc", "--non_ipsc_samples", dest="non_ipsc_samples", metavar='<non_ipsc_samples.txt>', help="list of sample names that indicate samples that are not iPSC", required=True)
  
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    
    
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to annotate extracted info and genotypes from lumpy/melt with various things such as gene intersections, minor allele frequency etc.')
    add_arguments_to_parser(parser)
    return parser


# In[ ]:

def run_from_args(args):

    unrelated_fn = args.unr_samples
    non_ipsc_fn = args.non_ipsc_samples

    unr_samples = [line.rstrip() for line in open(unrelated_fn)]
    non_ipsc_samples = [line.rstrip() for line in open(non_ipsc_fn)]
    
    
    gender_fn = args.gender_map
    sex_dict = {line.rstrip().split()[0]:line.rstrip().split()[1] for line in open(gender_fn)}
    
    t = datetime.datetime.now().strftime('%c')
    print "extracting variants, annotating filters and MAFs: ", t
    
    add_annotations_extract_gts(args.vcf_file, args.output_dir, unr_samples, non_ipsc_samples, sex_dict)
    
    t = datetime.datetime.now().strftime('%c')
    print 'complete:', t
    


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))

