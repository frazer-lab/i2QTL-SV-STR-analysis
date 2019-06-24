
# coding: utf-8

# In[4]:

from __future__ import division
import numpy as np
import os
import sys

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


# In[5]:

import time


# In[1]:


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

def filter_inv_BND(x):
    """ apply filtering similar to GTeX on these 
    variant types (making sure the required amount of evidence is present)"""
    SU = int(x.SU)
    PE = int(x.PE)
    SR = int(x.SR)
    Qual = float(x.MSQ)
    svtype = x.SVTYPE
    

    
    if svtype in ['BND', 'INV']:
        if Qual < 100:
            return True    

        if svtype == 'INV':
            tot = PE + SR
            for s in [PE,SR]:
                frac = s/tot
                if frac < 0.1:
                    return True
                
        
        if svtype == 'BND':
            tot = PE + SR
            for s in [PE,SR]:
                frac = s/tot
                if frac < 0.25:
                    return True
        return False
        
        
    else:
        return False

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


# In[13]:

def add_annotations_extract_gts(fn, out_dir, unrelated_samples, non_ipsc_samples, sex_dict):
    
    header_end = find_header_end(fn)
    
    count = 0
    progress_level = 0
    progress_increment = 50000

    
    
    if fn.split('.').pop() == 'gz':
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
    fn_annotated_vcf = os.path.join(out_dir, 'speedseq_output_annot.vcf')
    
    gt_file = open(os.path.join(out_dir, 'lumpy_gt.tsv'), 'w')
    ab_file = open(os.path.join(out_dir, 'lumpy_ab.tsv'), 'w')
    info_file = open(os.path.join(out_dir, 'lumpy_info.tsv'), 'w')
    annotated_vcf = open(os.path.join(out_dir, 'speedseq_output_annot.vcf'), 'w')
    
    gt_filt_file = open(os.path.join(out_dir, 'lumpy_gt_filt.tsv'), 'w')
    ab_filt_file = open(os.path.join(out_dir, 'lumpy_ab_filt.tsv'), 'w')
    
    info_filt_file = open(os.path.join(out_dir, 'lumpy_info_filt.tsv'), 'w')
    
    chroms = [str(i) for i in range(1,23)] + ['X', 'Y']
    
    for line in F:
        
        if progress_level == progress_increment:
            d = datetime.datetime.now()
            ts = d.strftime('%D- %H:%M')
            print "processed {} variants {}".format(count, ts)
            progress_level = 0
        
        
        line = line.rstrip()
                         
        lin_spl = line.split('\t')
                         
        if count < header_end-1:
            annotated_vcf.write(line + '\n')
            
                         
        if count == header_end-1:
            # add filter information into the header
            lines = (['##FILTER=<ID=ac0,Description="no alternative alleles genotyped at this site">'] + 
         ['##FILTER=<ID=OTHER_CONTIG,Description="not on autosomes or sex chroms">'] +
         ['##FILTER=<ID=SHORT_LENGTH,Description="variant smaller than 50bp">'] + 
         ['##FILTER=<ID=SHORT_DEL,Description="del less than 418bp without SR">'] +
         ['##FILTER=<ID=MSQ,Description="MSQ score below 90 for BND or INV, below 100 for DUP, below 20 for DEL and MEI">'] +
         ['##FILTER=<ID=NMissing,Description="at least 10% of samples have a missing GT (iPSCORE fib/blood + HipSci fib)>'] + 
         ['##FILTER=<ID=PE/SR,Description="low PE/SR ratio for BND or INV">'] +
                     ['##FILTER=<ID=QUAL,Description="quality score below 100 for BND or INV">'] + 
                     ['##FILTER=<ID=iPSC,Description="present only in iPSC">'] )
            
            annotated_vcf.write('\n'.join(lines)+ '\n')
            annotated_vcf.write(line + '\n')
            
        
            header = lin_spl
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            
            samples = header[9:]
            info_cols = header[:9]
            cols_maf = ['NNREF', 'NNREF_UUIDs','ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'UUIDS_REF', 'NMissing']
            cols_extracted = ['SVTYPE','MSQ', 'SU', 'PE', 'SR', 'SVLEN', 'AF', 'NSAMP', 'CIPOS', 'CIEND', 'END'] 
            cols_maf_unr = [i+ '_unr' for i in cols_maf]
            cols_maf_non_ipsc = [i+ '_non_ipsc' for i in cols_maf]
            
            
            #remove the info column for space savings
            
            
            info_header = ['NAME'] +info_cols[:-2] + cols_extracted + ['Missing_GTs', 'Missing_AB'] + cols_maf + cols_maf_unr + cols_maf_non_ipsc
            
#             info_header.remove('INFO')
            
            genotypes_header = ['NAME'] + samples
            
            gt_file.write("\t".join(genotypes_header) + '\n')
            gt_filt_file.write("\t".join(genotypes_header) + '\n')
            
            ab_file.write("\t".join(genotypes_header) + '\n')
            ab_filt_file.write("\t".join(genotypes_header) + '\n')
            
            info_file.write("\t".join(info_header) + '\n')
            info_filt_file.write("\t".join(info_header) + '\n')
            
            

            
            
        elif count > header_end:
            
            if count == header_end+1:
            
                format_fields = lin_spl[header_dict['FORMAT']].split(':')
                format_dict = {l:i for l,i in zip(format_fields, range(0, len(format_fields)))}
                  
            
            sample_gts = []
            sample_ab = []
            # original first 9 columns
            
            info_col = header_dict['INFO']
            
            
         
            
            POS = lin_spl[header_dict['POS']]
            
            ID = lin_spl[header_dict['ID']]
            
            info = lin_spl[header_dict['INFO']]
            chrom = lin_spl[header_dict['#CHROM']]
            
      
            
            cols_to_parse = ['SVTYPE','MSQ', 'SU', 'PE', 'SR', 'SVLEN', 'AF', 'NSAMP', 'CIPOS', 'CIEND']
            
            cols_dict = dict(zip(cols_to_parse, range(0, len(cols_to_parse))))
            
            info_str = lin_spl[header_dict['INFO']]
            
            lin_spl[header_dict['INFO']] = info_str.replace('MEI', 'rMEI') 
            info = lin_spl[header_dict['INFO']]
            
            info_out = []
            for l in cols_to_parse:
                c = parse_info_col(lin_spl[header_dict['INFO']], l)
                
                info_out.append(c)
            
            
            try:
                MSQ = float(info_out[cols_dict['MSQ']])
               
            except:
                MSQ = info_out[cols_dict['MSQ']]
                
            
            
            SVTYPE = info_out[cols_dict['SVTYPE']]
            
            # correct naming for downstream
#             if SVTYPE == 'MEI':
#                 info_out[cols_dict['SVTYPE']] = 'rMEI'
#                 #correct the SVTYPE in vcf out_line- make it easier to combine all the vcfs if needed
                
#                 lin_spl[header_dict['INFO']] = lin_spl[header_dict['INFO']].replace('MEI', 'rMEI') 
#                 print l
            
            
            SU = int(info_out[cols_dict['SU']])
            SR = int(info_out[cols_dict['SR']])
            PE = int(info_out[cols_dict['PE']])
            QUAL = float(lin_spl[header_dict['QUAL']])
            SVLEN = info_out[cols_dict['SVLEN']]
            NSAMP = info_out[cols_dict['NSAMP']]
            
            
   
                
          
        
            # give the variants more useful IDs 
            END = int(POS) + 1
            
            if SVTYPE != 'BND':
                END = parse_info_col(info, 'END')
                
                ID_MOD = "{}_{}_{}_{}".format(SVTYPE, chrom, POS, END)
            else:
                # these need to be different, because they have no end, and will be non-unique
            
                ID_MOD = 'BND_{}_{}_{}'.format(chrom, POS, ID)
            
            lin_spl[header_dict['ID']] = ID_MOD
            
            # add minor allele frequency:
            
            gts_out = []
            ab_out = []
            
            # all samples in this case
            
            for s in samples:
                d = lin_spl[header_dict[s]].split(':')
                gt = d[format_dict['GT']]
                ab = d[format_dict['AB']]
                gts_out.append(gt)
                ab_out.append(ab)
            
            
            missing_gts = ('./.' in gts_out)
            
            missing_ab = ('.' in ab_out)
            
 
            
            # minor allele frequencies
            maf_data_all = calculate_minor_allele_freq(gts_out, samples, chrom, sex_dict)
        
            samples_unr_gts = [gts_out[samples.index(i)] for i in unrelated_samples]
            maf_data_unr= calculate_minor_allele_freq(samples_unr_gts, unrelated_samples, chrom, sex_dict)
            
            s = [gts_out[samples.index(i)] for i in non_ipsc_samples] 
            maf_data_non_ipsc = calculate_minor_allele_freq(s, non_ipsc_samples, chrom, sex_dict)
            
            # number of non ref samples in non-ipsc samples 
            non_ipsc_nnref = int(maf_data_non_ipsc[0])
            num_missing_non_ipsc = int(maf_data_non_ipsc[8])
            
            
            
            
            # Variant Filtering - Annotate the FILTER column with detailed filters
            
            FILTER = filter_variants(NSAMP, MSQ, QUAL, SVTYPE, PE, SR, SVLEN, SU, chrom,chroms, non_ipsc_nnref, num_missing_non_ipsc)
            msq_filt = MSQ_filter_specific_classes(NSAMP, MSQ, QUAL, SVTYPE, PE, SR, SVLEN, SU, chrom,chroms, non_ipsc_nnref)
            if msq_filt:
                if FILTER != 'PASS':
                    FILTER += ";MSQ"
                else:
                    FILTER = 'MSQ'
            
 
            # replace the FILTER with my FILTER column
            lin_spl[header_dict['FILTER']] = FILTER
            
            
        
            
            
            
            # write out data
            annotated_vcf.write("\t".join(lin_spl) + '\n')
            
            
            first_cols = lin_spl[:7]
            
            out_line_info =  [ID_MOD] + first_cols + info_out + [END] + [missing_gts] + [missing_ab] + maf_data_all + maf_data_unr + maf_data_non_ipsc
            
            out_line_info = map(str, out_line_info)
            info_file.write("\t".join(out_line_info) + '\n')
            
            
            out_line_gt = [ID_MOD] + gts_out
            gt_file.write("\t".join(out_line_gt) + '\n')
            
            
            out_line_ab = [ID_MOD] + ab_out
            ab_file.write("\t".join(out_line_ab) + '\n')
            
            
            
            # check specific filter tags
            def check_filters(x):
                filters = ['PE/SR', 'iPSC', 'ac0','OTHER_CONTIG', 'SHORT_LENGTH', 'SHORT_DEL']
                for i in filters:
                    if x.find(i) != -1:
                        return True
                return False
            
            
            # check for specific flat filters (no MSQ)
            if not check_filters(FILTER):
                info_filt_file.write("\t".join(out_line_info) + '\n')
                gt_filt_file.write("\t".join(out_line_gt) + '\n')
                ab_filt_file.write("\t".join(out_line_ab) + '\n')
                
            
                
            
#             if FILTER == 'PASS':
            
#                 info_filt_file.write("\t".join(out_line_info) + '\n')
            
#                 gt_filt_file.write("\t".join(out_line_gt) + '\n')

#                 ab_filt_file.write("\t".join(out_line_ab) + '\n')
                      
        
        count +=1
        progress_level += 1  
    
    # close out files
    F.close()
    gt_file.close()
    ab_file.close()
    info_file.close()
    annotated_vcf.close()
    gt_filt_file.close()
    ab_filt_file.close()
    info_filt_file.close()
    return fn_annotated_vcf


# In[3]:

def MSQ_filter_specific_classes(NSAMP, MSQ, QUAL, SVTYPE, PE, SR, SVLEN, SU, chrom, chroms, nnref_non_ipsc, msq_inv = 90, msq_mei = 20, msq_dup = 100, msq_del = 20, msq_bnd = 90):
    
    if type(MSQ) == float:
        if SVTYPE == 'BND':
            if MSQ < msq_bnd:
                return True
        elif SVTYPE == 'DUP':
            if MSQ < msq_dup:
                return True
            
        elif SVTYPE == 'DEL':
            if MSQ < msq_del:
                return True
        
        elif SVTYPE == 'rMEI':
            if MSQ < msq_mei:
                return True
        
        elif SVTYPE == 'INV':
            if MSQ < msq_inv:
                return True
            
        else:
            return False
    else:
        # MSQ is missing
        return True
    
    


# In[14]:

def filter_variants(NSAMP, MSQ, QUAL, SVTYPE, PE, SR, SVLEN, SU, chrom, chroms, nnref_non_ipsc, num_missing_non_ipsc):
        FILTER = []
        
        if int(NSAMP) == 0:
            FILTER.append('ac0')
        

        if ((num_missing_non_ipsc/478) > 0.1):
            FILTER.append('NMissing')

    
#         if SU < 8:
#             FILTER.append('SU<8')
            
        if chrom not in chroms:
            FILTER.append('OTHER_CONTIG')

        if SVTYPE != 'BND':
            SVLEN = abs(int(SVLEN))
            if SVLEN < 50:
                FILTER.append('SHORT_LENGTH')

            if (SVTYPE == 'DEL') & (SVLEN < 418) & (int(SR) == 0):
                FILTER.append('SHORT_DEL')
            if SVTYPE == 'INV':
                # require qual of 100 and ratio of PE/SU at least 10%
                if (QUAL < 100) | (QUAL == '.'):
                    FILTER.append('QUAL')

                if filter_PE_SR_Ratio(PE, SR, recip_ratio=0.1):
                    FILTER.append('PE/SR')



        if SVTYPE == 'BND':
            if QUAL == '.':
                FILTER.append('QUAL')
            # require qual of 100 and ratio of PE/SU at least 25%
            if QUAL < 100:
                FILTER.append('QUAL')

            if filter_PE_SR_Ratio(PE, SR, recip_ratio=0.25):
                FILTER.append('PE/SR')
                
        if nnref_non_ipsc == 0:
            FILTER.append('iPSC')


        if len(FILTER) == 0:
            FILTER = ['PASS']

        return ";".join(FILTER)
            


# In[181]:

def filter_PE_SR_Ratio(PE, SR, recip_ratio = 0.25):
    """ apply filtering similar to GTeX on these 
    variant types (making sure the required amount of evidence is present)"""

    tot = PE + SR
    for s in [PE,SR]:
        frac = s/tot
        if frac < recip_ratio:
            return True
    return False


# In[33]:

def safe_div(x, y, alt=0):
    try:
        return x/y
    except:
        return alt


# In[214]:

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

        if c == '0/0':
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
    
    out_names = ['NNREF', 'NNREF_UUIDs','ALTAF', 'REFAF', 'MAF', 'Minor_Allele', 'NREF', 'UUIDS_REF', 'Missing']
    out_data = [non_ref, ",".join(non_ref_uuids), alt_af, ref_af, maf, minor_allele, ref, ",".join(ref_uuids), missing]
    
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


# In[12]:

def run_from_args(args):
        
    unrelated_fn = args.unr_samples
    non_ipsc_fn = args.non_ipsc_samples

    unr_samples = [line.rstrip() for line in open(unrelated_fn)]
    non_ipsc_samples = [line.rstrip() for line in open(non_ipsc_fn)]
    
    
    gender_fn = args.gender_map
    sex_dict = {line.rstrip().split()[0]:line.rstrip().split()[1] for line in open(gender_fn)}
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')

    print "extracting variants, annotating filters and MAFs, start {}".format(ts)
    
    fn_annotated_vcf = add_annotations_extract_gts(args.vcf_file, args.output_dir, unr_samples, non_ipsc_samples, sex_dict)
    
    d = datetime.datetime.now()
    ts = d.strftime('%D- %H:%M')
    
    print 'complete', ts
    print 'annotated vcf written'
    print fn_annotated_vcf
    


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    
    sys.exit(args.entry_point(args))


# In[11]:

import datetime

