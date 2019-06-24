
# coding: utf-8

# In[3]:

from __future__ import division
import numpy as np
import os
import sys
import datetime
from subprocess import call
import subprocess
import glob
import djPyBio as DJ
import pandas as pd
# import csv
# import seaborn as sns
# import matplotlib as mpl
# import matplotlib.pyplot as plt
import copy 
# import pybedtools as pbt
import ciepy
import cardipspy as cpy
import itertools
import tempfile
import gzip
import six
import networkx as nx
import scipy.stats as stats
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse


from djPyBio import Common as CM


# In[9]:

def find_header_end(fn):
    with gzip.open(fn) as F:
        count = 0
        header = False
        count_out = 0
        for line in F:
            count +=1
            try: 
                spl = line.split('\t')
                spl0 = spl[0] 
                if spl[0]=="#CHROM":
                    header = True
                    count_out  = count
                    
                    return count_out, line.split()
                if count > 250:
                    print count
                    break

                    
            except:
                if count < 10:
                    continue
                else:
                    break
    
def length_classifier(REF_len, ALT_len_max, ALT_len_min, variant_class):

    out_ins = ALT_len_max - REF_len
    out_del = REF_len - ALT_len_min
    
    if variant_class == 'INS':
        return out_ins
    elif variant_class == 'DEL':
        return out_del
    
    else:
        return REF_len



def concordance_lambda(x, uuid1,uuid2):
    gts = [x[uuid1], x[uuid2]]
    
    if x['variant_class']=='SNP':
        if x[uuid1]==x[uuid2]:
            return 'concordant'

        elif (['0/1','1/1'] == gts) or  (['1/1','0/1'] == gts):
            return 'het_homo'
        
        else:
            return 'discordant'
    else:
        if x[uuid1]==x[uuid2]:
            return 'concordant'
        else:
            return 'discordant'


# In[61]:


def column_len_counts(x):
    data = x
    spl = data.split(',')
    
    lengths = []
    for i in spl:
        lengths.append(int(len(i)))

    
    return lengths
    


# In[62]:

def indel_classifier_2(REF, ALT):
    "given a list of reference and alternate lengths, classify each alternative allele as a SNV or indel"
    
    
    length_REF = column_len_counts(REF)[0]
    lengths_ALT = column_len_counts(ALT)
    
    
    types_out = []
    lengths_out = []
    for i in lengths_ALT:
        if length_REF == i:
            types_out.append('SNV')
            lengths_out.append(0)
        
        elif length_REF < i :
            types_out.append('INS')
            length = i - length_REF
            lengths_out.append(length)
        
        elif length_REF > i :
            types_out.append('DEL')
            length = length_REF - i
            lengths_out.append(length)
        
        else:
            return 'broken case'
        
    return types_out, lengths_out
        
    


# In[ ]:

def flatten_list(x):
    merged = list(itertools.chain.from_iterable(x))
    return merged


# In[60]:

def parse_info_col(t, lab):
    t = t.split(';')
    
    cols_cleaned = []
    for l in t:
        if len(l.split('=')) == 1:
            pass
        else:
            cols_cleaned.append(l)
            
   

    cols = [i.split('=')[0] for i in cols_cleaned]
    
    vals = [i.split('=')[1] for i in cols_cleaned]    
    
    try:
        ind = cols.index(lab)
        v = vals[ind].split(',')
        v = [float(i) for i in v]
#         out = sum(v)
        return v
    except:
#         print t
#         print cols
#         print vals
        
        print 'irregular parse'
        return ['None']


# In[63]:

def clean_hist_dfs(df):
    
    # convert the hist df into a plottable format for 

    df.fillna(0, inplace=True)

    df.index = df.variant_types

    df = df.T.iloc[:-1,:]

    df =df.reset_index()

    df.rename(columns={'index':'MAF'}, inplace=True)
    
    
#     df = df.groupby('MAF').sum()
#     df = df.reset_index()
    return df


# In[94]:

def get_chrom(fn , chrom):
    bcf_tools_path = '/frazer01/software/bcftools-1.2/bcftools'
    command = '{} view -r {} -H {}'.format(bcf_tools_path, chrom, fn)
    
    output = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    stdout =  output.stdout
    
    if stdout.readline() == '':
        print 'No results for this chromosome'
        return False
    else:
        return stdout


# In[5]:

def calculate_MAF_hist(fn, samples, unrel_samples, Chroms, bed_exclude = False, prefix = False, maf_col = 'MAF_UNREL', add_unrel_to_info = False):
    
    
    
    header_end,header_line = find_header_end(fn)
   
    cols_dict = {i:header_line.index(i) for i in header_line}
  
    print "Starting MAF Histogram Extraction..."
    print "{} variants_processed {}".format(0, CM.datestring(hour=True, minute=True))
    variants_processed = 1000000
    
    if bed_exclude:
        command = ['bedtools','intersect','-a', fn, '-b', bed_exclude, '-v', '-header']
        F = subprocess.Popen(command, stdout=subprocess.PIPE)
        F = F.stdout
    else:
        F = gzip.open(fn)
        
    
    var_classes = ['SNV', 'INS', 'DEL']
    hist_dicts = [{},{},{}]
    
    var_classes_simple = ['SNV', 'INDEL']
    hist_dicts_simple = [{},{}]
    
    hist_dicts_var_lengths = [{}, {}, {}]
    hist_dicts_var_lengths_simple = [{}, {}]
    
    count = 0
    
    nr_var_dict = [{},{},{}]
    nr_var_dict_simple = [{},{}]

    for line in F:
        conc = 0
        disc = 0
        non_ref_conc = 0
        non_ref_disc = 0
        count +=1
        
        if count == header_end:
            
            # add shit to header to print
            
            if add_unrel_to_info:
                print '##INFO=<ID=NREF_UNREL,Number=A,Type=Integer,Description="number of non-reference samples per allele in unrelated individuals">'
                print '##INFO=<ID=NREF,Number=A,Type=Integer,Description="number of non-reference samples per allele in all samples">'
        if count > header_end:
        
            
            
            var_num = count - header_end
            
            nr_at_site_unrel = []
            nr_at_site_all = []
            
            
            line = line.rstrip()
            lin_spl = line.split()
            FORMAT = lin_spl[cols_dict['FORMAT']]
            INFO = lin_spl[cols_dict['INFO']]
            REF = lin_spl[cols_dict['REF']]
            ALT = lin_spl[cols_dict['ALT']]
            CHROM  = lin_spl[cols_dict['#CHROM']]
            ID = lin_spl[cols_dict['ID']]
            POS= lin_spl[cols_dict['POS']]
            
            if var_num == variants_processed:
                msg = "{} variants_processed {}".format(variants_processed, CM.datestring(hour=True, minute=True))
                print msg

                
                variants_processed += 1000000
                
            
            if CHROM in Chroms:

                # Classify indels and snvs based on the alternate and reference allele column

                type_convert = {'SNV':'SNV', 'INS':'INDEL', 'DEL':'INDEL'}
                
#                 ref_len = column_len_counts(REF)[0]
#                 alt_len_max, alt_len_min, var_lengths = column_len_counts(ALT)
                
                
                # fix this to split out the variants that are multi-allelic and count them all separately
                
                
                var_types = []
                
                
                var_types, var_lengths = indel_classifier_2(REF, ALT)
                
                mafs = parse_info_col(INFO, maf_col)
           
                
                allele_nums = range(1, len(mafs) + 1)
                
                for vt, length, maf, an in zip(var_types, var_lengths, mafs, allele_nums):
                    
                
                # add to hist dicts
                    ind = var_classes.index(vt)

                    hist_dicts[ind][maf] = hist_dicts[ind].get(maf, 0) +1
                    
                    hist_dicts_var_lengths[ind][length] = hist_dicts_var_lengths[ind].get(length, 0) +1


                    # convert to simple var class and add to other set of dicts
                    var_simple= type_convert[vt]
                    
                    ind = var_classes_simple.index(var_simple)
                    hist_dicts_simple[ind][maf] = hist_dicts_simple[ind].get(maf, 0) +1
                    hist_dicts_var_lengths_simple[ind][length] = hist_dicts_var_lengths_simple[ind].get(length, 0) +1
                    
                    
                    
         
                    nr_ur_samps = 0
                    nr_all_samps = 0
                    for samp in samples:
                        gt = lin_spl[cols_dict[samp]].split(':')[0]
                        if gt not in ['0/0', './.']:
                            alleles = [int(l) for l in gt.split('/')]
                            if alleles.count(an) > 0:
                                nr_all_samps +=1
                                if samp in unrel_samples:
                                    nr_ur_samps +=1


                    nr_at_site_unrel.append(nr_ur_samps)
                    nr_at_site_all.append(nr_all_samps)
                    
                    if nr_ur_samps == 0:
                        nr_at_site_unrel.append(0)
                        
                    if nr_all_samps == 0:
                        nr_at_site_all.append(0)
                        
                        
                        # mark samples that are singleton in unrelateds (in 1 of the selection of UR)
                        if nr_ur_samps == 1:
                            ind = var_classes.index(vt)                                     
                            nr_var_dict[ind][maf] = nr_var_dict[ind].get(maf, 0) +1
                            
                            ind = var_classes_simple.index(var_simple)
                            nr_var_dict_simple[ind][maf] = nr_var_dict_simple[ind].get(maf, 0) +1
                            
                    
                    
                    # add the nref_unrel_column if desired to vcf
                    # print the new line with that annotated onto it
                if add_unrel_to_info:
                    nref_formatted_unr = ";NREF_UNREL=" + ','.join([str(v) for v in nr_at_site_unrel])
                    nref_formatted_all = ";NREF=" + ','.join([str(v) for v in nr_at_site_all])
                    
                    lin_spl[cols_dict['INFO']] += nref_formatted_unr
                    lin_spl[cols_dict['INFO']] += nref_formatted_all
                    
                    print "\t".join(lin_spl)
                        
        else:
            if add_unrel_to_info:
                print line.rstrip()
                    
                    
                    
                    
                    
                    
                    
                
                
                
    # return as df for saving later
    
    df = pd.DataFrame(hist_dicts)
    df['variant_types'] = var_classes
    df = clean_hist_dfs(df)
    
    df2 =  pd.DataFrame(hist_dicts_simple)
    df2['variant_types'] = var_classes_simple
    df2 = clean_hist_dfs(df2)
    
    df3 =  pd.DataFrame(hist_dicts_var_lengths)
    df3['variant_types'] = var_classes
    df3 = clean_hist_dfs(df3)
    
    df4 =  pd.DataFrame(hist_dicts_var_lengths_simple)
    df4['variant_types'] = var_classes_simple
    df4 = clean_hist_dfs(df4)
    
    
    df5 =  pd.DataFrame(nr_var_dict)
    df5['variant_types'] = var_classes
    df5 = clean_hist_dfs(df5)
    
        
    df6 =  pd.DataFrame(nr_var_dict_simple)
    df6['variant_types'] = var_classes_simple
    df6 = clean_hist_dfs(df6)
    
    
    
    return df, df2, df3, df4, df5, df6


# In[3]:

len(['','',''])


# In[4]:

range(1,3)


# # Command Line Argument Handling

# In[ ]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-vcf", "--input", dest="vcf", metavar='<VCF>', help="vcf with snvs and indels", required=True)
    
    parser.add_argument("-Chrom", "--Chr", dest="chroms", metavar='<chromosomes>', help="chromosomes to include in calculation, comma separated list no spaces", default= '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y')
    
    parser.add_argument("-maf_col", "--maf_col", dest="maf_col", metavar='<maf_column_name>', help="name of the column with MAF that you want to get the histogram of", default= 'MAF_UNREL')
        
    parser.add_argument("-s", "--samples", dest="samples", metavar='<samples>', help="set of samples to annotate", required=True)
    
    parser.add_argument("-add_nref", "--add_nref", dest="add_nref", action = 'store_true', help="print out vcf, adding in changes to the info column with nref per unrelated", default = False)
    
    
    
    parser.add_argument("-ur", "--unrelated", dest="unrelated", metavar='<unrelated_samples>', help="unrelated samples to annotate", required=True)
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    parser.add_argument("-e", "--exclude_bed", dest="exclude_bed", metavar='<bed_file>', help="bedfile_to_exclude", default = False)
    
    parser.add_argument("-suff", "--suff", dest="suffix", metavar='<suffix>', help="suffix to name files", default = False)
    

    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'a command line utility to produce a MAF and lengths histogram at SNV and indel sites')
    add_arguments_to_parser(parser)
    return parser


# In[ ]:

def run_from_args(args):

    unrelated_file = args.unrelated
    unr_samples = [line.rstrip() for line in open(unrelated_file)]
    
    samples_file = args.samples
    samples = [line.rstrip() for line in open(samples_file)]
    
    
    
    
    # prepare chromosome list
    chroms = args.chroms
    chroms = chroms.split(',')
    
    df1, df2, df3, df4, df5, df6 = calculate_MAF_hist(args.vcf ,samples, unr_samples, 
                                                      chroms, bed_exclude = args.exclude_bed, 
                                                      maf_col = args.maf_col, add_unrel_to_info= args.add_nref)
    
    
    if args.suffix:
        fn_df1 = "snv_indel_maf_full" + args.suffix
        fn_df2 = "snv_indel_maf_simple" + args.suffix
        fn_df3 = "snv_indel_lengths_full" 
        fn_df4 = "snv_indel_lengths_simple"
        
        fn_df5 = "snv_indel_sing_maf_full" + args.suffix
        fn_df6 = "snv_indel_sing_maf_simple" + args.suffix 
        
        
        
        
    else:
  
        fn_df1 = "snv_indel_maf_full" 
        fn_df2 = "snv_indel_maf_simple"
    
    CM.save_dataframe(fn_df1, df1, args.output_dir)
    CM.save_dataframe(fn_df2, df2, args.output_dir)
    
    CM.save_dataframe(fn_df3, df3, args.output_dir)
    CM.save_dataframe(fn_df4, df4, args.output_dir)
    
    CM.save_dataframe(fn_df5, df5, args.output_dir)
    CM.save_dataframe(fn_df6, df6, args.output_dir)
   
    
    


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

