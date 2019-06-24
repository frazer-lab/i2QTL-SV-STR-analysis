
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

import copy 

import itertools
import tempfile
import six
import networkx as nx
import scipy.stats as stats

import datetime
import argparse

import warnings
warnings.filterwarnings("ignore")


# In[1]:

def rr_frame_collapsed_stats(df_in):
    df = df_in[df_in.replication_rate <> 'None']
    df['replication_rate'] = df.replication_rate.astype(float)
    tdf = df.groupby(['identifier','SVTYPE']).replication_rate.mean().to_frame()
    tdf2 = df.groupby(['identifier','SVTYPE']).size().to_frame()
    tdf2.rename(columns={0:'num_unique_sites'}, inplace=True)
    tdf = tdf.join(tdf2)
    tdf = tdf.reset_index()
    return tdf
    


# In[3]:

def replication_lm(fn, info, twin_list_uuid, identifier_col = 'None', lq_adjust = False):
    
    '''read in the genotype matrix file and parse it line by line for twin concordance stats''' 
    
    flattened_twins = list(set([i for sublist in twin_list_uuid for i in sublist]))
    
    
    data = []
    cols_dict = {}
    cols_twins_flattened = []
    
    cnv_classes = info['cnv_class'].to_dict()
#     print cnv_classes.keys()[0:2]
    
    
    starts = info['POS'].to_dict()

#         print gs_info[lq_adjust].head().to_dict()
 
    
    # collect data per twin pair per site- output a matrix for each index
    per_pair_data = {"_".join(tp): [] for tp in twin_list_uuid}
    indexes = []
    
    with open(fn, 'rU') as F:
        count = 0
        
        for line in F:            
            line = line.rstrip()
            lin_spl = line.split()
            
            if count == 0: 
                top_line = copy.deepcopy(lin_spl)
                cols_dict = {i:lin_spl.index(i) for i in flattened_twins}
    
            
            tps_discordant = []
            tps_concordant = []
            if count == 1:
                # the index doesn't have a column header- increment the indexes by one to fix
              
                if len(lin_spl) > len(top_line):
    
                    cols_dict = {i:k+1 for i,k in cols_dict.iteritems()}
                
            
            if count > 0:
                
                name_in = lin_spl[0]
                spl_name = name_in.split('_')
#                 print cnv_classes.keys()[0]
                sv_type = cnv_classes[name_in]
#                 print cnv_classes.keys()[0]
       

          
                
                chrom = spl_name[1]
                
                if sv_type != 'BND':
                    start = int(spl_name[2])
                    end =  int(spl_name[3])
                    Length = end - start

                else:
                    start = int(starts[name_in])
                    end =  start + 1
                    Length = end - start
                    
                    
                    
                number_with_var = 0
                number_concordant_with_var = 0
                number_discordant_with_var = 0
    
               
                
                if chrom not in ['X', 'Y']:
                    indexes.append(name_in)
                    for tp in twin_list_uuid:
                            
                        tp_id = "_".join(tp)
                        
                        twins = [cols_dict[tp[0]], cols_dict[tp[1]]]
                 
                        try:
                            vals = [lin_spl[i] for i in twins]
                        except:
                            
                            print name_in, tp, [lin_spl[i] for i in twins]
                            return 
                        
                        # check one way- first twin- for variant
                        if vals[0] not in ['0/0', './.']:
                        
                            number_with_var +=1
                            if vals[0] == vals[1]:
                                tps_concordant.append("_".join(tp))
                                number_concordant_with_var +=1
                                per_pair_data[tp_id].append('CONCORDANT')
                            else:
                                tps_discordant.append("_".join(tp))
                                number_discordant_with_var +=1 
                                per_pair_data[tp_id].append('DISCORDANT')
                        else:
                            per_pair_data[tp_id].append('REF')
                            
                            
                    if number_with_var > 0:
                        # replication rate- ranges from 0-100% 100% means every variant is a match
                        replication_rate = 1 - (number_discordant_with_var/number_with_var)
                    else:
                        replication_rate = 'None'
                    
                    if lq_adjust: 
                        out_line = [name_in, chrom, start, end, Length, number_with_var, number_concordant_with_var, 
                                number_discordant_with_var, replication_rate, sv_type, tps_discordant, tps_concordant,
                               num_lq_pairs, lq_samps_at_site, lq_pairs]
                    else:
                        out_line = [name_in, chrom, start, end, Length, number_with_var, number_concordant_with_var, 
                                number_discordant_with_var, replication_rate, sv_type, tps_discordant, tps_concordant]
                    
                    data.append(out_line)
              
            count +=1     
            
    
    if lq_adjust:
        
        out_df = pd.DataFrame(data, columns=['name','Chr','Start' ,'End', 'Length',
                                         'number_pairs_with_var', 'number_concordant_with_var',
                                         'number_discordant_with_var', 'replication_rate', 'SVTYPE', 
                                             'pairs_disc', 'pairs_conc', 'number_lq_pairs', 'lq_samps', 'lq_pairs'])
        
    else:
        out_df = pd.DataFrame(data, columns=['name','Chr','Start' ,'End', 'Length',
                                         'number_pairs_with_var', 'number_concordant_with_var',
                                         'number_discordant_with_var', 'replication_rate', 'SVTYPE', 
                                             'pairs_disc', 'pairs_conc'])
    per_twin_df = pd.DataFrame(per_pair_data)
    per_twin_df['ID'] = indexes
    per_twin_df.index = per_twin_df.ID
    
    out_df['identifier'] = identifier_col
    disc_collapsed = rr_frame_collapsed_stats(out_df)
    out_df.index = out_df['name']
    out_df.index.name = 'index'
    return out_df, disc_collapsed, per_twin_df


# In[4]:

def generate_per_pair_summary(df):
    df = df.copy()
    if 'ID' in df.columns:
        df = df.drop('ID', axis = 1)
        
    df = df.apply(pd.value_counts).T.stack().to_frame('counts').reset_index()
    names = ['pair', 'category', 'counts']
    df.columns = names
    df = df.pivot_table(values = 'counts', index = 'pair', columns = ['category'])
    
    for c in df.columns:
        df[c] = df[c].fillna(0).astype(int)
    
    df['total_variants'] = df.CONCORDANT + df.DISCORDANT
    df['percent_concordant'] = df.CONCORDANT/ df.total_variants
    return df


# In[2]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-gt", "--gt", dest="gt_tsv", metavar='<gt>', help="genome strip cns tsv", required=True)
    
    
    parser.add_argument("-info", "--info", dest="info_pkl", metavar='<info>', help="genome strip info pickle", required=True)
    
    parser.add_argument("-pairs", "--pairs", dest="pairs", metavar='<pair_fn>', help="file of pairs of samples, likely genetically duplicate samples (twins, replicates, biological replicates) (tsv with no header)", required=True)
    
    parser.add_argument("-id", "--id", dest="id", metavar='<id_col>', help="name for identifier column if desired (useful if several data sets are combined from different discovery efforts)", required=True, default = 'None')
    
    parser.add_argument("-suff_df", "--suffix_dataframe", dest="suff_df", metavar='<suffix for df columns>', help=" suffix for the df column", required=False, default = False)
    
    
    
    
    parser.add_argument("-o", "--output_dir", dest="output_dir", metavar='<out_dir>', help="output directory for summary output", required=True)
    

    
    parser.add_argument("-pre", "--prefix", dest="prefix", metavar='<prefix>', help="prefix to name files", default = False)
    
    parser.add_argument("-suff", "--suffix", dest="suffix", metavar='<suffix>', help="prefix to name files", default = False)
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to calculate and annotate the info data frame with discordance rate in twins')
    add_arguments_to_parser(parser)
    return parser


# In[5]:

def run_from_args(args):
    gt_fn = args.gt_tsv
    info = pd.read_table(args.info_pkl, index_col=0)
    
    
    pair_fn = args.pairs
    id_col = args.id
    
    
    pairs = [line.rstrip().split() for line in open(pair_fn)]
    
    
    
    print 'calculating replication rate statistics', CM.datestring(hour=True, minute=True)
    
    # this is generally the correct column name in VCF files so I'll add it in and use it everywhere
    
    try:
        info['SVTYPE'] = info['cnv_class']
        
    except:
        info['cnv_class'] = info.SVTYPE
        
        
    rr_df, collapsed_stats_df, per_pair_df = replication_lm(gt_fn, info, pairs, identifier_col=id_col)
    per_pair_summary = generate_per_pair_summary(per_pair_df)
    
    #cols = ['number_twins_with_var', 'number_concordant_with_var','number_discordant_with_var', 'discordance_score', 'identifier']
    
    
    if args.suff_df:
        suff_df = str(args.suff_df)
        rename = {i:i + suff_df for i in cols}
        cols = [i + suff_df for i in cols]
        
        # might need to have some suffixes if using this script on multiple sets of data
        rr_df.rename(columns=rename, inplace=True)
    
    #gs_info = gs_info.join(disc_frame[cols])
    
    output_location = args.output_dir
    
    if args.suffix:
        fn_info = os.path.join(output_location, 'gs_info' + args.suffix)
        fn_cns = os.path.join(output_location, 'gs_cns' + args.suffix)
        
        var_name_info = 'rr_per_site' + args.suffix
        var_name_pair_rr_summary = 'rr_summary' + args.suffix
        var_name_per_pair = 'rr_per_pair' + args.suffix
        var_name_per_pair_summary = 'rr_per_pair_summary' + args.suffix
    
    else:
        fn_info = os.path.join(output_location, 'gs_info')
        fn_cns = os.path.join(output_location, 'gs_cns')
        var_name_info = 'rr_per_site'
        var_name_pair_rr_summary = 'rr_summary'
        var_name_per_pair = 'rr_per_pair'
        var_name_per_pair_summary = 'rr_per_pair_summary'
        
    
    
    print 'calculation complete',  CM.datestring(hour=True, minute=True)
    
    CM.save_dataframe(var_name_info, rr_df, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_pair_rr_summary, collapsed_stats_df, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_per_pair, per_pair_df, output_location, print_vars_recorded_loc=False)
    CM.save_dataframe(var_name_per_pair_summary, per_pair_summary, output_location, print_vars_recorded_loc=False)
    
    


# In[5]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

