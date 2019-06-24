
# coding: utf-8

# In[ ]:

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
import csv
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy 
import pybedtools
pd.set_option('display.max_columns', 500)


def read_merged_lumpy_vcf(fn):
    """Read and parse a merged VCF of lumpy results."""
    df = pd.read_table(fn, skiprows=641)
    df.columns = [x.replace('#', '') for x in df.columns]
    df = df[df.ALT.apply(lambda x: x in ['<DEL>', '<DUP>'])]
    gcols = df.columns[8:]
    genotypes = df[gcols]
    genotypes = genotypes.apply(lambda x: [y.split(':')[0] for y in x], axis=1)
    df = df.drop(gcols, axis=1)
    cols = [x.split('=')[0] for x in df.INFO[0].split(';')]
    imprecise = []
    rows = []
    for i in df.index:
        vals = list(df.ix[i, 'INFO'].split(';'))
        if 'IMPRECISE' in vals:
            imprecise.append(True)
            vals.remove('IMPRECISE')
        else:
            imprecise.append(False)
        rows.append(dict([x.split('=') for x in vals]))
    df['imprecise'] = imprecise
    tdf = pd.DataFrame(rows, index=df.index)
    df = df.join(tdf)
    df = df.drop('INFO', axis=1)
    df.CHROM = 'chr' + df.CHROM.astype(str)
    # cols = df.FORMAT[0].split(':')
    # ds = df.apply(lambda x: pd.Series(dict(zip(x['FORMAT'].split(':'), x[df.columns[8]].split(':')))), axis=1)
    # ds = ds.drop(set(df.columns) & set(ds.columns), axis=1)
    # df = df.join(ds)
    # df = df.drop(['FORMAT', df.columns[8]], axis=1)
    df.ALT = df.ALT.apply(lambda x: x[1:4])
    df = df[df.END.isnull() == False]
    for c in ['POS', 'END', 'PE', 'SR', 'SU', 'SVLEN']:
        df[c] = df[c].astype(int)
    return df, genotypes


# In[1]:

def has_variant_Gt(row):
    Num_Var= 0
    range_col = range(0,len(row))
    
    for x in range_col:
        if row[x] =='0/1':
            Num_Var +=1
        if row[x]=='1/1':
            Num_Var+=1
    return Num_Var
        


# In[1]:

def get_genomestrip_vcf_copy_number_states(vcf):
    """Parse a genomestrip VCF file and get the diploid copy number states
    for each sample."""
    f = open(vcf)
    line = f.readline()
    while line[0:2] == '##':
        line = f.readline()
    header = line[1:].strip().split('\t')

    ind = []
    copy_numbers = []
    line = f.readline().strip()
    while line != '':
        t = line.split('\t')
        ind.append(t[2])
        copy_numbers.append([int(x.split(':')[1]) for x in t[9:]])
        line = f.readline().strip()
    cns = pd.DataFrame(copy_numbers, index=ind, columns=header[9:])
    return cns



def header_reorder_vcf(df, type_data = 'lumpy'):
    header_dict=False
    
    if type_data=='lumpy':
        try:
            desired_headers = ['CHROM', 'POS', 'END']
            head = df.columns.tolist()
            trunc_head = [x for x in head if x not in desired_headers]
            out_head = desired_headers + trunc_head
            ordered_df = df[out_head]
            header_dict = {'CHROM':'Chr', 'POS':'Start', 'END': 'End'}
            for key in header_dict:
                try:
                    ordered_df = ordered_df.rename(columns={key:header_dict[key]})
                except:
                    pass

            return ordered_df
        except:
            return df
                       
    elif type_data=='gs':
        try:
            desired_headers = ['chrom', 'start', 'end']
            head = df.columns.tolist()
            trunc_head = [x for x in head if x not in desired_headers]
            out_head = desired_headers + trunc_head
            ordered_df = df[out_head]
            header_dict = {'chrom': 'Chr', 'start':'Start', 'end': 'End'}
            for key in header_dict:
                try:
                    ordered_df = ordered_df.rename(columns={key:header_dict[key]})
                except:
                    pass
            return ordered_df
            
        except:
            return df

def dist_lambda(x,dict_in):
    Chr = x.Chr
    Start = x.Start
    End = x.End
    
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
    

def diff_sets(x, a, b):
    o = x[a]
    c = x[b]
    out = o.difference(c)
    return out


# Lambda for the deletion calls
def non_ref_lambda(x):
    NNon_Ref = 0 
    for l in x:
        if l =='0/1':
            NNon_Ref +=1
        elif l == '1/1':
            NNon_Ref +=1
    return NNon_Ref



def alleles_Discovered(x, list_cols):
    unique_gts = []
    for l in list_cols:
        unique_gts.append(x[l])
    
    out = set(unique_gts)
    
    return out

def variants_Discovered(x, list_cols):
    unique_gts = []
    for l in list_cols:
        if x[l] <> 2:
            unique_gts.append(x[l])
    
    out = set(unique_gts)
    
    return out

def alleles_CNDIST(x, list_cols):
    unique_gts = {}
    for l in list_cols:
        unique_gts[x[l]] = unique_gts.get(x[l], 0) +1
    
    return unique_gts

def is_mixed(x):
    duplication = False
    deletion = False
    for l in x: 
        if l<2:
            deletion=True
        elif l>2:
            duplication=True
    if deletion== True and duplication==False:
        return 'DEL'
    elif duplication==True and deletion==False:
        return 'DUP'
    elif duplication==True and deletion==True:
        return 'MIXED'

def calls_between_coords(x, start, end):
    
    """ 
    lambda function to subset calls between coordinates
    given that columns are named Start and End for 
    Chromosomal positions
    
    """
    if x.Start > start and x.Start < end:
        return True
    elif x.End < end and x.End > start:
        return True
    
    else:
        return False
    
def get_vars_overlapping(df, Chr, Start, End):
    
    out_df = df[df.Chr==Chr]
    sub = out_df.apply(lambda x: calls_between_coords(x, Start, End), axis=1)
    out_df = out_df.ix[sub[sub==True].index]
    return out_df

def nvariants_2(x):
    chrom = x.Chr
    start = x.Start
    end = x.End
    name = x['name']
    nvariants= 0
    for i in UUID_Dicts.keys():
        if chrom == 'X':
            if UUID_Sex_Dict[i]=='F':
                if x[i] >2 or x[i] < 2:
                    nvariants +=1
            if UUID_Sex_Dict[i]=='M':
                if x[i]<1 or x[i] > 1:
                    nvariants +=1
        if chrom == 'Y':
            if UUID_Sex_Dict[i]=='M':
                if x[i]<1 or x[i] > 1:
                    nvariants +=1
            if UUID_Sex_Dict[i]=='F':
                pass
        else: 
            if x[i] >2 or x[i] < 2:
                nvariants +=1
    return nvariants

