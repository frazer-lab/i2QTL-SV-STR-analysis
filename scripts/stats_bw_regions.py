
# coding: utf-8

# In[14]:

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
import argparse

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
from scipy.stats import mode
pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 100)
from mpl_toolkits.axes_grid1 import  make_axes_locatable
import datetime

from scipy.stats import mode

import scipy.stats as stats


from djPyBio import Common as CM


# In[1]:

def lambda_mod_hist(x):
    # fix the hist data, make a col of lists of floats for stats
    x = x.split(',')
    x = np.where(np.array(x)=='NA',0.0,x)
    x = [float(i) for i in x] 
    return x


# In[13]:

def add_summary_stats(fn, uuid):
    dtype_dict = dict(zip(range(0,5), [object,int, int, object, int, object]))
    df = pd.read_table(fn, names=['Chr', 'Start', 'End', 'ID','Length', 'Histogram'], dtype= dtype_dict)
    df['Histogram_Mod'] = df.Histogram.apply(lambda x: lambda_mod_hist(x))
    funcs = [np.mean, np.max, np.median, stats.mstats.mode]
    col_names = ['mean', 'max', 'median', 'mode']
    
    for f, cn in zip(funcs, col_names):
        df[cn] = df.Histogram_Mod.apply(f)
        
    df['mode'] = df['mode'].apply(lambda x: x[0][0])
    df['uuid'] = uuid
    
    df.drop(["Histogram", "Histogram_Mod"], axis =1, inplace = True)
    
    return df


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'generate stats for bw regions for SVs')
    add_arguments_to_parser(parser)
    return parser


# In[15]:

def add_arguments_to_parser(parser):
    parser.add_argument("-fn", "--fn", dest="fn_in", metavar='<fn_bw_regions.txt>', help="tab delimited bedgraph regions for SVs", required=True)
    parser.add_argument("-SN", "--sample", dest="sample", metavar='<sample_name>', help="sample_name", required=True)
    parser.add_argument("-o", "--output", dest="output", metavar='<output_fn>', help="output_fn", required=True)
    
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def run_from_args(args):
   
    df = add_summary_stats(args.fn_in, args.sample)
    df.to_csv(args.output, sep= '\t')


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

