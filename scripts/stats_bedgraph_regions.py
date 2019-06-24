
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


# In[13]:

def calc_hist_stats(fn, uuid):
    cols = ['ChrA', 'StartA', 'EndA', 'ID',  'ChrB', 'StartB', 'EndB', 'rd']
    dtype_dict = dict(zip(range(0,8), [object,int, int, object, object, int, int, int]))

    
    df_in = pd.read_table(fn, names=cols, dtype = dtype_dict)
    df_out = df_in.groupby('ID').rd.agg({'mean': np.mean, 'median': np.median, 'mode': stats.mstats.mode, 
                                      'min_rd': min, 'max_rd': max})
    
    df_out['mode'] = df_out['mode'].apply(lambda x: x[0][0])
    df_out['uuid']=uuid
    
    return df_out


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'generate stats for bedgraph regions for SVs')
    add_arguments_to_parser(parser)
    return parser


# In[15]:

def add_arguments_to_parser(parser):
    parser.add_argument("-fn", "--fn", dest="fn_in", metavar='<fn_bedgraph_regions.txt>', help="tab delimited bedgraph regions for SVs", required=True)
    parser.add_argument("-SN", "--sample", dest="sample", metavar='<sample_name>', help="sample_name", required=True)
    parser.add_argument("-o", "--output", dest="output", metavar='<output_fn>', help="output_fn", required=True)
    
    
    parser.set_defaults(entry_point=run_from_args)


# In[ ]:

def run_from_args(args):
    df = calc_hist_stats(args.fn_in, args.sample)
    df.to_csv(args.output, sep= '\t')


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

