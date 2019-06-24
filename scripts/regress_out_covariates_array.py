#!/usr/bin/env python
# coding: utf-8

# In[2]:


from __future__ import division
import numpy as np
import os
import sys
import datetime
import pandas as pd


import copy
# import ciepy
# import cardipspy as cpy


# dy_name = 'intensity_array_processing'

# private_out = os.path.join(DJ.root, 'private_output', dy_name)
# if not os.path.exists(private_out):
#     cpy.makedir(private_out)

import statsmodels.api as sm 


# In[3]:


predictors = pd.read_pickle('/frazer01/projects/CARDIPS/analysis/cardips-cnv-analysis/private_output/intensity_array_processing/predictors.pkl')


# In[8]:


cols = predictors.columns.tolist()


# In[4]:


samples = predictors.index.tolist()


# In[4]:


fn = '/frazer01/projects/CARDIPS/data/rsng/160823/extracted_arrays/combined_array_filt.txt'
# df = pd.read_table(fn)


# In[22]:


with open(fn, 'rU') as F:
    count = 0
    for line in F:
 
        line = line.rstrip()
        lin_spl = line.split(",")

        if count == 0:
            header = lin_spl
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            sample_order = lin_spl[5:]
            # ensure sample order
            predictors = predictors.loc[sample_order]
            covar = predictors.values
            # print out the header line
            print "\t".join(lin_spl[1:])

        else:
          
            data = lin_spl[5:]
            id_cols = lin_spl[1:5]
            chrom = lin_spl[header_dict['CHROM']]
            
            dtypes = [str, str, float, float]
            dtypes2 = [str, str, int, int]
            
            id_cols = [t(i) for t, i in zip(dtypes, id_cols)]
            id_cols = [t(i) for t, i in zip(dtypes2, id_cols)]
            id_cols = [str(i) for i in id_cols]
            
             
            try:
                chrom = float(chrom)
                chrom = str(int(chrom))
            except:
                pass

            id_cols[1] = chrom
            
            
            data = np.array([float(i) for i in data])
            l = sm.regression.linear_model.OLS(data, covar)
            f = l.fit()


            # regress out covars 
            data_original = copy.deepcopy(data)
            
            # remove all covariates except for the intercept
            for coef in f.params[:-1]:
                data = data - (data_original * coef)

            data = list(data)
            data = [str(round(i, 4)) for i in data]
            # print out the corrected data
            
            print "\t".join(id_cols + data)
        
        count +=1

