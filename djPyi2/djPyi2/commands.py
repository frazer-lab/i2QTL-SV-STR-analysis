
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
import pandas as pd
import csv
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy 
import pybedtools
import itertools
import glob
pd.set_option('display.max_columns', 500)


# In[1]:

def datestring(month=True, day=True, year=True, hour=False, minute=False):
    d = datetime
    
    y = str(d.datetime.now().timetuple().tm_year)
    m =  str(d.datetime.now().timetuple().tm_mon)
    D = str(d.datetime.now().timetuple().tm_mday)
    h = str(d.datetime.now().timetuple().tm_hour)
    mins = str(d.datetime.now().timetuple().tm_min)
    
    if int(month) in range(0,10):
        m = "0" + m
    
    
    dt_array = [y, m, D, h, mins]
    t_f_array = [year, month, day, hour, minute]
    
    Date_Array = []
    for s, i, in zip(dt_array, t_f_array):
        if i == True:
            Date_Array.append(s)
    
    return '_'.join(Date_Array)


# In[8]:


def qstat_u_df(user):
    """ Takes a qstat for a user, parses, and returns
    pandas dataframe"""
    
    command = ['qstat', '-u', user]
    out = subprocess.Popen(command, stdout=subprocess.PIPE)
    lines = []
    for line in out.stdout:
        lines.append(line.rstrip())
    head = lines[0]
    head = head.split(' ')
    head = filter(None, head)
    lines = lines[2:]
    for x,y in enumerate(lines):
 
        lines[x] = y.split(' ')
        lines[x] = filter(None, lines[x])
        if len(lines[x]) < len(head):
            lines[x].append(np.NaN)
    
    for x, y in enumerate(lines):
        state = y[4]
    
        if state == 'qw' or state =='hqw':
        
            out = y[0:7]
            slots= y[7:]
            lines[x] = out + ['Unknown'] + slots
        
    df = pd.DataFrame(lines, columns=head)
    
    return df


# In[9]:

# qstat_u_df('joreyna')


# In[1]:

def qstat_u_r_df(user):
    
    """ Returns Data_Frame with full info about job_names and resources
    needs to be developed for special categories, works with genome strip 
    jobs"""
    
    command = ['qstat', '-r', '-u', user]
    out = subprocess.Popen(command, stdout=subprocess.PIPE)
    lines = []
    for line in out.stdout:
        lines.append(line.rstrip())
        
    jobs = True
    try:
 
        head = lines[0]
        head = head.split(' ')
        head = filter(None, head)
        lines = lines[2:]
    except:
        jobs = False
        
        print 'No Jobs Running'
        
    
    if jobs==True:
        
    
        for x,y in enumerate(lines):
     
            lines[x] = y.split(' ')
            lines[x] = filter(None, lines[x])
        
        counter = 0
        data = []
        headers = ['job_ID', 'priority', 'user',
                  'status', 'date', 'time_start',
                  'queue', 'slots', 'job_name', 
                  'h_vmem']
        
        OUT = []
        
        
        # collapse the lines into their states
        
        queued = []
        running = []
        
        for x,y in enumerate(lines):
            line_pairs =[]
                 
            try:
                if lines[x][4] == 'r':
                    
    #  add in support for Requested PE/Granted PE here
                    running_line = True
                    running.append(list(itertools.chain(*lines[x:x+6])))
                elif lines[x][4]=='qw':
                    queued.append(list(itertools.chain(*lines[x:x+5])))
            except:
                pass
            
        for x in running:
            job_ID = x[0]
            prior = x[1]
            user = x[3]
            status = x[4]
            d = x[5]
            time_start = x[6]
            queue = x[7]
            slots = x[8]
            job_name = x[11]
            h_vmem=x[17].split('=')[1]
            row = [job_ID,prior,user,status,d,time_start, queue,
            slots,job_name,h_vmem]
            data.append(row)
      
            
        for x in queued:
            job_ID = x[0]
            prior = x[1]
            user = x[3]
            status = x[4]
            d = x[5]
            time_start = x[6]
            queue = np.NaN
            slots = x[7]
            job_name = x[10]
            h_vmem = x[13]
            row = [job_ID,prior,user,status,d,time_start, queue,
            slots,job_name,h_vmem]
            data.append(row)
            
        
    

        
        df = pd.DataFrame(data,columns=headers)
                
        return df
    


# In[ ]:

def SGE_Header(jobname, queue, h_vmem, cores, out, err):
    """
    USAGE: SGE_Header(jobname, queue, h_vmem, cores, out, err)'
    """
 
        
    if queue <> "all":
        Header = "#!/bin/bash" + "\n"         + "#$ -N " + jobname + "\n"         + "#$ -l h_vmem=" + h_vmem + "\n"         + "#$ -l " + queue + "\n"         + "#$ -pe smp " + cores + "\n"         + "#$ -V" + "\n"         + "#$ -e " + err + "\n"         + "#$ -o " + out
    else:
        Header = "#!/bin/bash" + "\n"         + "#$ -N " + jobname + "\n"         + "#$ -l h_vmem=" + h_vmem + "\n"         + "#$ -pe smp " + cores + "\n"         + "#$ -V" + "\n"         + "#$ -e " + err + "\n"         + "#$ -o " + out
    return Header


# In[37]:

def add_OL_Column_RO_Calls(df1, df2, overlap=0.5, vcf_format=True):
    """ Takes data frames with chr, start, end as input, converts to pybedtools
    performs reciprocal overlap analysis
    
    returns df of true/false and overlap percent with index of df1 (df)
    returns dataframe of bedtools overlap output (Z)
    
    if vcf_format = True, converts df1 from VCF format to properly formatted
    pandas df
    
    overlap: sets reciprocal overlap thresholds (default = 0.5)

    """
    
    Header_original = df1.columns.tolist()
    
    Remove_col = ['CHROM', 'POS', 'END','ID']
    
    if vcf_format:
        
## Fix the dataframe order if in VCF file format to make sure it can be 
## turned into a pybedtool

        for x in Remove_col:
            try:
                Header_original.remove(x)
            except:
                print x
        HEAD = Remove_col+ Header_original
    if not vcf_format:
        HEAD = Header_original

    df1=df1[HEAD]
    
## Modify the header names to attribute to A or B bedtools 


    Heading_1 = Header_list(df1.columns.tolist(),'A')
    Heading_2 = Header_list(df2.columns.tolist(), 'B')
    
## New Header for bedtool intersect -wo output dataframes

    Out_Head = Heading_1 + Heading_2 + ['Overlap']

    
    BT_B = pybedtools.BedTool.from_dataframe(df2)
    BT_A = pybedtools.BedTool.from_dataframe(df1)

    
    Z = pd.read_table(BT_A.intersect(BT_B, f=overlap, F=overlap, wo=True).fn, names=Out_Head)
    Z = Z.convert_objects(convert_numeric=True)
    df = pd.DataFrame(0, index = df1.Coords, columns=['Overlap', 'Percent_A', 'Percent_B'])
    Z['Percent_A'] = Z.Overlap/(Z.End_A - Z.Start_A)
    Z['Percent_B'] = Z.Overlap/(Z.End_B - Z.Start_B)
    df['Overlap'] = 'False'

    ind = 0        
    for x in Z.Coords_A:
        try:
            
            
            df.ix[x,'Overlap']= 'True'
            df.ix[x, 'Percent_A']= Z.ix[ind,'Percent_A']
            df.ix[x, 'Percent_B']= Z.ix[ind,'Percent_B']
            ind += 1 
        except:
            print x
    return df, Z, Out_Head


# In[1]:

def svviz_command(family, svtype, region, out_file):
    Fam = family
    Fam_Spl = Fam.split('_')
    bams = []
    
    Coords = ''
    
    reg_spl = region.split('_')
    for x in reg_spl:
        Chr = reg_spl[1]
        Start = reg_spl[2]
        End = reg_spl[3]
        Coords = " ".join([Chr,Start,End])
    
    
    
    for x in Fam_Spl:
        UUID = Subj_dict[x]
        bam_loc = '/frazer01/projects/CARDIPS/pipeline/WGS/alignment_bwa/' + UUID + '/' + UUID + '.mdup.bam'
        
        bams.append(bam_loc)
        
        
    child = bams[0]
    father = bams[1]
    mother = bams[2]
    
    reference = '/frazer01/publicdata/gatk_bundle_2.8/b37/human_g1k_v37_decoy_Sendai.fasta'
    command = 'svviz -t ' + svtype + ' -b ' + child +     ' -b ' + father + ' -b ' + mother + " " + reference + " " + Coords +     ' --export ' + out_file
     
    
    print command

