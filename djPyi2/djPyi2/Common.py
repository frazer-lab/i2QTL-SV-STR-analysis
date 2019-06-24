
# coding: utf-8

# In[13]:

import os
import datetime
import sys
import glob
import pandas as pd
import numpy as np
import datetime

def flatten_list(l):
    """flatten a list of lists"""
    x = [i for z in l for i in z]
    return x


def split_seq(seq, num_pieces):
    """ split a list into a specified number of equally sized pieces (last piece may be a different size) 
    returns a generator of each chunk of the original list"""
    start = 0
    for i in xrange(num_pieces):
        stop = start + len(seq[i::num_pieces])
        yield seq[start:stop]
        start = stop

def datestring(month=True, day=True, year=True, hour=False, minute=False):
    d = datetime
    
    y = str(d.datetime.now().timetuple().tm_year)
    m =  str(d.datetime.now().timetuple().tm_mon)
    D = str(d.datetime.now().timetuple().tm_mday)
    h = str(d.datetime.now().timetuple().tm_hour)
    mins = str(d.datetime.now().timetuple().tm_min)
    h_original = h
    
    if int(month) in range(0,10):
        m = "0" + m
    
    if hour == True:
        if int(h) > 12:
            h = str(int(h) - 12) + 'PM'
        
        elif int(h)<= 12:
            h = h + 'AM'
    
    dt_array = [y, m, D, h, mins]
    t_f_array = [year, month, day, hour, minute]
    
    Date_Array = []
    for s, i, in zip(dt_array, t_f_array):
        if i == True:
            Date_Array.append(s)
    
    return '_'.join(Date_Array)
    
normal_Chrs = range(1,23)
normal_Chrs.extend(['X','Y'])
normal_Chrs = [str(x) for x in normal_Chrs]
autosomes =[str(x) for x in range(1,23)]

def save_dataframe(var_name, df, folder_loc, alt_name = False, pickle= True,  csv = True, 
                      sep = "\t", raw_sep = r"\t", csv_suffix = '.tsv', rewrite=False, return_record_fn = False,
                      print_records =True, print_only_csv = False, print_only_pickle=False, print_fn_eq = False,
                      print_vars_recorded_loc=True,reset_index=False,**kwargs):
    """  
    purpose:
    save a pandas dataframe as pickle and tsv (or any combination thereof), while logging the save to a file
    with the needed code to load the variable elsewhere
    print the location of saves and commands to load them (by default)"""
    

    date_str = datestring(hour = True)
    record_file = os.path.join(folder_loc, 'load_saved_nb_variables.py')
    pickle_record_file = os.path.join(folder_loc, 'load_pickled_nb_variables.py')
    
    if return_record_fn == True:
        
        for f in [record_file, pickle_record_file]:
            if os.path.exists(f):
                print "# Record File Exists: {}".format(f)
    
    
    rw_kw = 'a'
    if not os.path.exists(record_file):
        rw_kw = 'w'
        
    elif rewrite == True:
        rw_kw = 'w'
    else:
        rw_kw = 'a'
        
    if alt_name==False:
        file_name = var_name
    else:
        file_name = alt_name
    
    
    with open(record_file, rw_kw) as F:
        with open(pickle_record_file, rw_kw) as P:
            
            out_fns = []
            out_pickles = []
            out_pickles.append("# " + date_str)
            out_fns.append("# " + date_str)
        
            fn = os.path.join(folder_loc, file_name)
            pickle_fn = fn + '.pkl'
            csv_fn = fn + csv_suffix
            csv_load = "{} = pd.read_csv('{}', sep='{}')".format(var_name, csv_fn, raw_sep)
            pickle_load = "{} = pd.read_pickle('{}')".format(var_name, pickle_fn)
            
            csv_fn_eq = "fn = '{}'".format(csv_fn)

            if csv == True:
                
                if reset_index ==True:
                    df.reset_index().to_csv(csv_fn, sep=sep, **kwargs)
                    out_fns.append(csv_load)
                else:

                    df.to_csv(csv_fn, sep=sep, **kwargs)
                    out_fns.append(csv_load)

                
            if pickle==True:
                df.to_pickle(pickle_fn)
                
                out_fns.append(pickle_load)
                out_pickles.append(pickle_load)
                
                P.write("\n".join(out_pickles))
                
            
            if print_records == True:
                if not print_only_csv:
                    if pickle == True:
                        print pickle_load
                if not print_only_pickle:
                    
                    if csv == True:
                        print csv_load
                    if print_fn_eq:
                        print csv_fn_eq

                else:
                    print ''
            
            out_string ="\n".join(out_fns)
            
            F.write(out_string)
    
    if print_vars_recorded_loc ==True:
        
    
        print '# all vars recorded: ' + record_file
        print '# pickled vars recorded:' + pickle_record_file

import gzip

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

