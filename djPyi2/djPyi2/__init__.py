import os
import sys
import datetime

import gzip 
from os.path import join
import subprocess
  
import pandas as pd
  

def _get_project_dir():
    return os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-3])
  
root = _get_project_dir()
  
## functions


  
def makedir(p):
    
    """Make a directory if it doesn't already exist"""
    try:
        os.makedirs(p)
    except OSError:
        pass


def convert_ipynb_to_script(fn):
    """ convery some ipynb to a script call this from within another notebook, where you are
    likely calling the script to make a .py after doing edits to .ipynb"""
    fn_script = fn.replace('ipynb', 'py')
    command = "jupyter nbconvert --to script {} --output {}".format(fn,fn_script[:-3])
    subprocess.call(command, shell=True)
    #! jupyter nbconvert --to script {fn} --output {fn_script[:-3]}
    return fn_script
    




def ucsc_format_cnv_id(cnv_id):
    spl = cnv_id.split('_')
    
    out = "chr{}:{}-{}".format(spl[1], spl[2], spl[3])
    return out
