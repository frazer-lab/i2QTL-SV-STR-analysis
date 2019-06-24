
# coding: utf-8

# In[2]:

from __future__ import division
import numpy as np
import os
import sys
import datetime
from subprocess import call

import itertools
import tempfile


import argparse
import warnings
# silence pandas warnings...
warnings.filterwarnings("ignore")


# In[39]:


def find_header_end(fn):
    """ find header end of vcf or gzipped vcf"""
    
    if fn.split('.').pop() == 'gz':
        import gzip
        F = gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
        

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
                F.close()
                return count_out
        except:
            break
    
    F.close()


# In[38]:

def annotate_vcf_info(fn, annotation_ID, description, info_type, annotation_file, out_file_name, number):
    
    mappings = {}
    
    if fn.split('.')[-1] == '.gz':
        F =  gzip.open(fn, 'rU')
    else:
        F = open(fn, 'rU')
    
    with open(annotation_file, 'rU') as A:
        for line in A:
            line = line.rstrip()
            lin_spl = line.rstrip().split()
            mappings[lin_spl[0]] = lin_spl[1]    
    
    annotated_vcf = open(out_file_name, 'w')



    count = 0
    header_end = find_header_end(fn)
    add_info = False
    add_filters = False

    for line in F:
        line = line.rstrip()
        lin_spl = line.split('\t')

        if count < header_end:
            if not add_info:

                if line.find('##INFO'.format('')) == 0:
                    add_info = True

                    l1 = '##INFO=<ID={},Number={},Type={}, Description="{}">'.format(annotation_ID, number, info_type, description)

                    annotated_vcf.write(l1 + '\n')


            annotated_vcf.write(line + '\n')

        if count == header_end:
            annotate_vcf.write(line + '\n')

        if count == header_end-1:

            annotated_vcf.write(line + '\n')
            header = lin_spl
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            samples = header[9:]
            info_cols = header[:9]


        elif count > header_end:

            # first format line, but this might be different for different lines so I'll calculate it on every line
            info = lin_spl[header_dict['INFO']]
            ID = lin_spl[header_dict['ID']]

            # if annotation is in the mappings add it
            try:
                annot = mappings[ID]

                add = ";{}={}".format(annotation_ID, annot)
                info += add
                lin_spl[header_dict['INFO']] = info

            except:
                pass


            annotated_vcf.write("\t".join(lin_spl) + '\n')


        count +=1     
        
    annotated_vcf.close()
    F.close()


# In[ ]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-vcf", "--vcf", dest="vcf", metavar='<vcf_fn>', help="vcf filename", required=True)
    parser.add_argument("-id", "--id_name", dest="id", metavar='<id>', help="id for info column- name you want the info column to be called- (eg: SVTYPE, MAF etc)", required=True)
    parser.add_argument("-af", "--annot", dest="annotation_fn", metavar='<annotation fn>', help="tsv file that in the form: ID  Annotation for annotating the vcf- IDs are expected to be unique", required=True)
    
    parser.add_argument("-number", "--number", dest="number", metavar='<number>', help="number of records in annotation record, for example, if only one number is expected, use 1, if unknown use '.', ", required=True, default = '.')
    
    parser.add_argument("-t", "--type", dest="type", metavar='<String/Float/Integer/Flag>', help="data type of the annotation", required=True, default = 'String')
    
    parser.add_argument("-d", "--description", dest="description", metavar='<description of the annotation>', help="Description of annotation for the header of the vcf", required=True, default = 'String')

    
    parser.add_argument("-o", "--output_fn", dest="output_fn", metavar='<out_fn>', help="output vcf filename", required=True)
    
    parser.set_defaults(entry_point=run_from_args)
    


# In[ ]:




# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to add an annotation to a vcf file given unique identifiers')
    add_arguments_to_parser(parser)
    return parser


# In[ ]:

def run_from_args(args):
    vcf_fn = args.vcf
    id_name = args.id
    num = args.number
    description = args.description
    out_fn = args.output_fn
    info_type = args.type
    annotation_fn = args.annotation_fn
    
    annotate_vcf_info(vcf_fn, id_name, description, info_type, annotation_fn, out_fn, num)
    


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

