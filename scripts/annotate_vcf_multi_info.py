
# coding: utf-8

# In[27]:

from __future__ import division
import numpy as np
import os
import sys
import datetime
from subprocess import call

import itertools
import tempfile
import pandas as pd

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


# In[31]:

def annotate_vcf_info(fn, annotation_file, header_info_fn, out_file_name, add_info=False, add_filter=False, desired_samples= False, hf_all=False):
    
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
    
    header_info = pd.read_table(header_info_fn)
#     print header_info
    
    annotate_for_all = False
    if hf_all:
        hf_all_df = pd.read_table(hf_all)
        
        annotate_for_all = hf_all_df.ANNOTATION.tolist()
        header_info = pd.concat([header_info, hf_all_df], ignore_index=True)
        
      
    
    count = 0
    header_end = find_header_end(fn)
    add_info = False
    add_filters = False

    for line in F:
        line = line.rstrip()
        lin_spl = line.split('\t')

        if count < (header_end - 1):
            if desired_samples:

                if line.find('##SAMPLE') == 0:
                    samp = line.split("=").pop()[:-1]
                    if samp in desired_samples:
                        annotated_vcf.write(line + '\n')
                        count +=1
                        continue
                    else:
                        count +=1
                        continue
            
                
                
            if not add_info:

                if line.find('##INFO') == 0:
                    add_info = True
                    
                    for i, x in header_info.iterrows():
                        annotation_id = x['ID_COL']
                        number = x['NUMBER']
                        info_type = x['TYPE'] 
                        description = x['DESCRIPTION']
                        

                        l1 = '##INFO=<ID={},Number={},Type={},Description="{}">'.format(annotation_id,
                                                                                        number, info_type,
                                                                                        description)

                        annotated_vcf.write(l1 + '\n')
                    
            

            annotated_vcf.write(line + '\n')



        if count == header_end-1:

       
            header = lin_spl
            header_dict = {l:i for l,i in zip(lin_spl, range(0, len(lin_spl)))}
            samples = header[9:]
            info_cols = header[:9]
            sample_order_dict = {s:header.index(s) for s in samples}
            if desired_samples:
                try:
                    reodered_samples = [header[sample_order_dict[i]] for i in desired_samples]
                except:
                    print sample_order_dict
                    print desired_samples
                    print samples
                    break
                header_mod = info_cols + reodered_samples
                annotated_vcf.write("\t".join(header_mod) + '\n')
            else:
                annotated_vcf.write(line + '\n')

        elif count > (header_end - 1):

            # first format line, but this might be different for different lines so I'll calculate it on every line
            info = lin_spl[header_dict['INFO']]
            ID = lin_spl[header_dict['ID']]

            # if annotation is in the mappings add it
            

            try:
                annot = mappings[ID]
                info_out =  info + ";" + annot                    
                lin_spl[header_dict['INFO']] = info_out
                

            except:
                pass
            
            
            if annotate_for_all != False:
                info = lin_spl[header_dict['INFO']]
                info_out =  info + ";" + ";".join(annotate_for_all)
                lin_spl[header_dict['INFO']] = info_out
            
            if desired_samples:
                info_cols = lin_spl[:9]
                reodered_samples = [lin_spl[sample_order_dict[i]] for i in desired_samples]
                out_line_mod = info_cols + reodered_samples
         
                annotated_vcf.write("\t".join(out_line_mod) + '\n')
            
            else:
                
                annotated_vcf.write("\t".join(lin_spl) + '\n')


        count +=1     
        
    annotated_vcf.close()
    F.close()


# In[3]:

def add_arguments_to_parser(parser):
    
    parser.add_argument("-vcf", "--vcf", dest="vcf", metavar='<vcf_fn>', help="vcf filename", required=True)
    parser.add_argument("-af", "--annot_file", dest="annotation_fn", metavar='<annotation fn>', help="tsv file that in the form: ID  Annotation - for annotating the vcf- IDs are expected to be unique- no header", required=True)
    parser.add_argument("-hf", "--header_file", dest="header_annot_fn", metavar='<header info fn>', help="tsv file that in the form of: ID_COL  NUMBER  TYPE  DESCRIPTION - with header- to annotate multiple info annotations in one step - must have this header exactly- TYPE can be- String/Float/Integer/Flag NUMBER can be 1 or . if unknown (number of items)", required=True)
    parser.add_argument("-s", "--samples", dest="samples_fn", metavar='<samples fn>', help="file of samples desired in VCF output in the order desired, if not provided will output all samples in their current order ", required=False, default=False)
    
#     parser.add_argument("-tag_all", "--tag_all", dest="tag_all", action = 'store_true', help="tag_all_samples with flags", required=False)
    
    parser.add_argument("-hf_all", "--header_info_all", dest='header_info_all_fn',  help="tsv file that in the form of: ID_COL  NUMBER  TYPE  DESCRIPTION ANNOTATION- with header- to annotate multiple info annotations- ANNOTATION contains the annotation as it will be added into the info col for all sites", default=False)
                        
    
    parser.add_argument("-o", "--output_fn", dest="output_fn", metavar='<out_fn>', help="output vcf filename", required=True)
    
    parser.set_defaults(entry_point=run_from_args)
    


# In[ ]:

def command_parser():
    parser = argparse.ArgumentParser(description= 'command line utility to add an annotation to a vcf file given unique identifiers')
    add_arguments_to_parser(parser)
    return parser


# In[28]:

def run_from_args(args):
    vcf_fn = args.vcf
#     id_name = args.id
#     num = args.number
#     description = args.description
    out_fn = args.output_fn
#     info_type = args.type
    annotation_fn = args.annotation_fn
    header_info_fn = args.header_annot_fn
    samples_fn = args.samples_fn
    samples = [l.rstrip() for l in open(samples_fn)]
#     print samples
    
    print 'starting annotation:', str(datetime.datetime.utcnow())
    annotate_vcf_info(vcf_fn, annotation_fn, header_info_fn, out_fn, desired_samples=samples, 
                      hf_all = args.header_info_all_fn)
    
    
    print "vcf annotated:", str(datetime.datetime.utcnow()) 
    print out_fn
    


# In[ ]:

if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    sys.exit(args.entry_point(args))

