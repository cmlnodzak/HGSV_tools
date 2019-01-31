#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:51:52 2019

@author: mitchnodzak
"""
from math import floor, log10, fabs
from pybedtools import BedTool
import re


# return the order of magnitude of a number.
def orderOfMag(val):
    n = floor(log10(fabs(val)))
    return n

# URL = ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
    
#gtf ='gencode.v25.annotation.gtf'

# extract pandas dataframe from a GTF for a particular feature of interest.
def gencodeReader(gtf, feat):
    anno_df = BedTool(gtf).to_dataframe()
    anno_df = anno_df[anno_df.feature == feat]
    anno_df['gene_id'] = anno_df['attributes'].map(lambda x: re.search(r'gene_id "ENSG\d+.\d+', x)[0].split()[1].strip("\""))    
    anno_df['gene_type'] = anno_df['attributes'].map(lambda x: re.search(r'gene_type\s"\w+',x)[0].split()[1].strip("\""))
    anno_df['gene_name'] = anno_df['attributes'].map(lambda x: re.search(r'gene_name\s"\w+',x)[0].split()[1].strip("\""))
    anno_df['gene_status'] = anno_df['attributes'].map(lambda x: re.search(r'gene_status\s"\w+',x)[0].split()[1].strip("\""))
    anno_df = anno_df.drop('attributes',axis=1)
    return anno_df

attr = {'gene_id': None, 'gene_type': None , 
        'gene_name': None, 'gene_status' : None}