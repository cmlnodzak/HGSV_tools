#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:05:54 2019

@author: mitchnodzak

Revisions to ASE analysis based on python rather than command line tools.

"""
import sys
import numpy as np
import pandas as pd
from pybedtools import BedTool
from collections import OrderedDict

# grep only passing variants from integrated VCF.
# saved as file name below...
# 'PASS_Illumina_Integrate_20170206.ALL.vcf' as infile.
 
infile = sys.argv[1]
svtype = infile.split('/')[-1].split('.')[1]
pref = '.'.join(infile.split('/')[-1].split('.')[0].split("_")[0:2])
outfy = pref+'.'+svtype+'_integrate.bed'


df = pd.DataFrame()

# Extract the relevant information to bed file for safe keeping, and a pandas 
# DataFrame and generate a subsequent coordinate intersections 
# with genome annotations.
class VCFparser:
    def __init__(self,infile):
        with open(infile) as SVs, open(outfy,'w') as outfile:
                for line in SVs.readlines():
                        if line.startswith('#'):
                                continue
                        else:
                                line = line.rstrip().split('\t')
                                end = line[7].split(';')[1].split('=')[1]
                                line.pop(2)
                                line.insert(2,end)
                                data = []
                                info = line[7].split(';')[7].split(',')
                                for genos in info:
                                        geno = genos.split(':')
                                        sampid = str(geno[4])
                                        genetype = str(geno[3])
                                        triotype = sampid + '\t' + genetype
                                        SVdata=line[0:7]
                                        SVdata.append(triotype)
                                        keep="\t".join(SVdata)
                                        data.append(keep)
                                for call in data:
                                    df = df.append(call)
                                    outfile.write(call+'\n')




x = BedTool.from_dataframe(df)





                            
                            
                            
                                