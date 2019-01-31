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


col_names =  ['#chrom', 'start', 'stop', 'REF', 'ALT', 'QUAL', 'FILT', 'sample', 'GT']

df  = pd.DataFrame(columns = col_names)
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

# extract heterozygous deletions and insertions, group by samples.
df = df.loc[df['GT'] == '0/1']

df = df.groupby('sample')
samp512 = df.get_group('HG00512').drop_duplicates()
samp513 = df.get_group('HG00513').drop_duplicates()
samp514 = df.get_group('HG00514').drop_duplicates()
samp731 = df.get_group('HG00731').drop_duplicates()
samp732 = df.get_group('HG00732').drop_duplicates()
samp733 = df.get_group('HG00733').drop_duplicates()
samp238 = df.get_group('NA19238').drop_duplicates()
samp239 = df.get_group('NA19239').drop_duplicates()
samp240 = df.get_group('NA19240').drop_duplicates()

HG00512 = BedTool.from_dataframe(samp512)
HG00513 = BedTool.from_dataframe(samp513)
HG00514 = BedTool.from_dataframe(samp514)
HG00731 = BedTool.from_dataframe(samp731)
HG00732 = BedTool.from_dataframe(samp732)
HG00733 = BedTool.from_dataframe(samp733)
NA19238 = BedTool.from_dataframe(samp238)
NA19239 = BedTool.from_dataframe(samp239)
NA19240 = BedTool.from_dataframe(samp240)



# x = BedTool.from_dataframe(df)







                            
                            
                            
                                