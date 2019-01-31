#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 10:05:54 2019

@author: mitchnodzak

Revisions to ASE analysis based on python rather than command line tools.

"""
import pandas as pd
from pybedtools import BedTool

# grep only passing variants from integrated VCF.
# saved as file name below...
# 'PASS_Illumina_Integrate_20170206.ALL.vcf' as infile.

infile = 'PASS_Illumina_Integrate_20170206.ALL.vcf'
#infile = sys.argv[1]

svtype = infile.split('/')[-1].split('.')[1]
pref = '.'.join(infile.split('/')[-1].split('.')[0].split("_")[0:2])
outfy = pref+'.'+svtype+'_integrate.bed'


col_names =  ['#chrom', 'start', 'stop', 'REF', 'ALT', 'QUAL', 'FILT', 'sample', 'GT']
hgsv_list = ['HG00512','HG00513','HG00514', 'HG00731', 'HG00732', 'HG00733', 'NA19238', 'NA19239', 'NA19240']

# Extract the relevant information to bed file for safe keeping, and a pandas 
# DataFrame and generate a subsequent coordinate intersections 
# with genome annotations.
class VCFparser:
    def __init__(self,infile):
        self.bedtoolsList = []
        self.frame  = pd.DataFrame(columns = col_names)
        
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
                                    self.frame = self.frame.append(call)
                                    outfile.write(call+'\n')
        self.bedtoolsList = self._Frame2Bed(self.frame,hgsv_list)
        self._getBedTools(hgsv_list)
        
#pass in a list of sample names as they are written in the file
    def _Frame2Bed(df,samples):
        bedList = []
        # extract heterozygous deletions and insertions, group by samples.
        df = df.loc[df['GT'] == '0/1']
        df = df.groupby('sample')
        for i in samples:
            i = df.get_group(i).drop_duplicates()
            i_bt = BedTool.from_dataframe(i)
            bedList.append(i_bt)
        return bedList

    def _getBedTools(self,samples):
        for k in samples:
            k = self.bedtoolsList[k]
            print('creating bedtool for sample ' + k +' \n')
            return k


# x = BedTool.from_dataframe(df)

infile = 'PASS_Illumina_Integrate_20170206.ALL.vcf'

VCFparser(infile)



                            
                            
                            
                                