#!/usr/bin/Python3

import sys
input = sys.argv[1]

with open(input) as text, open('battle_2015.bed','w') as bed:
	for line in text.readlines():
			line = line.rstrip()
			cols = line.split('\t')
			stuff = cols[3]+'\t'+cols[0]+'\t'+cols[2]+'\t'+'eQTL'+'\t'+''+'\t'+'CHS'+'\t'+'Battle_2015'
			#name = 'gene_name=='+cols[1]+'rsid=='+cols[0]
			#name = 'rsid=='+cols[1]
			#name = 'gene_name=='+cols[2]
			#chr = cols[0]
			#pos = cols[1]
			#start = str(int(cols[1])-1)
			#bedform = [chr, start, pos, name]
			#mybed = '\t'.join(bedform)
			bed.write(stuff+'\n')