#!/usr/bin/python3
# Cleans up weird format of merged Pacbio vcf file for insertions and deletions.

import sys

input = sys.argv[1]
pref = input.split('/')[-1].split('.')[0]
outfy = pref+'_all.vcf'

with open(input) as hetero, open(outfy,'w') as fixed:
	for line in hetero.readlines():
		if line.startswith('#'):
			continue
		else:
			line = line.rstrip()
			cols = line.split('\t')
			end = col[7].split(';')[0].split('=')[1]
			cols.pop(2)
			cols.insert(2,str(end))
			keeps = cols[0:7] + cols[8:]
			keep = '\t'.join(keeps)
			fixed.write(keep+'\n')
			else:
				continue
			
