#!/usr/bin/python3


import sys

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]

d = {}

with open(input1) as vcf1, open(input2) as vcf2, open(input3) as vcf3, open("Merged_Trio_SVs.bed",'w') as merged:
	for line in vcf1.readlines():
		if line.startswith('#'):
			continue
		else:
			line = line.rstrip()
			cols = line.split('\t')
			end = cols[7].split(';')[0].split('=')[1]
			cols.pop(2)
			cols.insert(2,str(end))
			keeps = cols[0:7]
			keeps.append(cols[8])
			keytrio = '_'.join(keeps)
			d[keytrio] = [cols[9],'NA','NA']
	for line in vcf2.readlines():
		if line.startswith('#'):
			continue
		else:
			line = line.rstrip()
			cols = line.split('\t')
			end = cols[7].split(';')[0].split('=')[1]
			cols.pop(2)
			cols.insert(2,str(end))
			keeps = cols[0:7]
			keeps.append(cols[8])
			keytrio = '_'.join(keeps)
			if keytrio in d.keys():
				d[keytrio].pop(2)
				d[keytrio].insert(1,cols[9])
			else:
				d[keytrio] = ['NA',cols[9],'NA']
	for line in vcf3.readlines():
		if line.startswith('#'):
			continue
		else:
			line = line.rstrip()
			cols = line.split('\t')
			end = cols[7].split(';')[0].split('=')[1]
			cols.pop(2)
			cols.insert(2,str(end))
			keeps = cols[0:7]
			keeps.append(cols[8])
			keytrio = '_'.join(keeps)
			if keytrio in d.keys():
				d[keytrio].pop(2)
				d[keytrio].insert(2,cols[9])
			else:
				d[keytrio] = ['NA','NA',cols[9]]
	merged.write('CHROM\tSTART\tEND\tREF\tALT\tQUAL\tFILTER\tFORMAT\tNA19240\tHG00514\tHG00733\n')
	for k in sorted(d.keys()):
		SVlist = k.split('_')
		SVlist.append(d[k][0])
		SVlist.append(d[k][1])
		SVlist.append(d[k][2])
		SVs = '\t'.join(SVlist)
		merged.write(SVs+'\n')



