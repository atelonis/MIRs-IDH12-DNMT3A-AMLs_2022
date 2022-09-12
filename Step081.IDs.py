
import os
from Functions import *

ensg2own={}
for i in read_file('HumanGenome/GeneRanges.bed')[1:]:
	i=i.split('\t')[3]
	ensg2own[i.split('|')[0]]=i

affy2own={}
for i in read_file('Hematopoiesis/mart_export.txt'):
	# file from ensembl/biomart with fields:
	#	Gene stable ID version
	#	Gene name
	#	Gene type
	#	AFFY HG U133A probe
	i=i.split('\t')
	if i[2]!='protein_coding' or i[3]=='':
		continue
	try:
		affy2own[i[3]]=ensg2own[i[0].split('.')[0]]
	except KeyError:
		continue

Az09=map(chr,range(97,123))+map(chr,range(65,91))+[str(x) for x in range(10)]

s=''
for i in read_file('Hematopoiesis/GPL4685-15513.txt'):
	if i[0]=='#':
		continue
	i=i.split('\t')
	if i[0]=='ID':
		continue
	if i[8]=='' or i[9]=='' or i[8]=='NA':
		continue
	if False in [x in Az09 for x in list(i[8])+list(i[9])]:
		continue
	s+='%s\t%s|%s\n' % (i[0],i[8],i[9])

write_file('Hematopoiesis/LookupTable.txt',s)



