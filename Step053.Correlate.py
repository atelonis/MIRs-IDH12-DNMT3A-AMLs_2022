
import os
import sys
import time
import numpy
import scipy
import random
import itertools
from scipy import stats
from multiprocessing import Pool

from Functions import *

def corr_1(l): # [meth, ENSG|gene, TAD, distance, gene]
	
	try:
		cMeth=read_file('TEMP_METH/%s.txt' % (l[0]))
		cMrna=read_file('TEMP_MRNA/%s.txt' % (l[4]))
	except IOError:
		return '\t'.join([ l[0], l[1], 'NA', 'NA', l[2], l[3] ])+'\n'
	
	if 'NA' in cMeth or 'NA' in cMrna:
		return '\t'.join([ l[0], l[1], 'NA', 'NA', l[2], l[3] ])+'\n'
	
	l0=[float(x) for x in cMeth]
	l1=[float(x) for x in cMrna]
	
	cSpear=scipy.stats.spearmanr( l0,l1 )
	r='\t'.join([ l[0], l[1], str(cSpear[0]), str(cSpear[1]), l[2], l[3] ])+'\n'
	return r

create_dir('TEMP_MRNA')
create_dir('TEMP_METH')
create_dir('TCGA/output/Correlation')
create_dir('TCGA/output/Correlation/Step1')

mRNA=read_file('TCGA/RnaData.txt')
meth=read_file('TCGA/MethylationData.txt')

mrnaSamples=mRNA[0].split('\t')
methSamples=meth[0].split('\t')

wOrder=[mrnaSamples.index(x) for x in methSamples]

for i in mRNA[1:]:
	i=i.split('\t')
	cValues=[i[x+1] for x in wOrder]
	write_file('TEMP_MRNA/%s.txt' % (i[0].split('|')[0]), '\n'.join(cValues)+'\n')

for i in meth[1:]:
	if 'NA' in i:
		continue
	i=i.split('\t')
	write_file('TEMP_METH/%s.txt' % (i[0]),'\n'.join(i[1:])+'\n')

allCombs=[]
with open('output/TAD.pairs.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		i[0]=i[0].split('.')[0]
		i+=[i[1].split('|')[1]]
		allCombs+=[ i ]
		i=f.readline()
random.shuffle(allCombs)

p=Pool(117)
s=p.map(corr_1,allCombs)
write_file('TCGA/output/Correlation/Step1/Correlations.txt',''.join(s))


