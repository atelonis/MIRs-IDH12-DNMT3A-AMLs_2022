
import os
import numpy

from Functions import *

wholeD  = {}
regionD = {'Overlap':{}, 'Intermediate':{}, 'Long':{}}

n=0
checkN=[2**x for x in range(30)]
with open('TCGA/output/Correlation/Step1/Correlations.txt') as f:
	i=f.readline()
	while i:
		
		n+=1
		if n in checkN:
			print n
		
		i=i.replace('\n','').split('\t')
		if i[2]=='NA':
			i=f.readline()
			continue
		
		cGene = i[1].split('|')[0]
		if i[5]=='GeneBody':
			cDist=0
		else:
			cDist = int(i[5].split('_')[1])
		cCorr = abs(float(i[2]))
		
		try:
			wholeD[cGene]+=[cCorr]
		except KeyError:
			wholeD[cGene]=[cCorr]
		
		if cDist<2000:
			cRegion='Overlap'
		elif cDist<500000:
			cRegion='Intermediate'
		else:
			cRegion='Long'
		
		try:
			regionD[cRegion][cGene]+=[cCorr]
		except KeyError:
			regionD[cRegion][cGene]=[cCorr]
		
		i=f.readline()

s=''
for i in wholeD.keys():
	s+='%s\t%s\n' % (i,max(wholeD[i]))
write_file('TCGA/output/Correlation/Step3/RankTads.Whole.txt',s)

for r in regionD.keys():
	s=''
	for i in regionD[r].keys():
		s+='%s\t%s\n' % (i,max(regionD[r][i]))
	write_file('TCGA/output/Correlation/Step3/RankTads.%s.txt' % (r),s)


