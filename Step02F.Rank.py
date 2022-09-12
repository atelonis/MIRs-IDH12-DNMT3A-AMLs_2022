
import os
import numpy
from Functions import *

create_dir('output/Correlation/Step3/CpGi')

wholeD  = {}
regionD = {'Proximal':{}, 'Intermediate':{}, 'Long':{}}

n=0
checkN=[2**x for x in range(30)]
with open('output/Correlation/Step1/CpGiCorrelations.txt') as f:
	i=f.readline()
	while i:
		
		n+=1
		if n in checkN:
			print n
		
		i=i.replace('\n','').split('\t')
		if i[2]=='nan':
			i=f.readline()
			continue
		
		cGene   = i[1].split('|')[1]
		cCorr   = abs(float(i[2]))
		cRegion = i[4]
		
		try:
			wholeD[cGene]+=[cCorr]
		except KeyError:
			wholeD[cGene]=[cCorr]
		
		try:
			regionD[cRegion][cGene]+=[cCorr]
		except KeyError:
			regionD[cRegion][cGene]=[cCorr]
		
		i=f.readline()

s=''
for i in wholeD.keys():
	s+='%s\t%s\n' % (i,max(wholeD[i]))
write_file('output/Correlation/Step3/CpGi/RankTads.Whole.txt',s)

for r in regionD.keys():
	s=''
	for i in regionD[r].keys():
		s+='%s\t%s\n' % (i,max(regionD[r][i]))
	write_file('output/Correlation/Step3/CpGi/RankTads.%s.txt' % (r),s)


