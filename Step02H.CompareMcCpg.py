
import os
from Functions import *

i2cpg={}
cpg2i={}
with open('bed/BedIntersect.Cs.CpGi.bed') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		
		cpg2i['%s.%s.NA0' % (i[0],i[2])]=i[6].replace(':','')
		
		try:
			i2cpg[i[6].replace(':','')]+=['%s.%s.NA0' % (i[0],i[2])]
		except KeyError:
			i2cpg[i[6].replace(':','')]=['%s.%s.NA0' % (i[0],i[2])]
		
		i=f.readline()

### mCs with correlations
mcCorr={}
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	mcCorr[i.split('\t')[0]]='Corr'

### CpGi with correlations
cpgiCorr={}
for i in read_file('output/Correlation/Step2/CpGiCorrelations.txt')[1:]:
	cpgiCorr[i.split('\t')[0]]='Corr'

### mCs in background
mcBack=[]
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.split('\t')
		if i[2]!='nan':
			mcBack+=[i[0]]
		i=f.readline()
mcBack=list(set(mcBack))

### CpGi in background
cpgiBack=[]
with open('output/Correlation/Step1/CpGiCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.split('\t')
		if i[2]!='nan':
			cpgiBack+=[i[0]]
		i=f.readline()
cpgiBack=list(set(cpgiBack))

print 3

s=''
for i in mcBack:
	
	# is Correlated?
	mc_is_corr=get_from_dict(mcCorr,i,'NS')
	
	# find CpGi
	cIsland=get_from_dict(cpg2i,i,'NA')
	
	# is CpGi correlated?
	cpgi_is_corr=get_from_dict(cpgiCorr,cIsland,'NS')
	
	s+='\t'.join([i,mc_is_corr,cIsland,cpgi_is_corr])+'\n'

write_file('output/mC_CpGi_Correspond.txt',s)


s=''
for i in cpgiBack:
	
	# is Correlated?
	cpgi_is_corr=get_from_dict(cpgiCorr,i,'NS')

	# find mCs
	cMC=i2cpg[i]

	# is mC correlated?
	cMC=[get_from_dict(mcCorr,x,'NS') for x in cMC]
	if True in [x!='NS' for x in cMC]:
		cMC='Corr'
	else:
		cMC='NotCorr'
	
	s+='\t'.join([i,cpgi_is_corr,cMC])+'\n'

write_file('output/CpGi_mC_Correspond.txt',s)



