
import os
from Functions import *


###################################################
### PART 1 | Rank based on correlation strength ###
### PART 2 | Check for same/different mCs       ###
###################################################

signCorrs={}
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	signCorrs['%s--%s' % (i[0],i[1])]=i[0]

### For checking if genes correlate with same/different mCs:
G       = {}
tadG    = {}
tadC    = {}

### For ranking genes based on correlation stength
wholeD  = {}
regionD = {'Overlap':{}, 'Intermediate':{}, 'Long':{}}

n=0
checkN=[2**x for x in range(30)]
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		
		n+=1
		if n in checkN:
			print n
		
		i=i.replace('\n','').split('\t')
		
		### For checking same/different mCs:
		try:
			G[i[1]]+=[i[0]]
		except KeyError:
			G[i[1]]=[i[0]]
		i[4]=i[4].split('|')[0]
		try:
			tadG[i[4]]+=[i[1]]
		except KeyError:
			tadG[i[4]]=[i[1]]
		try:
			tadC[i[4]]+=[i[0]]
		except KeyError:
			tadC[i[4]]=[i[0]]
		i=f.readline()
		
		### For ranking:
		if i[2]=='nan':
			i=f.readline()
			continue
		
		cGene = i[1].split('|')[1]
		cDist = int(i[5])
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

### For ranking:
s=''
for i in wholeD.keys():
	s+='%s\t%s\n' % (i,max(wholeD[i]))
write_file('output/Correlation/Step3/RankTads.Whole.txt',s)

for r in regionD.keys():
	s=''
	for i in regionD[r].keys():
		s+='%s\t%s\n' % (i,max(regionD[r][i]))
	write_file('output/Correlation/Step3/RankTads.%s.txt' % (r),s)


### For checking same/different mCs:
tKeys=tadG.keys()

with open('output/Correlation/Step3/TADsets.GenesCsInTAD.txt','w') as f:
	f.write('TAD\tTadCs\tGene1\tGene2\tGene1Cs\tGene2Cs\tComCs\n')
	for t in tKeys:
		print t,tKeys.index(t),len(tKeys)
	
		tc=set(tadC[t])
		allGenes=list(set(tadG[t]))
		for i in range(len(allGenes)-1):
			gi=allGenes[i]
			ci=set(G[gi]).intersection(tc)
			sci=[get_from_dict(signCorrs,'%s--%s' % (x,gi),'NA') for x in ci]
			sci=set([x for x in sci if x!='NA'])
			for j in range(i+1,len(allGenes)):
				gj=allGenes[j]
				cj=set(G[gj]).intersection(tc)
				scj=[get_from_dict(signCorrs,'%s--%s' % (x,gj),'NA') for x in cj]
				scj=set([x for x in scj if x!='NA'])
				s=[t,len(list(tc)),gi,gj,len(list(sci)),len(list(scj)),len(list(sci.intersection(scj)))]
				s='\t'.join([str(x) for x in s])+'\n'
				f.write(s)





