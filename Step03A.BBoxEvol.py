
import os
import re
from scipy.stats import hypergeom

from Functions import *


#############
### B-Box ###
#############
# https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC7279448/

BBOX=re.compile('GTTC.A..C')
BBOX_RC=re.compile('G..T.GAAC')

s=''
with open('HumanGenome/MIRs.Ext.fa') as f:
	iName=f.readline()
	iSeq=f.readline()
	while iName:
		iName=iName.replace('>','').split('|')
		iName=iName[3:]
		iName='|'.join(iName[:6])
		
		if BBOX.search(iSeq) or BBOX_RC.search(iSeq):
			s+=iName+'\n'

		iName=f.readline()
		iSeq=f.readline()

write_file('HumanGenome/MIRs_with_BBox.txt',s)


##################
### Divergence ###
##################

s=''
with open('HumanGenome/hg19.fa.out') as f:
	i=f.readline()
	i=f.readline()
	i=f.readline()
	i=f.readline()
	while i:
		while '  ' in i:
			i=i.replace('  ',' ')
		if i[0]==' ':
			i=i[1:]
		i=i.split(' ')
		if i[8]=='+':
			cStrand='+'
		else:
			cStrand='-'
		cRep='|'.join([i[10],i[9],i[4],cStrand,i[5],i[6]])
		s+='\t'.join([cRep,i[1],i[2],i[3]])+'\n'
		i=f.readline()

write_file('Divergence.txt',s)



#######################################
### MIR subfamilies in correlations ###

geneSets=['%s.%s' % (x.split('.')[1],x.split('.')[2]) for x in os.listdir('output/Correlation/Step4/')]
geneSets=[x for x in geneSets if 'TADset' in x]

geneMIRs={}
for i in geneSets:
	cMIRs=read_file('output/Correlation/Step4/Reps.%s.ENSG.txt' % (i))
	cMIRs=[x.split('\t')[0] for x in cMIRs if 'SINE/MIR' in x]
	geneMIRs[i]=set(cMIRs)

mcMIRs={}
for sb in ['Significant','Background']:
	cMIRs=read_file('output/Correlation/Step4/Reps.Cs.%s.txt' % (sb))
	cMIRs=[x.split('\t')[0] for x in cMIRs if 'SINE/MIR' in x]
	mcMIRs['mCs.%s' % (sb)]=set(cMIRs)

# 1. Genes
backFamilies={}
for i in list(geneMIRs['TADsets.Background']):
	i=i.split('|')[1]
	try:
		backFamilies[i]+=1
	except KeyError:
		backFamilies[i]=1

Nh=len(geneMIRs['TADsets.Background'])

s=''
for g in geneSets:
	if g=='TADsets.Background':
		continue
	cFams=[x.split('|')[1] for x in list(geneMIRs[g])]
	for i in list(set(cFams)):
		kh=cFams.count(i)
		sh=len(cFams)
		Mh=backFamilies[i]

		cHG=hypergeom_test(kh,sh,Mh,Nh)

		s+='\t'.join([str(x) for x in [g,i,kh,sh,Mh,Nh]+cHG])+'\n'

# 2. mCs
backFamilies={}
for i in list(mcMIRs['mCs.Background']):
	i=i.split('|')[1]
	try:
		backFamilies[i]+=1
	except KeyError:
		backFamilies[i]=1

Nh=len(mcMIRs['mCs.Background'])

cFams=[x.split('|')[1] for x in list(mcMIRs['mCs.Significant'])]
for i in list(set(cFams)):
	kh=cFams.count(i)
	sh=len(cFams)
	Mh=backFamilies[i]

	cHG=hypergeom_test(kh,sh,Mh,Nh)

	s+='\t'.join([str(x) for x in ['mCs',i,kh,sh,Mh,Nh]+cHG])+'\n'

write_file('output/MIR_Subfamilies.txt',s)


##################
### Insulators ###

insulators=read_file('PMID26216945.MirInsulators/BI.Insul.Reps.bed')
insulators=set([x.split('\t')[6] for x in insulators if 'SINE/MIR' in x])

# 1. Genes
backIns=insulators.intersection(geneMIRs['TADsets.Background'])

Mh=len(backIns)
Nh=len(geneMIRs['TADsets.Background'])

s=''
for g in geneSets:
	if g=='TADsets.Background':
		continue
	kh=len(geneMIRs[g].intersection(insulators))
	sh=len(geneMIRs[g])
	
	s+='\t'.join([str(x) for x in [g,kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'


# 2. mCs
backIns=insulators.intersection(mcMIRs['mCs.Background'])

Mh=len(backIns)
Nh=len(mcMIRs['mCs.Background'])
kh=len(mcMIRs['mCs.Significant'].intersection(insulators))
sh=len(mcMIRs['mCs.Significant'])
s+='\t'.join([str(x) for x in ['mCs',kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'

write_file('output/MIR_EnrichInsulators.txt',s)


#####################
### Histone marks ###

histones=['H3K27ac','H3K27me3','H3K4me1','H3K4me3']

s=''
for hi in histones:
	cHI=read_file('PMID31085557.Adelman/BI.%s.Repeats.bed' % (hi))
	cHI=[x.split('\t')[7] for x in cHI]
	cHI=set([x for x in cHI if 'SINE/MIR' in x])
	
	# 1. Genes
	backMIR=geneMIRs['TADsets.Background'].intersection(cHI)
	Mh=len(backMIR)
	Nh=len(geneMIRs['TADsets.Background'])
	
	for g in geneSets:
		if g=='TADsets.Background':
			continue
		kh=len(geneMIRs[g].intersection(cHI))
		sh=len(geneMIRs[g])
		
		s+='\t'.join([str(x) for x in [hi,g,kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'
	
	# 2. mCs
	backMIR=mcMIRs['mCs.Background'].intersection(cHI)
	Mh=len(backMIR)
	Nh=len(mcMIRs['mCs.Background'])
	kh=len(mcMIRs['mCs.Significant'].intersection(cHI))
	sh=len(mcMIRs['mCs.Significant'])
	s+='\t'.join([str(x) for x in [hi,'mCs',kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'

write_file('output/MIR_EnrichHistones.txt',s)


##################
### ENCODE TFs ###

encode=[x.split('.')[0] for x in os.listdir('HumanGenome/Encode')]

s=''
for e in encode:
	cE=read_file('HumanGenome/BI.%s.RM.bed' % (e))
	cE=[x.split('\t')[6] for x in cE]
	cE=set([x for x in cE if 'SINE/MIR' in x])
	
	# 1. Genes
	backMIR=geneMIRs['TADsets.Background'].intersection(cE)
	Mh=len(backMIR)
	Nh=len(geneMIRs['TADsets.Background'])
	
	for g in geneSets:
		if g=='TADsets.Background':
			continue
		kh=len(geneMIRs[g].intersection(cE))
		sh=len(geneMIRs[g])
		
		s+='\t'.join([str(x) for x in [e,g,kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'
	
	# 2. mCs
	backMIR=mcMIRs['mCs.Background'].intersection(cE)
	Mh=len(backMIR)
	Nh=len(mcMIRs['mCs.Background'])
	kh=len(mcMIRs['mCs.Significant'].intersection(cE))
	sh=len(mcMIRs['mCs.Significant'])
	s+='\t'.join([str(x) for x in [e,'mCs',kh,sh,Mh,Nh]+hypergeom_test(kh,sh,Mh,Nh)])+'\n'

write_file('output/MIR_EnrichEncode.txt',s)




