
import os
from Functions import *
from multiprocessing import Pool


gID={}
with open('Microarrays/LookupIDs.txt') as f:
	i=f.readline()
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		gID[i[5]]='%s|%s' % (i[3],i[1])
		gID[i[3]]=i[5]
		i=f.readline()

samGenes={}
for i in read_file('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.IncludedInSAM.txt'):
	samGenes[gID[i]]='hack'

tss={}
### File downloaded from Biomart of ensembl.org with fields:
###		Chromosome/scaffold name
### 	Strand
### 	Transcription start site (TSS)
### 	Transcript stable ID version
### 	Gene stable ID version
### 	Gene name
for i in read_file('HumanGenome/mart_export.TSS.txt')[1:]:
	i=i.split('\t')
	cGene='%s|%s' % (i[4].split('.')[0],i[5])
	try:
		cHack=samGenes[cGene]
	except KeyError:
		continue
	try:
		tss['chr%s' % (i[0])]+=['%s|%s' % (i[2],cGene)]
	except KeyError:
		tss['chr%s' % (i[0])]=['%s|%s' % (i[2],cGene)]

cytos=[]
with open('Datasets/Dataset.Rownames.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('.')
		if i[2]!='NA0':
			i=f.readline()
			continue
		cytos+=['%s.%s.NA0' % (i[0],i[1])]
		i=f.readline()

########################
def process_one(myCyto):
	cPos=int(myCyto.split('.')[1])
	try:
		cTSS=tss[myCyto.split('.')[0]]
	except KeyError:
		return '\t'.join([ myCyto, 'NA', 'NA', 'NA', 'NA' ]) + '\n'
	cTSS.sort(key=lambda x:abs(int(x.split('|')[0])-cPos))
	cTSS=cTSS[0].split('|')
	return '%s\t%s\t%s\t%s|%s\t%s\n' % (myCyto, cTSS[0], abs(int(cTSS[0])-cPos), cTSS[1], cTSS[2], gID[cTSS[1]])
########################

p=Pool(150)
s=p.map(process_one,cytos)
write_file('output/ClosestTSS.txt','cytosine\tTSSpos\tDistance\tGeneID1\tGeneID2\n'+''.join(s))


