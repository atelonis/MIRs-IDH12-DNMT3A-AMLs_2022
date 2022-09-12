
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

create_dir('NonMutantDatasets')

exclTCGA=[]
for i in read_file('TCGA/MutationData.txt')[1:]:
	i=i.split('\t')
	exclTCGA+=['TCGA-AB-%s' % (i[0])]

# Get the samples
d0=read_file('TCGA/LAML.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt')
mrnaFiles=d0[0].split('\t')[1:]

lamlDir='TCGA/LAML.HumanMethylation450.Level_3/jhu-usc.edu_LAML.HumanMethylation450.Level_3.2.3.0'
methFiles=[x.split('.')[5][:12] for x in os.listdir(lamlDir) if 'TCGA-AB-' in x]
wIDs=[x for x in methFiles if x in mrnaFiles and x not in exclTCGA]
wDCols=[x for x,y in enumerate(mrnaFiles) if y in wIDs]


# Make mRNA dataset

nArray=[]
geneNames=[]
for i in d0[1:]:
	i=i.split('\t')
	cGene=i[0]
	if cGene.count('|')>1 or '?' in cGene:
		continue
	cGene=cGene.replace('_calculated','')
	cValues=[float(i[x+1]) for x in wDCols]
	nArray+=[cValues]
	geneNames+=[cGene]
	
nArray=numpy.array(nArray)

nMedians = numpy.median(nArray,1)
wMedian  = numpy.median(nMedians)

wIndexes=[x for x,y in enumerate(nMedians) if y>wMedian]

wGenes = [geneNames[x] for x in wIndexes]
wArray = nArray[wIndexes,:]

d='\t'.join([d0[0].split('\t')[x+1] for x in wDCols])+'\n'
for i in range(len(wGenes)):
	d+='%s\t%s\n' % (wGenes[i],'\t'.join([str(x) for x in wArray[i,]]))

write_file('NonMutantDatasets/RnaData.txt',d)



# Make methylation dataset
toExclude=[]
s=[''] * 485579
for i in wIDs:
	print wIDs.index(i),i
	cFile=[x for x in os.listdir(lamlDir) if i in x]
	with open('%s/%s' % (lamlDir,cFile[0])) as f:
		i=f.readline()
		i=f.readline()
		i=f.readline()
		n=0
		while i:
			i=i.split('\t')
			if s[n]=='':
				s[n]+=i[0]
			s[n]+='\t%s' % (i[1])
			n+=1
			i=f.readline()

s2='%s\n' % ('\t'.join(wIDs))
for i in s:
	i=i.split('\t')
	if 'NA' in i:
		continue
	cVals=[x for x in i[1:] if float(x)>0.3]
	if float(len(cVals)) < float(len(i)-1)*0.03:
		continue
	s2+='\t'.join(i)+'\n'

write_file('NonMutantDatasets/MethylationData.txt',s2)





