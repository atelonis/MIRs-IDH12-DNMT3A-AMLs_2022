
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

def corr_1(l): # [ name1, name2, annotation]
	
	l0=[float(x) for x in read_file('TEMP/%s.txt' % (l[0]))]
	l1=[float(x) for x in read_file('TEMP/%s.txt' % (l[1]))]
	cSpear=scipy.stats.spearmanr( l0,l1 )
	r='\t'.join([ l[0], l[1], str(cSpear[0]), str(cSpear[1]), l[2] ])+'\n'
	return r
create_dir('TEMP')

i2cpg={}
cpg2i={}
cpgiBed={}
with open('bed/BedIntersect.Cs.CpGi.bed') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		
		cpg2i['%s.%s.NA0' % (i[0],i[2])]=i[6].replace(':','')
		
		cpgiBed[i[6].replace(':','')]=[float(i[4]),float(i[5])]
		
		try:
			i2cpg[i[6].replace(':','')]+=['%s.%s.NA0' % (i[0],i[2])]
		except KeyError:
			i2cpg[i[6].replace(':','')]=['%s.%s.NA0' % (i[0],i[2])]
		i=f.readline()

meths={}
wGenes={}
with open('Datasets/CombinedData.Rownames.txt') as fr:
	with open('Datasets/CombinedData.txt') as fd:
		r=fr.readline()
		d=fd.readline()
		while r:
			r=r.replace('\n','')
			wGenes[r]='hack'
			if '|' in r: # gene
				d=d.replace('\t','\n')
				write_file('TEMP/%s.txt' % (r),d)
			else: # methylation
				meths[r]=[float(x) for x in d.replace('\n','').split('\t')]
			r=fr.readline()
			d=fd.readline()

noSamples=len(meths[meths.keys()[0]])

for i in i2cpg.keys():
	cCs=i2cpg[i]
	cMeth=[ str(numpy.mean([meths[y][x] for y in cCs])) for x in range(noSamples) ]
	write_file('TEMP/%s.txt' % (i), '\n'.join(cMeth)+'\n')

geneBed={}
with open('HumanGenome/GeneRanges.bed') as f:
	i=f.readline()
	while i:
		i=i.split('\t')
		geneBed[i[3]]=[float(i[1]),float(i[2])]
		i=f.readline()

ensg2myname={}
for i in read_file('Microarrays/LookupIDs.txt')[1:]:
	i=i.split('\t')
	try:
		cHack=wGenes[i[5]]
	except KeyError:
		continue
	try:
		ensg2myname[i[3]]+=[i[5]]
	except KeyError:
		ensg2myname[i[3]]=[i[5]]

for i in ensg2myname.keys():
	ensg2myname[i]=list(set(ensg2myname[i]))

n=0
tadPairs=[]
with open('output/TAD.pairs.txt') as f:
	i=f.readline()
	while i:
		
		n+=1
		if n in [2**x for x in range(30)]:
			print n,len(tadPairs)
		
		i=i.replace('\n','').split('\t')
		#- mC
		try:
			cC=cpg2i[i[0].split('\t')[0]]
		except KeyError:
			i=f.readline()
			continue
		#- Gene
		try:
			cGenes=ensg2myname[i[1].split('|')[0]]
		except KeyError:
			i=f.readline()
			continue
		cGenes.sort(key=lambda x:int(x.split('|')[1]))
		#- Anotation
		iBed=cpgiBed[cC]
		gBed=geneBed[i[1]]
		if iBed[0]<gBed[1]+2000 and iBed[1]>gBed[0]-2000:
			cAnno='Proximal'
		elif iBed[0]<gBed[1]+500000 and iBed[1]>gBed[0]-500000:
			cAnno='Intermediate'
		else:
			cAnno='Long'
		
		tadPairs+=['%s--%s--%s' % (cC,cGenes[0],cAnno) ]
		i=f.readline()
tadPairs=list(set(tadPairs))
tadPairs=[ [x.split('--')[0],x.split('--')[1],x.split('--')[2]] for x in tadPairs]
random.shuffle(tadPairs)

print 'pooling'
p=Pool(117)
s=p.map(corr_1,tadPairs)
write_file('output/Correlation/Step1/CpGiCorrelations.txt',''.join(s))


