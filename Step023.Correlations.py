
import os
import sys
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
	
	if 'GeneBody' in l[2]:
		cDist='0'
	else:
		cDist=l[2].split('_')[4]
	
	r='\t'.join([ l[0], l[1], str(cSpear[0]), str(cSpear[1]), l[2], cDist ])+'\n'
	return r

create_dir('output')
create_dir('output/Correlation')
create_dir('output/Correlation/Step1')
create_dir('output/Correlation/Step2')
create_dir('output/Correlation/Step3')

create_dir('TEMP')

wGenes={}
with open('Datasets/CombinedData.Rownames.txt') as fr:
	with open('Datasets/CombinedData.txt') as fd:
		r=fr.readline()
		d=fd.readline()
		while r:
			r=r.replace('\n','')
			wGenes[r]='hack'
			d=d.replace('\t','\n')
			write_file('TEMP/%s.txt' % (r),d)
			r=fr.readline()
			d=fd.readline()

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
		cC=i[0].split('\t')[0]
		#- Gene
		try:
			cGenes=ensg2myname[i[1].split('|')[0]]
		except KeyError:
			i=f.readline()
			continue
		cGenes.sort(key=lambda x:int(x.split('|')[1]))
		#- Anotation
		cAnno='%s|%s' % (i[2],i[3])
		
		tadPairs+=[ [cC,cGenes[0],cAnno] ]
		i=f.readline()
random.shuffle(tadPairs)

print 'pooling'
p=Pool(117)
s=p.map(corr_1,tadPairs)
write_file('output/Correlation/Step1/TadCorrelations.txt',''.join(s))



