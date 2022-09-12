
import os
import sys
import time
import numpy
import scipy
import random
import itertools
from scipy import stats
from multiprocessing import Pool

def read_file(myfile):
	f=open(myfile,'r')
	fl=[x.replace('\n','') for x in f.readlines()]
	f.close()
	return fl

def write_file(MyFile,MyString):
	f=open(MyFile,'w')
	f.write(MyString)
	f.close()
	return

def create_dir(myDIR):
	if not os.path.isdir(myDIR):
		os.mkdir(myDIR)
	return

def get_from_dict(myD,myKey,myKE):
	try:
		return myD[myKey]
	except KeyError:
		return myKE

def corr_1(l): # [ name1, name2]
	
	try:
		l0=[float(x) for x in read_file('TEMP/%s.txt' % (l[0]))]
		l1=[float(x) for x in read_file('TEMP/%s.txt' % (l[1]))]
	except IOError:
		return '\t'.join([ l[0], l[1], 'NA','NA' ])+'\n'
	
	if numpy.std(l0)==0 or numpy.std(l1)==1:
		return '\t'.join([ l[0], l[1], 'NA','NA' ])+'\n'
	
	cSpear=scipy.stats.spearmanr( l0,l1 )
	r='\t'.join([ l[0], l[1], str(cSpear[0]), str(cSpear[1]) ])+'\n'
	return r

create_dir('TEMP')

wCorrs={}
corrPairs=[]
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	wCorrs[i[0]]='hack'
	wCorrs[i[1]]='hack'
	corrPairs+=[ [i[0],i[1]] ]

wGenes={}
with open('WtWtDatasets/CombinedData.Rownames.txt') as fr:
	with open('WtWtDatasets/CombinedData.txt') as fd:
		r=fr.readline()
		d=fd.readline()
		while r:
			r=r.replace('\n','')
			try: # Keep the genes and mCs with significant correlations
				cHack=wCorrs[r]
			except KeyError:
				r=fr.readline()
				d=fd.readline()
				continue
			d=d.replace('\t','\n')
			write_file('TEMP/%s.txt' % (r),d)
			r=fr.readline()
			d=fd.readline()

random.shuffle(corrPairs)

print 'pooling'
p=Pool(117)
s=p.map(corr_1,corrPairs)
write_file('output/Correlation/Step2/WtWtTadCorrelations.txt',''.join(s))



