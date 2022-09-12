
import os
import numpy
from scipy.stats import hypergeom
from scipy.stats import chisquare

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

for be in ['Beck','Encode']:
	print be
	s='GeneSet\tTF\tZScore\tObserved\tMeanGenespace\tSdGenespace\n'
	TFs=[x.split('.')[2] for x in os.listdir('output/%sMonteCarlo' % (be)) if 'BI.Background' in x]
	for tf in TFs:
		print tf
		for nw in ['Negat','W']:
			for g in ['Overlap','Intermediate','Long']:
				
				# Observed
				obs=read_file('%s/%s_%s.%s.bed' % (be,nw,g,tf))
				obs=float(len(list(set([x.split('\t')[3] for x in obs]))))
				
				# Expected
				genespaceExp=[]
				for i in range(1,1001):
					cGenespace=read_file('output/%sMonteCarlo/Shuffling/BI.Genespace.%s_%s.%s.%s.bed' % (be,nw,g,tf,i))
					cGenespace=list(set([x.split('\t')[3] for x in cGenespace]))
					genespaceExp+=[len(cGenespace)]
				
				gsMean = numpy.mean(genespaceExp)
				gsSD   = numpy.std(genespaceExp)
				
				newS=[ '%s_%s' % (nw,g), tf, round((obs-gsMean)/gsSD,4), obs, gsMean, gsSD]
				s+='\t'.join([str(x) for x in newS])+'\n'
	
	write_file('output/TAD.%sTFs.txt' % (be),s)

