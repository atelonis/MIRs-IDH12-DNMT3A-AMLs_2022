
import os
import numpy
import scipy
from scipy import stats
from scipy.stats import hypergeom
from Functions import *

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

for i in ['ECO','RED']:
	create_dir('output/%s' % (i))

e2m={} # ENSG to myName
for i in read_file('HumanGenome/GeneRanges.bed'):
	i=i.split('\t')[3]
	e2m[i.split('|')[0]]=i

entrez2myname={}
n2m={} # name in microarrays  to   myName
for i in read_file('HumanGenome/LookupIDs.txt')[1:]:
	i=i.split('\t')
	try:
		n2m[i[5]]=e2m[i[3]]
	except KeyError:
		pass
	try:
		entrez2myname[i[5].split('|')[1]]=e2m[i[3]]
	except KeyError:
		pass

backGenes=[]
entrez2affy={}
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.split('\t')[1]
		backGenes+=[i]
		entrez2affy[i.split('|')[1]]=i
		i=f.readline()

backGenes=[get_from_dict(n2m,x,'NA') for x in list(set(backGenes))]
backGenes=[x for x in backGenes if x!='NA']
write_file('output/Correlation/Step3/TADsets.Background.ENSG.txt','\n'.join(backGenes)+'\n')

tadSets=[x.replace('.txt','') for x in os.listdir('output/Correlation/Step3/') if x.split('.')[0]=='TADsets']
tadSets=[x for x in tadSets if x.split('.')[0]=='TADsets' and x!='TADsets.GenesCsInTAD' and 'ENSG' not in x and 'Names' not in x]
for i in tadSets:
	print i
	cGenes=read_file('output/Correlation/Step3/%s.txt' % (i))
	cGenes=[get_from_dict(n2m,x,'NA') for x in list(set(cGenes))]
	cGenes=[x for x in cGenes if x!='NA']
	if len(cGenes)<50:
		continue
	write_file('output/Correlation/Step3/%s.ENSG.txt' % (i),'\n'.join(cGenes)+'\n')
	write_file('output/Correlation/Step3/%s.Names.txt' % (i),'\n'.join(list(set([x.split('|')[1] for x in cGenes])))+'\n')
	
	eco_analysis(cGenes,backGenes,'HumanGenome','output/ECO/%s' % (i))
	red_analysis(cGenes,backGenes,'HumanGenome','output/RED/%s' % (i))

wSets=[x.split('.')[1] for x in os.listdir('output/DAVID') if 'TAD_W' in x]
for i in wSets:
	print i
	cGenes=read_file('output/DAVID/TAD_W.%s/Genes.txt' % (i))
	write_file('output/Correlation/Step3/TADsets_W.%s.txt' % (i),'\n'.join([entrez2affy[x] for x in cGenes])+'\n')
	
	cGenes=[get_from_dict(entrez2myname,x,'NA') for x in list(set(cGenes))]
	cGenes=[x for x in cGenes if x!='NA']
	write_file('output/Correlation/Step3/TADsets_W.%s.ENSG.txt' % (i),'\n'.join(cGenes)+'\n')
	
	eco_analysis(cGenes,backGenes,'HumanGenome','output/ECO/TadW_%s' % (i))
	red_analysis(cGenes,backGenes,'HumanGenome','output/RED/%s' % (i))


