
import os
import numpy

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

ensgD={}
for i in read_file('MouseGenome/GeneRanges.bed'):
	i=i.split('\t')[3].split('|')
	ensgD[i[0]]='%s|%s' % (i[0],i[1])

mapFiles=[x.replace('.genes.results','') for x in os.listdir('GSE60055/StarMapping/') if 'genes.results' in x]
mapFiles.sort()

d={}
countD={}
allGenes=[]
inSamples=[]
for i in mapFiles:
	d[i]={}
	countD[i]={}
	inSamples+=[i]
	for j in read_file('GSE60055/StarMapping/%s.genes.results' % (i))[1:]:
		j=j.split('\t')
		try:
			j[0]=ensgD[j[0]]
		except KeyError:
			continue
		d[i][j[0]]=j[6]
		countD[i][j[0]]=j[4]
		allGenes+=[j[0]]

inSamples=sorted(inSamples)
allGenes=sorted(list(set(allGenes)))

s='\t'.join(inSamples)+'\n'
m='\t'.join(inSamples)+'\n'
for i in allGenes:
	
	s+=i
	
	newL=[]
	for j in inSamples:
		newL+=[float(d[j][i])]
		s+='\t%s' % (countD[j][i])
	s+='\n'
	
	if numpy.median(newL)>1:
		strToAdd='\t'.join([str(x) for x in newL])
		if strToAdd in m: # for some genes, RSEM calculated the exact same expression (?!)
			print i
			continue
		m+='%s\t%s\n' % (i,strToAdd)

write_file('GSE60055/ExprMatrix.Counts.txt',s)
write_file('GSE60055/ExprMatrix.FPKM_Filter.txt',m)



