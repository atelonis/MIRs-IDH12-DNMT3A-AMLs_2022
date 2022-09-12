
import os
import time
import numpy

def read_file(myfile):
	f=open(myfile,'r')
	fl=[x.replace('\n','').replace('\r','') for x in f.readlines()]
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
		return myD[myK]
	except KeyError:
		return myKE

create_dir('WtWtDatasets')

idSAMPLES=[
	'methylcall.Sample_2268.mincov10.txt', ### IDH
	'methylcall.Sample_4335.mincov10.txt',
	'methylcall.Sample_5289.mincov10.txt',
	'methylcall.Sample_6887.mincov10.txt',
	'methylcall.Sample_7059.mincov10.txt',
	'methylcall.Sample_7180.mincov10.txt',
	'methylcall.Sample_7301.mincov10.txt',
	'methylcall.Sample_7418.mincov10.txt',
	'methylcall.Sample_7420.mincov10.txt',
	'methylcall.Sample_2185.mincov10.txt', ### DNMT3A
	'methylcall.Sample_2195.mincov10.txt',
	'methylcall.Sample_2206.mincov10.txt',
	'methylcall.Sample_2234.mincov10.txt',
	'methylcall.Sample_3331.mincov10.txt',
	'methylcall.Sample_4334.mincov10.txt',
	'methylcall.Sample_5286.mincov10.txt',
	'methylcall.Sample_6239.mincov10.txt',
	'methylcall.Sample_6374.mincov10.txt',
	'methylcall.Sample_7131.mincov10.txt',
	'methylcall.Sample_7188.mincov10.txt',
	'methylcall.Sample_7313.mincov10.txt',
	'methylcall.Sample_7318.mincov10.txt',
	'methylcall.Sample_7322.mincov10.txt',
	'methylcall.Sample_7407.mincov10.txt',
	'methylcall.Sample_7411_repeat.mincov10.txt',
	'methylcall.Sample_7067.mincov10.txt', ### Double mutants
	'methylcall.Sample_7074.mincov10.txt',
	'methylcall.Sample_7122.mincov10.txt',
	'methylcall.Sample_7145.mincov10.txt',
	'methylcall.Sample_7168.mincov10.txt',
	'methylcall.Sample_7172.mincov10.txt',
	'methylcall.Sample_7316.mincov10.txt',
	'methylcall.Sample_7319.mincov10.txt',
	'methylcall.Sample_7324.mincov10.txt',
	'methylcall.Sample_7328.mincov10.txt',
	'methylcall.Sample_7408.mincov10.txt']

allFiles=[x for x in os.listdir('Datasets') if 'mincov10.txt' in x]
idFiles=[x for x in allFiles if '_'.join(x.split('_')[1:]) in idSAMPLES]
wtwtSAMPLES=[x for x in allFiles if x not in idFiles]

d={}
allBases=[]
for i in wtwtSAMPLES:
	print i
	t0=time.time()
	cName=i.split('.')[1]
	d[cName]={}
	for j in read_file('Datasets/%s' % (i))[1:]:
		j=j.split('\t')
		d[cName][j[0]]=j[5]
		allBases+=[j[0]]
	print time.time()-t0

dKeys=d.keys()
allBases=sorted(list(set(allBases)))
checks=[ int(float(len(allBases))*float(x)/100.0) for x in [10.0,20.0,40.0,60.0,80.0]]

counter=0
t0=time.time()
print 'Step 2...'

s=''
inBases=[]
for i in allBases:
	counter+=1
	if counter in checks:
		print '%s out of %s completed in %s sec' % (counter,len(allBases),round(time.time()-t0))
	n=[]
	for j in dKeys:
		try:
			n+=[ d[j][i] ]
		except KeyError:
			n+=[ 'NA' ]
	NAs=n.count('NA')
	if NAs > float(len(n))/2.0:
		continue
	inBases+=[ '%s.NA%s' % (i,NAs) ]
	s+='\t'.join(n)+'\n'
write_file('WtWtDatasets/Dataset.Data.txt',s)
write_file('WtWtDatasets/Dataset.Rownames.txt','\n'.join(inBases)+'\n')
write_file('WtWtDatasets/Dataset.Colnames.txt','\n'.join(dKeys)+'\n')

#### Run Step004 at this point
#### Now, combine the datasets

ensg2myname={}
for i in read_file('WtWtDatasets/MicroarraystLookupIDs.txt')[1:]:
	i=i.split('\t')
	try:
		ensg2myname[i[3]]+=[i[5]]
	except KeyError:
		ensg2myname[i[3]]=[i[5]]

# Methylation data
d0=read_file('WtWtDatasets/Dataset.Data.txt')
rownames0=read_file('WtWtDatasets/Dataset.Rownames.txt')
rownames=rownames0
d=[ [z for z in d0[x].split('\t')] for x,y in enumerate(rownames)]

# Expression data
e0=read_file('WtWtDatasets/Expression.txt')

rownames+=[x.split('\t')[0].replace('"','') for x in e0[1:]]
d='\n'.join(d0)+'\n'
d+=''.join([ '\t'.join([y for y in x.split('\t')[1:]])+'\n' for x in e0[1:] ])

write_file('WtWtDatasets/CombinedData.txt',d)
write_file('WtWtDatasets/CombinedData.Rownames.txt','\n'.join(rownames)+'\n')







