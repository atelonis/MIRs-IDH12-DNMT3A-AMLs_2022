
import os
from Functions import *

dSAMPLES={
'IDH':[
	'methylcall.Sample_2268.mincov10.txt',
	'methylcall.Sample_4335.mincov10.txt',
	'methylcall.Sample_5289.mincov10.txt',
	'methylcall.Sample_6887.mincov10.txt',
	'methylcall.Sample_7059.mincov10.txt',
	'methylcall.Sample_7180.mincov10.txt',
	'methylcall.Sample_7301.mincov10.txt',
	'methylcall.Sample_7418.mincov10.txt',
	'methylcall.Sample_7420.mincov10.txt'],
'DNMT3A':[
	'methylcall.Sample_2185.mincov10.txt',
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
	'methylcall.Sample_7411_repeat.mincov10.txt'],
'IDH_DNMT3A':[
	'methylcall.Sample_7067.mincov10.txt',
	'methylcall.Sample_7074.mincov10.txt',
	'methylcall.Sample_7122.mincov10.txt',
	'methylcall.Sample_7145.mincov10.txt',
	'methylcall.Sample_7168.mincov10.txt',
	'methylcall.Sample_7172.mincov10.txt',
	'methylcall.Sample_7316.mincov10.txt',
	'methylcall.Sample_7319.mincov10.txt',
	'methylcall.Sample_7324.mincov10.txt',
	'methylcall.Sample_7328.mincov10.txt',
	'methylcall.Sample_7408.mincov10.txt']}

allFiles=os.listdir('Datasets') # This directory contains the methylcall files downloaded from GEO GSE86952

### Create a dictionary with methylation percentages per sample
d={}
allBases=[]
for dii in dSAMPLES.keys():
	wSamples=dSAMPLES[dii]
	for i in wSamples:
		wF=[x for x in allFiles if i in x][0]
		cName='%s_%s' % (dii,wF.split('.')[1].split('_')[1])
		d[cName]={}
		for j in read_file('Datasets/%s' % (wF))[1:]:
			j=j.split('\t')
			d[cName][j[0]]=j[5]
			allBases+=[j[0]]

dKeys=d.keys()
allBases=sorted(list(set(allBases)))
checks=[ int(float(len(allBases))*float(x)/100.0) for x in [10.0,20.0,40.0,60.0,80.0]]

### Build matrix
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
write_file('Datasets/Dataset.Data.txt',s)
write_file('Datasets/Dataset.Rownames.txt','\n'.join(inBases)+'\n')
write_file('Datasets/Dataset.Colnames.txt','\n'.join(dKeys)+'\n')



