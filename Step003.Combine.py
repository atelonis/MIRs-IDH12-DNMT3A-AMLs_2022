
import os
from Functions import *

### Now, combine the expresison and methylation datasets

ensg2myname={}
for i in read_file('Microarrays/LookupIDs.txt')[1:]:
	i=i.split('\t')
	try:
		ensg2myname[i[3]]+=[i[5]]
	except KeyError:
		ensg2myname[i[3]]=[i[5]]

# Methylation data
d0=read_file('Datasets/Dataset.Data.txt')
rownames0=read_file('Datasets/Dataset.Rownames.txt')
rownames=[x for x in rownames0 if '.NA0' in x]
d=[ [float(z) for z in d0[x].split('\t')] for x,y in enumerate(rownames0) if '.NA0' in y]

# Expression data
e0=read_file('Microarrays/Dataset.50percMostExpr.txt')

rownames+=[x.split('\t')[0].replace('"','') for x in e0[1:]]
d+=[ [float(y) for y in x.split('\t')[1:]] for x in e0[1:] ]

d=numpy.array(d)
numpy.savetxt('Datasets/CombinedData.txt',d,delimiter='\t')
write_file('Datasets/CombinedData.Rownames.txt','\n'.join(rownames)+'\n')


