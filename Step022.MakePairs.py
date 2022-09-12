
import os
import itertools
from Functions import *

strDir={'+':'Upstream','-':'Downstream'}
geneCoord={}
for i in read_file('HumanGenome/GeneRanges.bed'):
	i=i.split('\t')
	geneCoord[i[3]] = [ int(i[1]), int(i[2]), i[5] ]

tadC={}
for i in read_file('output/TAD.Cs.bed')+read_file('output/TADnot.Cs.bed'):
	i=i.split('\t')
	try:
		tadC[i[3]]+=['%s.%s.NA0' % (i[0],i[6])]
	except KeyError:
		tadC[i[3]]=['%s.%s.NA0' % (i[0],i[6])]

tadG={}
for i in read_file('output/TAD.genes.bed')+read_file('output/TADnot.genes.bed'):
	i=i.split('\t')
	try:
		tadG[i[3]]+=[i[7]]
	except KeyError:
		tadG[i[3]]=[i[7]]

tKeys=list(set(tadG.keys()).intersection(set(tadC.keys())))
print len(tKeys)

with open('output/TAD.pairs.txt','w') as f:
	for i in tKeys:
		print tKeys.index(i)
		for c in tadC[i]:
			cCoord=int(c.split('.')[1])
			for g in tadG[i]:
				
				# Get position with respect to gene
				if cCoord >= geneCoord[g][0] and cCoord <= geneCoord[g][1]:
					cPosition='GeneBody'
				elif cCoord < geneCoord[g][0]:
					if geneCoord[g][2]=='+':
						cPosition='Upstream_%s' % (geneCoord[g][0]-cCoord)
					else:
						cPosition='Downstream_%s' % (geneCoord[g][0]-cCoord)
				elif cCoord > geneCoord[g][1]:
					if geneCoord[g][2]=='+':
						cPosition='Downstream_%s' % (cCoord-geneCoord[g][1])
					else:
						cPosition='Upstream_%s' % (cCoord-geneCoord[g][1])
				else:
					print g,c,geneCoord[g]
					opa # Nothing gets here, good!
				
				f.write('\t'.join([c,g,i,cPosition])+'\n')


