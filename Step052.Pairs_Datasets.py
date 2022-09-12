
import os
import itertools

from Functions import *


######################################################
### PART 1 | make gene-cytosine pairs to correlate ###
######################################################

strDir={'+':'Upstream','-':'Downstream'}
geneCoord={}
for i in read_file('HumanGenome/GeneRanges.bed'):
	i=i.split('\t')
	geneCoord[i[3]] = [ int(i[1]), int(i[2]), i[5] ]

tadC={}
for i in read_file('TCGA/bed/BedInters.cg.TADs.bed'):
	i=i.split('\t')
	try:
		tadC[i[7]]+=['%s.%s' % (i[3],i[2])]
	except KeyError:
		tadC[i[7]]=['%s.%s' % (i[3],i[2])]

tadG={}
for i in read_file('output/TAD.genes.bed') + read_file('output/TADnot.genes.bed'):
	i=i.split('\t')
	try:
		tadG[i[3]]+=[i[7]]
	except KeyError:
		tadG[i[3]]=[i[7]]

tKeys=list(set(tadG.keys()).intersection(set(tadC.keys())))
print len(tKeys)

with open('TCGA/output/TAD.pairs.txt','w') as f:
	for i in tKeys:
		print tKeys.index(i)
		for c in tadC[i]:
			try:
				cCoord=int(c.split('.')[1])
			except ValueError:
				continue
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



#############################################################
###  PART 2  |  Make expression and methylation datasets  ###
#############################################################

wTCGA=[]
for i in read_file('TCGA/MutationData.txt')[1:]:
	i=i.split('\t')
	if i[4]!='WT':
		continue
	mutStatus=''
	if i[1]!='WT':
		mutStatus+='DNMT3A'
	if i[2]!='WT' or i[3]!='WT':
		mutStatus+='IDH'
	if mutStatus=='DNMT3AIDH':
		mutStatus='Double'
	wTCGA+=['TCGA-AB-%s.%s' % (i[0],mutStatus)]

# Get the samples
d0=read_file('TCGA/LAML.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt')
mrnaFiles=d0[0].split('\t')[1:]

lamlDir='TCGA/LAML.Methylation.Level_3/jhu-usc.edu_LAML.HumanMethylation450.Level_3.2.3.0'
methFiles=[x.split('.')[5][:12] for x in os.listdir(lamlDir) if 'TCGA-AB-' in x]
wIDs=[x for x in methFiles if x in mrnaFiles and x in [y.split('.')[0] for y in wTCGA]]
wTCGA=[x for x in wTCGA if x.split('.')[0] in wIDs]

wDCols=[x for x,y in enumerate(mrnaFiles) if y in wIDs]
write_file('Mutations.txt','\n'.join(wTCGA)+'\n')


# Make mRNA dataset

nArray=[]
geneNames=[]
for i in d0[1:]:
	i=i.split('\t')
	cGene=i[0]
	if cGene.count('|')>1 or '?' in cGene:
		continue
	cGene=cGene.replace('_calculated','')
	cValues=[float(i[x+1]) for x in wDCols]
	nArray+=[cValues]
	geneNames+=[cGene]
	
nArray=numpy.array(nArray)

nMedians = numpy.median(nArray,1)
wMedian  = numpy.median(nMedians)

wIndexes=[x for x,y in enumerate(nMedians) if y>wMedian]

wGenes = [geneNames[x] for x in wIndexes]
wArray = nArray[wIndexes,:]

d='\t'.join([d0[0].split('\t')[x+1] for x in wDCols])+'\n'
for i in range(len(wGenes)):
	d+='%s\t%s\n' % (wGenes[i],'\t'.join([str(x) for x in wArray[i,]]))

write_file('TCGA/RnaData.txt',d)



# Make methylation dataset
toExclude=[]
s=[''] * 485579
for i in wIDs:
	print wIDs.index(i),i
	cFile=[x for x in os.listdir(lamlDir) if i in x]
	with open('%s/%s' % (lamlDir,cFile[0])) as f:
		i=f.readline()
		i=f.readline()
		i=f.readline()
		n=0
		while i:
			i=i.split('\t')
			if s[n]=='':
				s[n]+=i[0]
			s[n]+='\t%s' % (i[1])
			n+=1
			i=f.readline()

s2='%s\n' % ('\t'.join(wIDs))
for i in s:
	i=i.split('\t')
	if 'NA' in i:
		continue
	cVals=[x for x in i[1:] if float(x)>0.3]
	if float(len(cVals)) < float(len(i)-1)*0.03:
		continue
	s2+='\t'.join(i)+'\n'

write_file('TCGA/MethylationData.txt',s2)





