
import os
from Functions import *

def get_distance(a):
	a=a.split('\t')
	cRep= ( int(a[0].split('|')[4]) + int(a[0].split('|')[5]) ) / 2
	cBounds=[ int(a[1].split('_')[2]), int(a[1].split('_')[3]) ]
	cDistance=str(min([abs(x-cRep) for x in cBounds]))
	return '\t'.join(a+[cDistance])


create_dir('output/Correlation/Step4')

tadD={}
for tn in ['TAD','TADnot']:
	with open('output/%s.Cs.bed' % (tn)) as f:
		i=f.readline()
		while i:
			i=i.replace('\n','').split('\t')
			tadD['%s.%s.NA0' % (i[4],i[6])]=i[3]
			i=f.readline()
	with open('output/%s.genes.bed' % (tn)) as f:
		i=f.readline()
		while i:
			i=i.replace('\n','').split('\t')
			tadD[i[7]]=i[3] # but about 6% of the genes are exactly on a TAD boundary
			i=f.readline()

tadReps={}
with open('TADs.Repeats.bed') as f:
	i=f.readline()
	while i:
		i=i.split('\t')
		tadReps[i[7]]=i[3] # For MIRs, only 87 span TAD boundaries
		i=f.readline()


###########
### mCs ###
print 'mCs'

repC={}
with open('bed/BedIntersect.TadCs.Repeats.bed') as f:
	i=f.readline()
	while i:
		i=i.split('\t')
		cCoord=int(i[3].split('.')[1])
		if int(i[7])<cCoord and cCoord<=int(i[8]):
			repC['%s.NA0' % (i[3])]=[i[9]]
		i=f.readline()

for c in ['Background','Significant']:
	print c
	myReps=[]
	for i in read_file('output/Correlation/Step3/Cs.%s.txt' % (c)):
		myReps+=get_from_dict(repC,i,[])
	myReps=sorted(list(set(myReps)))
	myReps=['%s\t%s' % (x,tadReps[x]) for x in myReps]
	myReps=[get_distance(x) for x in myReps]
	write_file('output/Correlation/Step4/Reps.Cs.%s.txt' % (c),'\n'.join(myReps)+'\n')


#############
### GENES ###
print 'genes'

repG={}
with open('HumanGenome/BedInters.Intron.Reps.bed') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		if i[11]!=i[5]:
			i=f.readline()
			continue
		try:
			repG[i[3]]+=[i[9]]
		except KeyError:
			repG[i[3]]=[i[9]]
		i=f.readline()

myGenes=[x for x in os.listdir('output/Correlation/Step3/') if x[-8:]=='ENSG.txt' and x[:7]=='TADsets']

for g in myGenes:
	print g
	myReps=[]
	for i in read_file('output/Correlation/Step3/%s' % (g)):
		myReps+=get_from_dict(repG,i,[])
	myReps=sorted(list(set(myReps)))
	myReps=['%s\t%s' % (x,tadReps[x]) for x in myReps]
	myReps=[get_distance(x) for x in myReps]
	write_file('output/Correlation/Step4/Reps.%s' % (g),'\n'.join(myReps)+'\n')




