
import os
import sys
import numpy

from Functions import *

for hm in ['Human','Mouse']:
	allReps = read_file('%sGenome/RepeatMasker.bed' % (hm))
	allReps = sorted(list(set([x.split('\t')[3].split('|')[0] for x in allReps])))
	
	write_file('HumanGenome/AllRepeats.Family.txt', '\n'.join(allReps)+'\n' )
	
	EIs=['Exon','Intron']
	
	dRanges={}
	allGenes=[]
	for ei in EIs:
		dRanges[ei]={}
		for i in read_file('%sGenome/%sRanges.bed' % (hm,ei)):
			i=i.split('\t')
			try:
				dRanges[ei][i[3]]+=int(i[2])-int(i[1])
			except KeyError:
				dRanges[ei][i[3]]=int(i[2])-int(i[1])
			allGenes+=[i[3]]
	allGenes=sorted(list(set(allGenes)))
	
	for ei in EIs:
		print ei
		
		dRange={}
		dRange.clear()
		with open('%sGenome/%sRanges.bed' % (hm,ei)) as f:
			i=f.readline()
			while i:
				i=i.split('\t')
				try:
					dRange[i[3]]+=int(i[2])+int(i[1])+1
				except KeyError:
					dRange[i[3]]=int(i[2])+int(i[1])+1
				i=f.readline()
		
		# Define matrix
		famSense=numpy.array([[0]*len(allReps)]*len(allGenes))
		famAntis=numpy.array([[0]*len(allReps)]*len(allGenes))
		
		with open('bed/%s.BedInters.%s.Reps.bed' % (hm,ei)) as f:
			i=f.readline()
			while i:
				i=i.replace('\n','').split('\t')
				
				cRep=i[9].split('|')[0]
				cSubRep='%s|%s' % (cRep,i[9].split('|')[1])
				if cRep not in allReps:
					i=f.readline()
					continue
				
				cOverlap = min([int(i[2]),int(i[8])]) - max([int(i[1]),int(i[7])])
				
				if i[5]==i[11]:
					famSense[ allGenes.index(i[3]), allReps.index(cRep) ] += cOverlap
				else:
					famAntis[ allGenes.index(i[3]), allReps.index(cRep) ] += cOverlap
					cOrient='Antisense'
				i=f.readline()
		
		# Make strings and write file
		fCols='GeneBases\t'+'\t'.join(allReps)
		rownames=['%s\t%s' % (x,get_from_dict(dRanges[ei],x,0)) for x in allGenes]
		
		save_rep_matrix( famSense, rownames, fCols, '%sGenome/Number_Family.%s.Sense.txt' % (hm,ei) ) # Family sense
		save_rep_matrix( famAntis, rownames, fCols, '%sGenome/Number_Family.%s.Antisense.txt' % (hm,ei) ) # Family antisense
		
	for eip in ['Exon','Intron']:
		for sa in ['Sense','Antisense']:
			print eip,sa
	
			cTable=read_file('Number_Family.%s.%s.txt' % (eip,sa))
			
			s='\t'.join(cTable[0].split('\t')[1:])+'\n'
			for i in cTable[1:]:
				i=i.split('\t')
				s+=i[0]
				for j in i[2:]:
					try:
						s+='\t%s' % (float(j)/float(i[1]))
					except ZeroDivisionError:
						s+='\t0'
				s+='\n'
			write_file('%sGenome/Density_Family.%s.%s.txt' % (hm,eip,sa),s)

