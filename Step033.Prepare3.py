
import os
import numpy

for hm in['Human','Mouse']:
	allGenes=sorted(list(set([x.split('\t')[3] for x in read_file('%sGenome/GeneRanges.bed' % (hm))])))
	
	dRanges={}
	dCounts={}
	GCcontent={}
	EvolCons={}
	dRanges.clear()
	dCounts.clear()
	GCcontent.clear()
	EvolCons.clear()
	for i in ['Exon','Intron']:
		print i
		dRanges[i]   = get_ranges(i)
		dCounts[i]   = get_counts(i)
		GCcontent[i] = get_gc_content(i)
		EvolCons[i]  = get_evol_cons(i)
	
	s ='GeneName\tExonLength\tIntronLength\tExonContent\tExonCount\t'
	s+='GC_Exon\tGC_Intron\tExonEvolCons\tIntronEvolCons\n'
	for i in allGenes:
		s+=i
		s+='\t%s' % (numpy.log2(dRanges['Exon'][i]))
		s+='\t%s' % (numpy.log2(get_from_dict(dRanges['Intron'],i,1)))
		s+='\t%s' % (float(dRanges['Exon'][i]) / (float(get_from_dict(dRanges['Intron'],i,0))+float(dRanges['Exon'][i])))
		s+='\t%s' % (dCounts['Exon'][i])
		
		s+='\t%s' % (GCcontent['Exon'][i])
		s+='\t%s' % (get_from_dict(GCcontent['Intron'],i,'NA'))
		
		s+='\t%s' % (get_from_dict(EvolCons['Exon'],i,'NA'))
		s+='\t%s' % (get_from_dict(EvolCons['Intron'],i,'NA'))
		s+='\n'
	
	write_file('%sGenome/GeneParameters.txt' % (hm),s)


