
import os
import sys
import numpy
import scipy
from scipy import stats
from scipy.stats import rankdata
from scipy.stats import hypergeom

from Functions import *


###########################
### PART 1 | Conver IDs ###
###########################

dNames={}
for i in read_file('HumanGenome/GeneRanges.bed'):
	i=i.split('\t')[3]
	dNames[i.split('|')[1]]=i

rz='ZScored_RMA'
create_dir('Hematopoiesis/output/%s' % (rz))	
samComps=list(set([x.split('.')[0] for x in os.listdir('Hematopoiesis/output/SAM/') if rz==x[:len(rz)]]))

s=[]
for i in samComps:
	cGenes0=read_file('Hematopoiesis/output/SAM/%s.SAM.IncludedInSAM.txt' % (i))
	cGenes1=[x.split('|')[0] for x in cGenes0]
	cGenes=[get_from_dict(dNames,x,'NA') for x in cGenes1]
	s+=['%s\t%s\n' % (cGenes0[x],cGenes[x]) for x in range(len(cGenes0))]
	cGenes=list(set([x for x in cGenes if x!='NA']))
	write_file('Hematopoiesis/output/%s/%s.Background.txt' % (rz,i),'\n'.join(cGenes)+'\n')
	for pn in ['Posit','Negat']:
		print i,pn
		cGenes=read_file('Hematopoiesis/output/SAM/%s.SAM.%sSign.txt' % (i,pn))
		cGenes=[x.split('\t')[2].replace('"','').split('|')[0] for x in cGenes[1:]]
		cGenes=[get_from_dict(dNames,x,'NA') for x in cGenes]
		cGenes=list(set([x for x in cGenes if x!='NA']))
		write_file('Hematopoiesis/output/%s/%s.%s.txt' % (rz,i,pn),'\n'.join(cGenes)+'\n')

s=sorted(list(set(s)))
s=[x for x in s if x.split('\t')[1]!='NA\n']
write_file('Hematopoiesis/ID_LookupTable.%s.txt' % (rz),''.join(s))



###################################
### PART 2 | Check architecture ###
###################################

for i in ['ECO','RED']:
	create_dir('Hematopoesis/output/%s' % (i))

myComps=sorted(list(set([x.split('.')[0] for x in os.listdir('Hematopoiesis/output/ZScored_RMA/')])))

for i in myComps:
	cBack=read_file('output/ZScored_RMA/%s.Background.txt' % (i))
	for pn in ['Posit','Negat']:
		cGenes=read_file('output/ZScored_RMA/%s.%s.txt' % (i,pn))
		
		eco_analysis(cGenes,cBack,'HumanGenome','Hematopoiesis/output/ECO/%s.%s' % (i,pn))
		red_analysis(cGenes,cBack,'HumanGenome','Hematopoiesis/output/RED/%s.%s' % (i,pn))



#######################################################
### PART 3 | Overlaps of DE genes with AML/DE genes ###
#######################################################

def calculate_FDR(p_vals):
	ranked_p_values = rankdata(p_vals)
	fdr = numpy.array(p_vals) * len(p_vals) / ranked_p_values
	fdr[fdr > 1] = 1
	return fdr	

def check_enrichments(myGenes,myOutput):
	comps=sorted(list(set([x.split('_')[2] for x in os.listdir('Hematopoiesis/output/SAM/') if 'ZScored_RMA'==x[:11]])))
	
	s=[]
	for i in comps:
		for pn in ['Posit','Negat']:
			inSAM=read_file('Hematopoiesis/output/SAM/ZScored_RMA_%s_vs_HSC.SAM.IncludedInSAM.txt' % (i))
			de=read_file('Hematopoiesis/output/SAM/ZScored_RMA_%s_vs_HSC.SAM.%sSign.txt' % (i,pn))[1:]
			de=[x.split('\t')[2].replace('"','') for x in de]
	
			mySAM = float(len(list(set(inSAM).intersection(myGenes))))
			myDE  = float(len(list(set(de).intersection(myGenes))))
			de    = float(len(de))
			inSAM = float(len(inSAM))
			
			FE = (myDE/de)/(mySAM/inSAM)
			
			if FE > 1:
				pVal=1-hypergeom.cdf( myDE-1, inSAM, mySAM, de )
			else:
				pVal=hypergeom.cdf( myDE, inSAM, mySAM, de )
			
			ipn='%s.%s' % (i,pn)
			s+=[ '\t'.join([str(x) for x in [ipn,myDE,de,mySAM,inSAM,numpy.log2(FE),pVal]]) ]
	
	pVALS = [float(x.split('\t')[-1]) for x in s]
	FDR   = calculate_FDR(pVALS)
	
	s=''.join(['%s\t%s\n' % (s[x],FDR[x]) for x in range(len(s))])
	
	s='Cell\tDE_IQA\tDE\tInIQA\tInSAM\tFE\tPValue\tFDR\n'+s
	
	write_file('Hematopoiesis/output/DE.EnrichedInAmlCorrs.txt',s)
	write_file(myOutput,s)
	
	return

create_dir('Hematopoiesis/output/AmlAnalyses')
samDir='output/SAM'
quant3=list(set([x.split('.')[1] for x in os.listdir(samDir) if 'Quantile_3'==x.split('.')[0]]))

for i in quant3:
	for pn in ['Posit','Negat']:
		cSign=read_file('%s/Quantile_3.%s.SAM.%sSign.txt' % (samDir,i,pn))[1:]
		cSign=[x.split('\t')[2].replace('"','') for x in cSign]
		if len(cSign)<10:
			continue
		check_enrichments(cSign,'Hematopoiesis/output/AmlAnalyses/%s.%s.txt' % (i,pn))






