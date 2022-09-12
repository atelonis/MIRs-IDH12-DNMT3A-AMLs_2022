
import os
import numpy
import scipy
from scipy import stats
from scipy.stats import hypergeom

from Functions import *

for i in ['ECO','RED']:
	create_dir('GSE60055/output/%s' % (i))

for mg in ['MEP']:
	for i in ['Double','KO3a','M140']:
		for pn in ['Posit','Negat']:
			print mg,i,pn
			
			cBack=read_file('GSE60055/output/SAM/%s.%s_vs_WT.SAM.IncludedInSAM.txt' % (mg,i))
			cGenes=read_file('GSE60055/output/SAM/%s.%s_vs_WT.SAM.%sSign.txt' % (mg,i,pn))
			cGenes=[x.split('\t')[2].replace('"','') for x in cGenes[1:]]
			
			eco_analysis(cGenes,cBack,'MouseGenome','GSE60055/output/ECO/%s.%s.%s' % (mg,i,pn))
			red_analysis(cGenes,cBack,'MouseGenome','GSE60055/output/RED/%s.%s.%s' % (mg,i,pn))


