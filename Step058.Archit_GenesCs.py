
import os
import numpy
import scipy
from scipy import stats
from scipy.stats import hypergeom
from Functions import *


#####################################
###  PART 1  |  Gene architecture ###
#####################################

for i in ['ECO','RED']:
	create_dir('TCGA/output/%s' % (i))

backGenes=read_file('TCGA/output/DAVID/Background.TADs_WholeNames.txt')

e2f={}
for i in backGenes:
	e2f[i.split('|')[0]]=i

tadSets=[x.replace('.txt','') for x in os.listdir('TCGA/output/Correlation/Step3/') if x.split('.')[0]=='TADsets']
tadSets=[x for x in tadSets if x.split('.')[0]=='TADsets' and x!='TADsets.GenesCsInTAD' and 'ENSG' not in x]
for i in tadSets:
	cGenes=read_file('TCGA/output/Correlation/Step3/%s.txt' % (i))
	if len(cGenes)<50:
		continue
	eco_analysis(cGenes,backGenes,'HumanGenome','TCGA/output/ECO/%s' % (i))
	red_analysis(cGenes,backGenes,'HumanGenome','TCGA/output/RED/%s' % (i))

wSets=[x.split('.')[1] for x in os.listdir('TCGA/output/DAVID') if 'TAD_W' in x]
for i in wSets:
	cGenes=read_file('TCGA/output/DAVID/TAD_W.%s/Genes.txt' % (i))
	cGenes=[e2f[x] for x in cGenes]
	write_file('TCGA/output/Correlation/Step3/TADsets_W.%s.txt' % (i),'\n'.join(cGenes)+'\n')
	
	eco_analysis(cGenes,backGenes,'HumanGenome','TCGA/output/ECO/TadW_%s' % (i))
	red_analysis(cGenes,backGenes,'HumanGenome','TCGA/output/RED/TadW_%s' % (i))



##########################################################################################
###  PART 2  |  Are genes correlated with the same or different cytosines in the TAD?  ###
##########################################################################################


signCorrs={}
for i in read_file('TCGA/output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	signCorrs['%s--%s' % (i[0],i[1])]=i[0]

G={}
tadG={}
tadC={}
with open('TCGA/output/Correlation/Step1/Correlations.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		try:
			G[i[1]]+=[i[0]]
		except KeyError:
			G[i[1]]=[i[0]]
		i[4]=i[4].split('|')[0]
		try:
			tadG[i[4]]+=[i[1]]
		except KeyError:
			tadG[i[4]]=[i[1]]
		try:
			tadC[i[4]]+=[i[0]]
		except KeyError:
			tadC[i[4]]=[i[0]]
		i=f.readline()

tKeys=tadG.keys()

with open('TCGA/output/Correlation/Step3/TADsets.GenesCsInTAD.txt','w') as f:
	f.write('TAD\tTadCs\tGene1\tGene2\tGene1Cs\tGene2Cs\tComCs\n')
	for t in tKeys:
		print t,tKeys.index(t),len(tKeys)
	
		tc=set(tadC[t])
		allGenes=list(set(tadG[t]))
		for i in range(len(allGenes)-1):
			gi=allGenes[i]
			ci=set(G[gi]).intersection(tc)
			sci=[get_from_dict(signCorrs,'%s--%s' % (x,gi),'NA') for x in ci]
			sci=set([x for x in sci if x!='NA'])
			for j in range(i+1,len(allGenes)):
				gj=allGenes[j]
				cj=set(G[gj]).intersection(tc)
				scj=[get_from_dict(signCorrs,'%s--%s' % (x,gj),'NA') for x in cj]
				scj=set([x for x in scj if x!='NA'])
				s=[t,len(list(tc)),gi,gj,len(list(sci)),len(list(scj)),len(list(sci.intersection(scj)))]
				s='\t'.join([str(x) for x in s])+'\n'
				f.write(s)






