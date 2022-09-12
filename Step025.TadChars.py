
import os
from multiprocessing import Pool

from Functions import *

signCorrs={}
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	i[5]=i[5].split('|')[0]
	signCorrs['correlation|%s--%s' % (i[0],i[1])]=i[5]
	i[0]='cytosine|%s' % (i[0])
	i[1]='gene|%s' % (i[1])
	signCorrs[i[0]]=i[5]
	signCorrs[i[1]]=i[5]

tadCorrs={}
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		cTAD=i[4].split('|')[0]
		try:
			tadCorrs[cTAD]+=['correlation|%s--%s' % (i[0],i[1])]
		except KeyError:
			tadCorrs[cTAD]=['correlation|%s--%s' % (i[0],i[1])]
		tadCorrs[cTAD]+=['cytosine|%s' % (i[0])]
		tadCorrs[cTAD]+=['gene|%s' % (i[1])]
		i=f.readline()

dedm=[]
for pn in ['Posit','Negat']:
	for i in read_file('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.%sSign.txt' % (pn))[1:]:
		i=i.replace('"','').split('\t')
		dedm+=['gene|%s' % (i[2])]

for h in ['Hyper','Hypo']:
	for i in read_file('output/MethylSig/%s.txt' % (h)):
		dedm+=['cytosine|%s.NA0' % (i.split('-')[0])]

dedm=set(dedm)


allReps=list(set(allReps))
s= 'TADspan\t'
s+='BackCorrs\tBackGenes\tBackCytos\t'
s+='SignCorrs\tSignGenes\tSignCytos\t'
s+='DeGenes\tDmCytos\n'
for i in tadCorrs.keys():
	s+='%s\t%s' % (i,int(i.split('_')[3])-int(i.split('_')[2])+1)
	
	# Background
	bTAD=get_from_dict(tadCorrs,i,[])
	bCorrs=[x for x in bTAD if 'correlation' in x]
	bGenes=list(set([x for x in bTAD if 'gene' in x]))
	bCytos=list(set([x for x in bTAD if 'cytosine' in x]))
	s+='\t'+'\t'.join([ str(len(bCorrs)), str(len(bGenes)), str(len(bCytos)) ])
	
	# Significant
	sCorrs=[get_from_dict(signCorrs,x,'NA') for x in bCorrs]
	sGenes=[get_from_dict(signCorrs,x,'NA') for x in bGenes]
	sCytos=[get_from_dict(signCorrs,x,'NA') for x in bCytos]
	sCorrs=[x for x in sCorrs if x==i]
	sGenes=[x for x in sGenes if x==i]
	sCytos=[x for x in sCytos if x==i]
	s+='\t'+'\t'.join([ str(len(sCorrs)), str(len(sGenes)), str(len(sCytos)) ])
	
	# DE Genes / DM cytosines
	s+='\t%s' % (len(list(set(bGenes).intersection(dedm))))
	s+='\t%s' % (len(list(set(bCytos).intersection(dedm))))


write_file('TADs.GenesDeDm.txt',s)




