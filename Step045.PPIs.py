
import os
from multiprocessing import Pool
from scipy.stats import hypergeom
from Functions import *

ppi={}
for i in read_file('Pickle/v3.3/ENSG.GeneNormalized.txt')[1:]:
	i=i.split('\t')
	i0=i[0].split('|')[1]
	i1=i[1].split('|')[1]
	try:
		ppi[i0]+=[i1]
	except KeyError:
		ppi[i0]=[i1]
	try:
		ppi[i1]+=[i0]
	except KeyError:
		ppi[i1]=[i0]

geneSets=['Negat_Overlap','Negat_Intermediate','Negat_Long','W_Overlap','W_Intermediate','W_Long']

beTFs={}
for be in ['Beck','Encode']:
	beTFs[be]={}
	
	beTFs[be]['mCs']={}
	beTFs[be]['Gene']={}
	TFs=[x.split('.')[2] for x in os.listdir(be) if x[:5]=='BI.Cs']
	
	for tf in TFs:
		print be,tf
		### cytosines
		with open('%s/BI.Cs.%s.bed' % (be,tf)) as f:
			i=f.readline()
			while i:
				i=i.replace('\n','').split('\t')
				cC='%s.%s.NA0' % (i[0],i[2])
				try:
					beTFs[be]['mCs'][cC]+=[tf]
				except KeyError:
					beTFs[be]['mCs'][cC]=[tf]
				i=f.readline()
		
		### Genes
		for g in geneSets:
			with open('%s/%s.%s.bed' % (be,g,tf)) as f:
				i=f.readline()
				while i:
					i=i.split('\t')
					i3=i[3].split('|')[1]
					try:
						beTFs[be]['Gene'][i3]+=[tf]
					except KeyError:
						beTFs[be]['Gene'][i3]=[tf]
					i=f.readline()

print 2

######################
def list_to_string(a):
	if a==[]:
		return '\t-'
	else:
		return '\t%s' % (','.join(a))
#####################################

checks=[2**x for x in range(50)]

n=0
with open('output/TAD.PPI_Pairs.txt','w') as fw:
	with open('output/TAD.pairs.txt') as f:
		i=f.readline()
		while i:
			n+=1
			if n in checks:
				print n
			i=i.replace('\n','').split('\t')
			s='%s\t%s' % (i[0],i[1])
			for be in ['Beck','Encode']:
				mcTF = list(set(get_from_dict( beTFs[be]['mCs'], i[0], [] )))
				geneTF = set(get_from_dict( beTFs[be]['Gene'], i[1].split('|')[1], [] ))
				
				cPPI=[]
				for j in mcTF:
					jPPI=set(get_from_dict(ppi,j,[]))
					for k in list(jPPI.intersection(geneTF)):
						cPPI+=['--'.join([j,k]) ] #  !! unsorted, so with localization
				
				# PPIs
				s+=list_to_string(list(geneTF))
				s+=list_to_string(mcTF)
				s+=list_to_string(list(set(cPPI)))
			
			s+='\n'
			fw.writelines(s)
			i=f.readline()


################
###  PART 2  ###
################

create_dir('TEMP')

ensg2myname={}
for i in read_file('Microarrays/LookupIDs.txt')[1:]:
	i=i.split('\t')
	ensg2myname[i[3]]=[i[5]]

for i in ensg2myname.keys():
	ensg2myname[i]=sorted(list(set(ensg2myname[i])))[0]

n=0
checks=[2**x for x in range(40)]

print 1

signCorrs={'Overlap':[], 'Intermediate':[], 'Long':[]}
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	signCorrs[i[7]]+=['%s--%s' % (i[0],i[1])]
signCorrs['TADwide'] = signCorrs['Overlap'] + signCorrs['Intermediate'] + signCorrs['Long']

for i in signCorrs.keys():
	signCorrs[i]=set(signCorrs[i])

allCorrs={'Overlap':[], 'Intermediate':[], 'Long':[]}
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		if i[5]=='0':
			cBin='Overlap'
		elif int(i[5])<500000:
			cBin='Intermediate'
		else:
			cBin='Long'
		allCorrs[cBin]+=['%s--%s' % (i[0],i[1])]
		i=f.readline()
allCorrs['TADwide'] = allCorrs['Overlap'] + allCorrs['Intermediate'] + allCorrs['Long']

for i in allCorrs.keys():
	allCorrs[i]=set(allCorrs[i])

########################
def get_one_ppi(myPPI):
	
	cPPIs=set(read_file('TEMP/%s.txt' % (myPPI)))
	
	hN=len(allCorrs)
	
	r=''
	for si in signCorrs.keys():
		hs=len(cPPIs.intersection(allCorrs[si]))
		hk=len(cPPIs.intersection(signCorrs[si]))
		hM=len(signCorrs[si])
		hN=len(allCorrs[si])
		
		r+='\t'.join([str(x) for x in [myPPI,si,hk,hs,hM,hN]+hypergeom_test(hk,hs,hM,hN)])+'\n'
	
	return r
#############

print 2

PPIs=[x.replace('.txt','') for x in os.listdir('TEMP')]

p=Pool(117)
s=p.map(get_one_ppi,PPIs)
s=''.join(s)
write_file('output/TAD.EnrichPPIs.txt',s)



