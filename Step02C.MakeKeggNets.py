
import os
from Functions import *

create_dir('KEGG')
create_dir('output/ForCytoscape')

lipidPaths=['hsa00565:EtherLipidMetabolism',
			'hsa00561:GlycerolipidMetabolism',
			'hsa00564:GlycerophospholipidMetabolism',
			'hsa00600:SphingolipidMetabolism',
			'hsa00601:GlycosphingolipidBiosynthesisLactoAndNeolactoSeries',
			'hsa00603:GlycosphingolipidBiosynthesisGloboAndIsogloboSeries',
			'hsa00604:GlycosphingolipidBiosynthesisGanglioSeries',
			'hsa00061:FattyAcidBiosynthesis',
			'hsa00062:FattyAcidElongation',
			'hsa00071:FattyAcidDegradation',
			'hsa01040:BiosynthesisOfUnsaturatedFattyAcids',
			'hsa01212:FattyAcidMetabolism']

entrez2name={}
for i in read_file('HumanGenome/mart_export.Name_IDs.txt')[1:]:
	i=i.split('\t')
	if i[3]=='':
		continue
	entrez2name[i[3]]=i[2]


### Convert KGML to Network

for cPath in lipidPaths:
	print cPath
	
	cXML='HumanGenome/KGML/%s.kgml' % (cPath.split(':')[0])
	data=read_file(cXML)
	
	substrNet={} # Dictionary with the reactions that each compound is as a substrate
	prodNet={}  # Dictionary with the reactions that each compound is a product
	react2gene={} # which genes catalyze one reaction
	substrNet.clear()
	prodNet.clear()
	react2gene.clear()
	isReaction=False
	for i in data:
		while True: # remove starting spaces
			if i[0]==' ':
				i=i[1:]
			else:
				break
		if i[:6]=='<entry': # this is a new entry
			i=i.replace('" ','"_')
			i=i.split('_')
			cID=i[0].split('"')[1]
			if 'type="gene"' in i:
				cGenes=i[1].split('"')[1].replace('hsa:','').split(' ')
				react2gene[cID]=cGenes
		elif i[:9]=='<reaction': # this is a new reaction
			cReactID=i.split(' ')[1].split('"')[1]
			isReaction=True
			if 'type="reversible"' in i:
				isRev=True
			else:
				isRev=False
		elif isReaction:
			if i=='</reaction>':
				isReaction=False
				continue
			i=i.split(' ')
			addSusbtr=False
			addProd=False
			if i[0] in ['<substrate','<product']:
				cComp=i[1].split('"')[1]
				if isRev:
					addSusbtr=True
					addProd=True
				if i[0]=='<substrate':
					addSubstr=True
				elif i[0]=='<product':
					addProd=True
				
				if addSubstr:
					try:
						substrNet[cComp]+=[cReactID]
					except KeyError:
						substrNet[cComp]=[cReactID]
				
				if addProd:
					try:
						prodNet[cComp]+=[cReactID]
					except KeyError:
						prodNet[cComp]=[cReactID]
	
	s=[]
	for i in substrNet.keys():
		for j in substrNet[i]:
			cGenes=react2gene[j]
			for k in cGenes:
				s+=['%s\tuses\t%s\n' % (k,i)]
	for i in prodNet.keys():
		for j in prodNet[i]:
			cGenes=react2gene[j]
			for k in cGenes:
				s+=['%s\tproduces\t%s\n' % (k,i)]
	
	s='Source\tInteraction\tTarget\n' + ''.join(list(set(s)))
	write_file('KEGG/%s.EntrezMetaboliteNetwork.txt' % (cPath),s)
	
	d={}
	d.clear()
	for i in read_file('KEGG/%s.EntrezMetaboliteNetwork.txt' % (cPath))[1:]:
		i=i.split('\t')
		try:
			d[i[0]]+=[i[2]]
		except KeyError:
			d[i[0]]=[i[2]]
	for i in d.keys():
		d[i]=sorted(list(set(d[i])))
	
	dKeys=d.keys()
	for i in range(len(dKeys)-1):
		iK=dKeys[i]
		for j in range(i+1,len(dKeys)):
			jK=dKeys[j]
			if d[iK]==d[jK]:
				continue
			elif len([x for x in d[iK] if x in d[jK]])>=1:
				iName=entrez2name[iK]
				jName=entrez2name[jK]
				ij=sorted([iName,jName])
				s+=['%s\tmetab\t%s\n' % (ij[0],ij[1])]
	s='Source\tInteraction\tTarget\n' + ''.join(list(set(s)))
	write_file('KEGG/%s.GeneNetwork.txt' % (cPath),s)
	
	Genes=[x.split('\t')[0] for x in s.split('\n')[1:-1]] + [x.split('\t')[2] for x in s.split('\n')[1:-1]]
	write_file('KEGG/%s.GeneNames.txt' % (cPath),'\n'.join(Genes)+'\n')


### Add information from correlations

tadCorrs=read_file('output/Correlation/Step2/TadCorrelations.txt')
corrD={}
for i in tadCorrs[1:]:
	i=i.split('\t')
	try:
		corrD[i[1].split('|')[0]]+=[i[7]]
	except KeyError:
		corrD[i[1].split('|')[0]]=[i[7]]

s='Source\tInteraction\tTarget\tSourceAttribute\tTargetAttribute\n'
for cPath in lipidPaths:
	cNet=read_file('KEGG/%s.GeneNetwork.txt' % (cPath))[1:]
	for i in cNet:
		i=sorted([i.split('\t')[0],i.split('\t')[2]])
		c=[','.join(list(set(get_from_dict(corrD,x,['notInCorrs'])))) for x in i]
		newS='%s\tmetab\t%s\t%s\t%s\n' % (i[0],i[1],c[0],c[1])
		if newS in s:
			continue
		s+=newS
write_file('output/ForCytoscape/LipidNetwork.txt',s) ### FIGURE 1H --> beautify on cytoscape



