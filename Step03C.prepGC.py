
import os

def read_file(myfile):
	f=open(myfile,'r')
	fl=[x.replace('\n','') for x in f.readlines()]
	f.close()
	return fl

def write_file(MyFile,MyString):
	f=open(MyFile,'w')
	f.write(MyString)
	f.close()
	return

def create_dir(myDIR):
	if not os.path.isdir(myDIR):
		os.mkdir(myDIR)
	return

def get_from_dict(myD,myKey,myKE):
	try:
		return myD[myKey]
	except KeyError:
		return myKE

create_dir('output/GC_Ccount')

backGenes=[]
with open('output/Correlation/Step1/TadCorrelations.txt') as f:
	i=f.readline()
	while i:
		i=i.replace('\n','').split('\t')
		if i[2]=='nan' or 'GeneBody' not in i[4]:
			i=f.readline()
			continue
		backGenes+=[i[1].split('|')[0]]
		i=f.readline()

cCount={} # Remember, these are in CpG context
for i in list(set(backGenes)):
	cCount[i]=backGenes.count(i)

gc={}
for i in read_file('HumanGenome/GeneParameters.txt'):
	i=i.split('\t')
	if i[0]=='GeneName':
		gcExon   = i.index('GC_Exon')
		gcIntron = i.index('GC_Intron')
		eLength  = i.index('ExonLength')
		iLength  = i.index('IntronLength')
	else:
		gc[i[0].split('|')[1]]=[i[x] for x in [gcExon,gcIntron,eLength,iLength]]

aluD={}
mirD={}
for ei in ['Exon','Intron']:
	for i in read_file('HumanGenome/Number_Family.%s.Sense.txt' % (ei)):
		i=i.split('\t')
		if i[0]=='GeneBases':
			alu=i.index('SINE/Alu')
			mir=i.index('SINE/MIR')
		else:
			if ei=='Exon':
				aluD[i[0].split('|')[1]] = float(i[alu+1])
				mirD[i[0].split('|')[1]] = float(i[mir+1])
			else:
				aluD[i[0].split('|')[1]] += float(i[alu+1])
				mirD[i[0].split('|')[1]] += float(i[mir+1])

s='C_Count\tGC_Exons\tGC_Introns\tExon_Length\tIntron_Length\tSINE/Alu\tSINE/MIR\n'
for i in cCount.keys():
	try:
		alu=str(aluD[i]/((2**float(gc[i][2]))+(2**float(gc[i][3]))))
		mir=str(mirD[i]/((2**float(gc[i][2]))+(2**float(gc[i][3]))))
		s+='\t'.join([i,str(cCount[i])]+gc[i]+[alu,mir])+'\n'
	except KeyError:
		continue

write_file('output/GC_Ccount/TAD.GeneTable.txt',s)




