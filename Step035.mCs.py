
import os
from Functions import *

repD={}
for i in read_file('bed/BedIntersect.TadCs.Repeats.bed'):
	i=i.split('\t')
	cC=i[3]
	cPos=int(i[3].split('.')[1])
	if int(i[7])<=cPos and cPos<=int(i[8]):
		try:
			#repD[cC]+=[i[9].split('|')[0]]
			repD[i[9].split('|')[0]]+=[cC]
		except KeyError:
			repD[i[9].split('|')[0]]=[cC]
			#repD[cC]=[i[9].split('|')[0]]

corD={}
for i in read_file('output/Correlation/Step2/TadCorrelations.txt')[1:]:
	i=i.split('\t')
	cC=i[0].replace('.NA0','')
	try:
		corD[i[7]]+=[cC]
	except KeyError:
		corD[i[7]]=[cC]
cKeys=sorted(corD.keys())

allCs=[x.split('\t')[3] for x in read_file('bed/ForTadCEvol.bed')]
repCs=[item for sublist in repD.values() for item in sublist]
corCs=[item for sublist in corD.values() for item in sublist]

s='Background\t'+'\t'.join(cKeys)+'\n'

s+='NoRep\t%s' % (len(list(set(allCs).difference(set(repCs)))))
for j in cKeys:
	s+='\t%s' % (len(list(set(corD[j]).difference(set(repCs)))))
s+='\n'

s+='InReps\t%s' % (len(list(set(allCs).intersection(set(repCs)))))
for j in cKeys:
	s+='\t%s' % (len(list(set(corD[j]).intersection(set(repCs)))))
s+='\n'

for i in sorted(repD.keys()):
	s+='%s\t%s' % (i,len(repD[i]))
	for j in cKeys:
		s+='\t%s' % (len(list(set(repD[i]).intersection(set(corD[j])))))
	s+='\n'
write_file('output/TAD.C_in_Repeats.txt',s)





