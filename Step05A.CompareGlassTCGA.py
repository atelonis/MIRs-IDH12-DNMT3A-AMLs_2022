
import os
import subprocess

from Functions import *

create_dir('ForComparingGlassTCGA')

d={}
with open('tadCytosines.fa') as f:
	name=f.readline()
	seq=f.readline()
	while name:
		
		name=name.replace('\n','').replace('>','').split('|')
		analysisName='%s.%s.NA0' % (name[0],name[2])
		seq=seq.replace('\n','')
		
		if seq[1]=='C':
			d[analysisName]='\t'.join([name[0],str(int(name[1])+1),name[2]])+'\n'
		else:
			d[analysisName]='\t'.join([name[0],name[1],str(int(name[2])-1)])+'\n'
		name=f.readline()
		seq=f.readline()

mCs=read_file('output/Correlation/Step3/Cs.Background.txt')
mCs=list(set([d[x] for x in mCs]))

write_file('ForComparingGlassTCGA/Glass.bed',''.join(mCs))
	
bash1="bedtools intersect -a ForComparingGlassTCGA/Glass.bed -b cg.bed -wa -wb > ForComparingGlassTCGA/BI.Glass.TCGA.bed"
subprocess.call(['bash','-c',bash1])

g=len(mCs)
t=len(read_file('TCGA/cg.bed'))
bi=len(read_file('ForComparingGlassTCGA/BI.Glass.TCGA.bed'))

s='Number for Venn Diagram\n'
s+='Total in Glass\t%s\n' % (g)
s+='Total in TCGA\t%s\n' % (t)
s+='Intersection\t%s\n' % (bi)
s+='Unique in Glass\t%s\n' % (g-bi)
s+='Unique in TCGA\t%s\n' % (t-bi)

write_file('ForComparingGlassTCGA/ForVenn.txt',s)


