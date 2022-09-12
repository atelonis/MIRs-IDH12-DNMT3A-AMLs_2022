
import os
from Functions import *

create_dir('output/GEAR')

for am in ['Alu','MIR','L1']:
	allS=read_file('output/SAM/Methylation.%s.IDH_vs_DNMT.SAM.IncludedInSAM.txt' % (am))
	allS=['%s.%s-%s' % (x.split('.')[0],x.split('.')[1],x.split('.')[1]) for x in allS]
	sigS=read_file('output/SAM/Methylation.%s.IDH_vs_DNMT.SAM.PositSign.txt' % (am))[1:]
	sigS=[x.split('\t')[2].replace('"','') for x in sigS]
	sigS=['%s.%s-%s' % (x.split('.')[0],x.split('.')[1],x.split('.')[1]) for x in sigS]
	
	write_file('output/GEAR/SAM_%s.Background.txt' % (am),'\n'.join(allS)+'\n')
	write_file('output/GEAR/SAM_%s.DiffMeth.txt' % (am),'\n'.join(sigS)+'\n')
	
	myFiles={'%s' % (am):'output/GEAR/SAM_%s.DiffMeth.txt' % (am)}
	GEAR('Human',myFiles,'output/GEAR/SAM_%s.Background.txt' % (am),'SAM_%s' % (am),False)



