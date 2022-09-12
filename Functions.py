
import os
import sys
import math
import numpy
import random
import subprocess

from multiprocessing import Pool

from __future__ import division

import scipy
import scipy.stats
from scipy.stats import hypergeom
from scipy.stats import chi2_contingency
from scipy.stats.contingency import margins

from statsmodels.distributions.empirical_distribution import ECDF

###############
###  BASIC  ###
###############

def read_file(myfile):
	f=open(myfile,'r')
	fl=[x.replace('\n','').replace('\r','') for x in f.readlines()]
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
		return myD[myK]
	except KeyError:
		return myKE

def save_rep_matrix(myMatrix,myRows,myCols,myName):
	numpy.savetxt('temp.%s.txt' % (myName),myMatrix,fmt='%s',delimiter='\t')
	cM=read_file('temp.%s.txt' % (myName))
	s=myCols+'\n'
	for i in range(len(myRows)):
		s+='%s\t%s\n' % (myRows[i],cM[i])
	write_file(myName,s)
	os.remove('temp.%s.txt' % (myName))
	return

def hypergeom_test(k,s,M,N):
	# k: number of successes
	# s: sample size
	# M: Number of successes in the population
	# N: population size
	
	try:
		myFE = (float(k)/float(s)) / (float(M)/float(N))
	except ZeroDivisionError:
		return ['NA','NA']
	
	if myFE>1:
		p=1-hypergeom.cdf( k-1, N, M, s)
	else:
		p=hypergeom.cdf( k, N, M, s)
	
	return [ myFE, p ]


#############
###  ECO  ###
#############

def compute_ECO_KS(myArgs):
	
	myParam = myArgs[0]
	myDict  = myArgs[1]
	myGenes = myArgs[2]
	myBack  = myArgs[3]
	myOut   = myArgs[4]
	
	myBack=[x for x in myBack if x in myDict.keys()]
	valueBack=[myDict[x] for x in myBack]
	
	myGenes=[x for x in myGenes if x in myDict.keys()]
	valueGENES=[myDict[x] for x in myGenes]
	
	cKS=scipy.stats.kstest(valueGENES,ECDF(valueBack))
	
	if numpy.median(valueGENES) > numpy.median(valueBack):
		myDirection='Higher_Median'
	elif numpy.median(valueGENES) < numpy.median(valueBack):
		myDirection='Lowrt_Median'
	else:
		myDirection='Equal_Median'
	
	write_file('%s.KS_Test.%s.txt' % (myOut,myParam), '%s\t%s\t%s\n' % (myParam,str(cKS[1]),myDirection) )
	
	return 


def eco_analysis(ecoGENES,ecoBack,ecoSpecies,ecoOUT):
	# ECO stands for Exonic COntent
	# (the script was inititially designed for exons/intron ratio but the expanded to include
	#  parameters like GC content, evolutionary conservation, etc)
	
	# File required in this dir: GeneParameters.txt
	ecoDir='%s/GeneParameters.txt' % (ecoSpecies)
	ecoFile=read_file(ecoDir)
	
	ecoPAR={}
	cPars=ecoFile[0].split('\t')[1:]
	for i in cPars:
		ecoPAR[i]={}
		cIndex=ecoFile[0].split('\t').index(i)
		
		for j in ecoFile[1:]:
			j=j.split('\t')
			if j[cIndex]=='NA':
				continue
			if i in ['ExonLength','IntronLength']: # do calculations at linear and not log-scale
				ecoPAR[i][j[0]]=2.0**float(j[cIndex])
			else:
				ecoPAR[i][j[0]]=float(j[cIndex])
	
	ecoArgs=[ [x,ecoPAR[x],ecoGENES,ecoBack,ecoOUT] for x in cPars]
	
	p=Pool(15)
	p.map(compute_Z_KS,ecoArgs)
	

def get_ranges(eip):
	r={}
	r.clear()
	with open('HumanGenome/%sRanges.bed' % (eip)) as f:
		i=f.readline()
		while i:
			i=i.split('\t')
			try:
				r[i[3]]+=int(i[2])-int(i[1])
			except KeyError:
				r[i[3]]=int(i[2])-int(i[1])
			i=f.readline()
	return r


def get_counts(eip):
	r={}
	r.clear()
	with open('HumanGenome/%sRanges.bed' % (eip)) as f:
		i=f.readline()
		while i:
			i=i.split('\t')
			try:
				r[i[3]]+=1
			except KeyError:
				r[i[3]]=1
			i=f.readline()
	return r


def get_evol_cons(eip):
	r={}
	r.clear()
	with open('bed/BedInters.%s.EvolCons.bed' % (eip)) as f:
		i=f.readline()
		while i:
			i=i.replace('\n','').split('\t')
			try:
				r[i[3]]+=[float(i[10])]
			except KeyError:
				r[i[3]]=[float(i[10])]
			i=f.readline()
	for i in r.keys():
		r[i]=numpy.median(r[i])
	return r
		

def get_gc_content(eip):
	r={}
	r.clear()
	with open('HumanGenome/%s.fa' % (eip)) as f:
		cGene = f.readline()
		cSeq  = f.readline()
		while cGene:
			cGene='%s|%s' % (cGene.split('|')[3],cGene.split('|')[4])
			GC=cSeq.count('G')+cSeq.count('C')
			try:
				r[cGene]+=[ [GC,float(len(cSeq))] ]
			except KeyError:
				r[cGene]=[ [GC,float(len(cSeq))] ]
			cGene = f.readline()
			cSeq  = f.readline()
	for i in r.keys():
		r[i]=float(sum([x[0] for x in r[i]]))/float(sum([x[1] for x in r[i]]))
	return r


#############
###  RED  ###
#############

def red_one_comb(multi_ARGS):
	
	multi_FILE  = multi_ARGS[0] # '%sNumber_Family.%s.txt' % (repDir,i)
	multi_BACK  = multi_ARGS[1] # redBack
	multi_GENES = multi_ARGS[2] # redGENES
	multi_REPS  = multi_ARGS[3] # reps
	
	dens={}
	dens.clear()
	multi_FILE=read_file(multi_FILE)
	repOrder=multi_FILE[0].split('\t')
	for j in multi_FILE[1:]:
		j=j.split('\t')
		dens[j[0]]={}
		for k in range(2,len(j)):
			try:
				dens[j[0]][repOrder[k-1]]=float(j[k])/float(j[1])
			except ZeroDivisionError:
				dens[j[0]][repOrder[k-1]]=0.0
	
	backGENES    = [x for x in multi_BACK if x in dens.keys()]
	i_redGENES = [x for x in multi_GENES if x in dens.keys()]
	
	multi_COMB = multi_ARGS[0].split('/')[-1].split('.')
	r='%s.%s\n' % ( multi_COMB[1], multi_COMB[2] ) ### string to return ###
	for cREP in multi_REPS:
		r+='%s' % (cREP)
		
		cGENES=[x for x in i_redGENES if dens[x][cREP]!=0.0]
		
		#=== Hypergeometric Test
		# how many genes of the background have the specific repeat?
		hyperBackRep=len([x for x in backGENES if dens[x][cREP]!=0.0])
		
		try:
			foldEnrich=(float(len(cGENES))/float(len(i_redGENES))) / (float(hyperBackRep)/float(len(backGENES)))
		except ZeroDivisionError:
			r+='\t0\tx_x_x_x\t99.0\t99.0\t99.0\n'
			continue
		
		r+='\t%s' % (foldEnrich)
		
		lGenes=len(cGENES)-1
		if lGenes==-1:
			lGenes=0
		r+='\t%s_%s_%s_%s' % (lGenes, len(backGENES), hyperBackRep, len(i_redGENES))
		
		if foldEnrich>1:
			pVal=1-hypergeom.cdf( lGenes, len(backGENES), hyperBackRep, len(i_redGENES) )
		elif foldEnrich==0:
			pVal=99.0
		else:
			pVal=hypergeom.cdf( lGenes, len(backGENES), hyperBackRep, len(i_redGENES) )
		r+='\t%s' % (pVal)
		
		#=== Now on the genes that have repeats
		if cGENES==[]:
			r+='\t99.0\t99.0\t99.0\n'
			continue
		
		#=== Kolmogorov-Smirnov test
		densGENES=[dens[x][cREP] for x in cGENES]
		densBACK=[dens[x][cREP] for x in backGENES if dens[x][cREP]!=0]
		cKS=scipy.stats.kstest(densGENES,ECDF(densBACK))
		r+='\t%s' % (cKS[1])
		r+='\t%s' % (numpy.mean(densGENES)/numpy.mean(densBACK))
		
		r+='\n'
	return r


def red_analysis(redGENES,redBack,redSpecies,redOUT):
	# RED stands for Repeat Element Densities
	
	# === Let's start
	repDir=redSpecies
	reps=read_file('%s/AllRepeats.Family.txt' % (repDir))
	reps.sort()
	
	# === Combinations of where/how to look
	geneFeats=['Exon','Intron']
	orients=['Sense','Antisense']
	combs=[['%s.%s' % (x,y) for x in geneFeats] for y in orients]
	combs=combs[0]+combs[1]
	
	# === Prepare and execute parallel runnings
	forMulti = []
	for i in combs:
		forMulti += [[ '%sNumber_Family.%s.txt' % (repDir,i), redBack, redGENES, reps ]]
	
	print 'pooling...'
	p=Pool(69)
	R=p.map(red_one_comb,forMulti)
	
	# === Put the results in the arrays
	
	myKS=numpy.array([[99.0]*len(combs)]*len(reps)) # rows=reps, columns=combs
	myED=numpy.array([[99.0]*len(combs)]*len(reps)) # rows=reps, columns=combs
	
	hyperV=numpy.array([['xxxxxx_xxxxxx_xxxxxx_xxxxxx']*len(combs)]*len(reps)) # rows=reps, columns=combs
	hyperF=numpy.array([[0.0]*len(combs)]*len(reps)) # rows=reps, columns=combs
	hyperP=numpy.array([[99.0]*len(combs)]*len(reps)) # rows=reps, columns=combs
	
	for i in R:
		iCOMB = i.split('\n')[0]
		for j in i.split('\n')[1:-1]:
			j=j.split('\t')
			jRep=reps.index(j[0]) # DNA/Kolobok DNA/MULE-MuDR
			hyperF[ jRep, combs.index(iCOMB) ] = float(j[1]) # foldEnrich
			hyperV[ jRep, combs.index(iCOMB) ] = j[2] # '%s_%s_%s_%s' % (lGenes, len(backGENES), hyperBackRep, len(i_redGENES))
			hyperP[ jRep, combs.index(iCOMB) ] = float(j[3]) # pVal
			myKS[ jRep, combs.index(iCOMB) ] = float(j[4]) # cKS[1]
			myED[ jRep, combs.index(iCOMB) ] = float(j[5]) # numpy.mean(densGENES)/numpy.mean(densBACK)
	
	# === Convert to strings and write the files
	# Hypergeometric
	hV='\t'.join(combs)+'\n'
	hP='\t'.join(combs)+'\n'
	hF='\t'.join(combs)+'\n'
	
	# KS
	cKS='\t'.join(combs)+'\n'
	cED='\t'.join(combs)+'\n'
	for i in range(len(reps)):
		hV+='%s\t%s\n' % (reps[i],'\t'.join([str(x) for x in hyperV[i,:]]))
		hF+='%s\t%s\n' % (reps[i],'\t'.join([str(x) for x in hyperF[i,:]]))
		hP+='%s\t%s\n' % (reps[i],'\t'.join([str(x) for x in hyperP[i,:]]))
		
		cKS +='%s\t%s\n' % (reps[i],'\t'.join([str(x) for x in myKS[i,:]]))
		cED +='%s\t%s\n' % (reps[i],'\t'.join([str(x) for x in myED[i,:]]))
		
	
	write_file(redOUT+'.Hypergeometric.ParametersUsed.txt',hV)
	write_file(redOUT+'.Hypergeometric.FoldEnrichment.txt',hF)
	write_file(redOUT+'.Hypergeometric.PValues.txt',hP)
	
	write_file(redOUT+'.KStest.Dens.txt',cKS)
	write_file(redOUT+'.EnrichmentDensity.txt',cED)
	
	return


#===========================================================
#======== Genomic Elements Association with Regions ========
#======== *       *        *                *       ========
#===========================================================

#-----------------------------------------------------------
def residuals(observed, expected):
	return (observed - expected) / numpy.sqrt(expected)

def stdres(observed, expected):
	n = observed.sum()
	rsum, csum = margins(observed)
	rsum = rsum.astype(numpy.float64)
	csum = csum.astype(numpy.float64)
	v = csum * rsum * (n - rsum) * (n - csum) / n**3
	return (observed - expected) / numpy.sqrt(v)

#-----------------------------------------------------------

def fill_row(myA,myI,myC,myDict,myNames):
	for i in myC:
		try:
			myA[ myI, myNames.index(myDict[i]) ] += 1
		except KeyError:
			myA[ myI, len(myNames) ] += 1
	return myA

def run_CHIs(myA, myRow, myCol, myOut):
	
	# Compute percentages
	myAF=myA.astype(numpy.float64)
	myP=numpy.rint( (100.0 * myAF.transpose() / myAF.sum(axis=1)).transpose() )
	myP=myP.astype(numpy.int)
	
	# Make one table of strings
	myW=numpy.core.defchararray.add( myA.astype(numpy.character), numpy.array([[' (']*myA.shape[1]]*myA.shape[0]) )
	myW=numpy.core.defchararray.add( myW, myP.astype(numpy.character) )
	myW=numpy.core.defchararray.add( myW, numpy.array([['%)']*myA.shape[1]]*myA.shape[0]) )
	
	# Run Chi-Squared tests
	resSTAR=['']*myA.shape[1]
	PValues=[]
	for i in range(1,len(myRow)):
		cTEST=myA[[0,i],:]
		cTEST[0,:]-=cTEST[1,:]
		cCHI = chi2_contingency(cTEST)
		cRES = residuals(cTEST, cCHI[3])
		
		PValues+=[cCHI[1]]
		if cCHI[1]<1e-4: # P value <10^(-4) and absolute residual >2
			resSTAR=numpy.vstack((resSTAR,[' *' if abs(x)>2 else '' for x in cRES[1,:]]))
		else:
			resSTAR=numpy.vstack((resSTAR,['']*myA.shape[1]))
	
	# Append the 'stars' in the matrix
	myW=numpy.core.defchararray.add( myW, resSTAR )
	
	# Add the Pvalues and row names 
	PValues=['-']+[str(x) for x in PValues]
	myW=numpy.hstack((myW,numpy.vstack((myRow,PValues)).transpose()))
	myW=myW[:,[-2]+range(len(myCol))+[-1]]
	
	# Write file
	s='-\t%s\tPValue\n' % ('\t'.join(myCol))
	s+='\n'.join(['\t'.join(myW[x,]) for x in range(myW.shape[0])])+'\n'
	write_file(myOut,s)

	return

def GEAR(species,forAnalysis,forBackground=True,gearPRE='',gearREP=True):
	
	#=========================================
	# 0 Creating output directories and prefix
	#=========================================
	create_dir('output')
	create_dir('output/GEAR')
	
	if gearPRE!='' and gearPRE[-1]!='.':
		gearPRE+='.'
	
	#=====================
	# 1 Getting Background
	#=====================
	
	cBack=read_file(forBackground)
	write_file('output/GEAR/%sBackground.txt' % (gearPRE),'\n'.join(cBack)+'\n')
	
	#==========================
	# 2 Finding the right files
	#==========================
	
	repeats_bed = "%s/RepeatMasker.bed" % (species)
	genes_bed   = "%s/GeneRanges.bed" % (species)
	gene_types  = "%s/mart_export.ForGEAR.txt" % (species)
	
	#===================
	# 3 Running bedtools
	#===================
	
	bash_Repeats = "bedtools intersect -a <(cat output/GEAR/%sBackground.txt | tr '.' '\t' | tr '-' '\t' | awk '$2=$2-1' | tr ' ' '\t') -b %s -wa -wb > output/GEAR/%sBedInters.C.Repeats.bed" % (gearPRE,repeats_bed,gearPRE) 
	bash_Genes   = "bedtools intersect -a <(cat output/GEAR/%sBackground.txt | tr '.' '\t' | tr '-' '\t' | awk '$2=$2-1' | tr ' ' '\t') -b <(cat %s | awk '{print $1,$2-2000,$3+2000,$4}' | awk '$2>0' | tr ' ' '\t') -wa -wb > output/GEAR/%sBedInters.C.GeneExt.bed" % (gearPRE,genes_bed,gearPRE)
	
	if gearREP:
		subprocess.call(['bash','-c',bash_Repeats])
	subprocess.call(['bash','-c',bash_Genes])
	
	
	#===============================
	# 4 Finding the overlaps with...
	#===============================
	
	###... Repeats
	if gearREP:
		repD={}
		catRepD={}
		allReps=[]
		allCatReps=[]
		with open('output/GEAR/%sBedInters.C.Repeats.bed' % (gearPRE)) as f:
			i=f.readline()
			while i:
				
				i=i.split('\t')
				
				cRep=i[6].split('|')[0]
				if cRep in ['Other','Unknown','Simple_repeat','Low_complexity'] or '?' in cRep:
					i=f.readline()
					continue
				
				repD['%s.%s-%s' % (i[0],int(i[1])+1,i[2])]    = i[6].split('|')[0]
				catRepD['%s.%s-%s' % (i[0],int(i[1])+1,i[2])] = i[6].split('|')[0].split('/')[0]
				allReps   +=[i[6].split("|")[0]]
				allCatReps+=[i[6].split("|")[0].split('/')[0]]
				
				i=f.readline()
		allReps=sorted(list(set(allReps)))
		allCatReps=sorted(list(set(allCatReps)))
	
	###... Genes
	geneTypes={}
	with open(gene_types) as f:
		i=f.readline()
		i=f.readline()
		while i:
			i=i.replace('\n','').split('\t')
			if i[2]=='protein_coding':
				geneTypes['%s|%s' % (i[0],i[1])] = 'proteincoding'
			else:
				geneTypes['%s|%s' % (i[0],i[1])] = 'noncoding'
			i=f.readline()
	
	geneD={}
	with open('output/GEAR/%sBedInters.C.GeneExt.bed' % (gearPRE)) as f:
		i=f.readline()
		while i:
			i=i.replace('\n','').split('\t')
			geneD['%s.%s-%s' % (i[0],int(i[1])+1,i[2])] = geneTypes[i[6]]
			i=f.readline()
	
	
	#==============
	# 5 Make tables
	#==============
	
	faKeys=sorted(forAnalysis.keys())
	
	G = numpy.array([[0]*3]*(len(faKeys)+1))
	if gearREP:
		R = numpy.array([[0]*(len(allReps)+1)]*(len(faKeys)+1))
		C = numpy.array([[0]*(len(allCatReps)+1)]*(len(faKeys)+1))
	
	# Background
	G = fill_row(G,0,cBack,geneD,['proteincoding','noncoding'])
	if gearREP:
		R = fill_row(R,0,cBack,repD,allReps)
		C = fill_row(R,0,cBack,catRepD,allCatReps)
	
	# Files
	for i in faKeys:
		cC=read_file(forAnalysis[i])
		G = fill_row( G, faKeys.index(i)+1, cC, geneD, ['proteincoding','noncoding'])
		if gearREP:
			R = fill_row( R, faKeys.index(i)+1, cC, repD,  allReps)
			C = fill_row( R, faKeys.index(i)+1, cC, catRepD,  allCatReps)
	
	# Exclude repeats with <50 represantation
	if gearREP:
		repIncl=[x for x,y in enumerate(numpy.max(R,0)) if y>=50]
		R=R[:,repIncl]
		allReps=[y for x,y in enumerate(allReps) if x in repIncl]
		
		catRepIncl=[x for x,y in enumerate(numpy.max(C,0)) if y>=50]
		C=C[:,catRepIncl]
		allCatReps=[y for x,y in enumerate(allCatReps) if x in catRepIncl]
	
	
	#======================
	# 6 Compute Chi-squared
	#======================
	run_CHIs(G, ['Background']+faKeys, ['ProteinCoding','Noncoding','Intergenic'], 'output/GEAR/%sGenes.txt' % (gearPRE))
	if gearREP:
		run_CHIs(R, ['Background']+faKeys, allReps+['NonRepeat'], 'output/GEAR/%sRepeats.txt' % (gearPRE))
		run_CHIs(C, ['Background']+faKeys, allCatReps+['NonRepeat'], 'output/GEAR/%sCategoryRepeats.txt' % (gearPRE))
	
	#=======================
	# 7 Hypergeometric tests
	#=======================
	
	h='File\tRepeat\tPopSize\tPopSucc\tSampSize\tSampSucc\tFoldEnrichment\tPValue\n'
	for i in faKeys:
		cC=read_file(forAnalysis[i])
		for j in allReps:
			
			popsiz = len(cBack)
			popsuc = R[ 0, allReps.index(j) ]
			samsiz = len(cC)
			samsuc = R[ faKeys.index(i)+1, allReps.index(j) ]
			
			foldEnrich= (samsuc/samsiz) / (popsuc/popsiz)
			
			if foldEnrich>1:
				pVal=1-hypergeom.cdf( samsuc-1, popsiz, popsuc, samsiz )
			else:
				pVal=hypergeom.cdf( samsuc, popsiz, popsuc, samsiz )
			
			h+='\t'.join([str(x) for x in [i, j, popsiz, popsuc, samsiz, samsuc, foldEnrich, pVal]])+'\n'
	
	write_file('output/GEAR/%sHypergeometrics.txt' % (gearPRE),h)
			
	return




