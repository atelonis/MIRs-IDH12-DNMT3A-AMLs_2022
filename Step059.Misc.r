
library(data.table)
library(VennDiagram)

hypergeom.test=function(k,s,M,N){
	# k: number of successes
	# s: sample size
	# M: Number of successes in the population
	# N: population size
	myFE = (k/s) / (M/N)
	if (is.na(myFE)){
		p=NA
	} else if (myFE>1){
		p=1-phyper( k-1, M, N-M, s )
	} else {
		p=phyper( k, M, N-M, s )
	}
	return(c(myFE,p))
}


########################################################################################
###  PART 1 | Are genes associated with the same or different cytosines in the TAD?  ###
########################################################################################

m0=read.table('TCGA/output/Correlation/Step3/TADsets.GenesCsInTAD.txt',header=TRUE,stringsAsFactors=FALSE)
m0=as.matrix(m0)
m=m0[which(as.numeric(m0[,5])>10 & as.numeric(m0[,6])>10 & as.numeric(m0[,7])>0),]
# Required at least one common mC - The ones that have no common mC, that depletion is not significant

mh=NULL
for (i in c(1:nrow(m))){
	cHT=hypergeom.test(as.numeric(m[i,7]),as.numeric(m[i,6]),as.numeric(m[i,5]),as.numeric(m[i,2]))
	mh=rbind(mh,c(m[i,],cHT[1],cHT[2]))
}
mh=cbind(mh,p.adjust(mh[,9],method='fdr'))
colnames(mh)=c(colnames(m),'FE','PValue','FDR')


uTads=sort(unique(mh[,1]))

vv=NULL
sv=NULL
for (i in uTads){
	wRows=which(mh[,1]==i)
	vv=c(vv,as.numeric(mh[wRows,'FE']))
	sv=c(sv,as.numeric(mh[wRows,'FDR']))
}
sv[sv==0]=1e-17
vCols=rep('#AAAAAA',length(vv))
vCols[which(sv<0.05)]='#B22222AA'

pdf('TCGA/output/TAD.Enrich_GeneC_pairs.pdf',height=4,width=4) ### FIGURE 1E
par(mar=c(5,5,2,2))
plot(log2(vv),-log10(sv),xlab='Fold enrichment in common mCs (log2)',
	 ylab='FDR (-log10)',col=vCols,pch=19)
dev.off()



################################################
###  PART 2 | Plot cumulative distributions  ###
################################################


#=========================
# 1  |  Gene Parameters
#=========================

dir.create('TCGA/output/CumDist',showWarnings=FALSE)

tadBack   = read.table('TCGA/output/DAVID/Background.TADs_WholeNames.txt',stringsAsFactors=FALSE,sep='\t')[,1]

fileNames=c('TCGA/output/Correlation/Step3/TADsets.Negat_Overlap.txt',
			'TCGA/output/Correlation/Step3/TADsets.Negat_Intermediate.txt',
			'TCGA/output/Correlation/Step3/TADsets.Negat_Long.txt',
			'TCGA/output/Correlation/Step3/TADsets.Posit_Overlap.txt',
			'TCGA/output/Correlation/Step3/TADsets.Posit_Intermediate.txt',
			'TCGA/output/Correlation/Step3/TADsets.Posit_Long.txt',
			'TCGA/output/Correlation/Step3/TADsets_W.Overlap.txt')
names(fileNames)=c('Posit_Proximal','Posit_Intermediate','Posit_Long',
				   'Negat_Proximal','Negat_Intermediate','Negat_Long',
				   'W set')

geneSets = NULL
for (i in names(fileNames)){
	geneSets[[i]]=read.table(fileNames[i],stringsAsFactors=FALSE,sep='\t')[,1]
}
geneCols=c( '#F8C1B8','#FF6655','#A6093D',  '#D2C2E6','#9063CD','#440099', 'forestgreen')
names(geneCols)=names(fileNames) #c('Overlap','Intermediate','Long','W set')

myScores=read.table('HumnaGenome/GeneParameters.txt',header=TRUE)
rownames(myScores)=myScores[,1]
myScores=myScores[,-1]

myREG=c('Length','GC','EvolCons')

xLAB=c('Length (log2 nt)','GC Content','Evolutionary conservation')
names(xLAB)=myREG

yLIM=c(0.2,0.4,0.2)
names(yLIM)=myREG

myRB1=c(0,0,-1)
names(myRB1)=myREG

myRB2=c(22,1,1)
names(myRB2)=myREG


for (reg in myREG){
	
	rb=seq(myRB1[reg],myRB2[reg],length.out=100)
	
	eiName = intersect( grep('Intron',colnames(myScores)), grep(reg,colnames(myScores)))
	eiBACK = myScores[tadBack ,eiName]
	b      = ecdf(eiBACK)
	bb     = b(rb)
	
	pdf(paste('TCGA/output/CumDist/TAD',reg,'pdf',sep='.'),width=4,height=4)  ### FIGURE 2A 3B
	par(mar=c(5,5,2,2))
	cY=-yLIM[reg]
	plot(rb,(-2*cY*bb)+cY,main='',type='l',lwd=2,pch=NA,col='gray',lty=2,
		 xlab=xLAB[reg],ylab='Cumulative fraction difference\nfrom background distribution',
		 xlim=c(myRB1[reg],myRB2[reg]),ylim=c(-yLIM[reg],yLIM[reg]))
	axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray',col.axis='gray')
	mtext('Cumulative background distribution',4,col='gray',line=3)
	for (i in names(geneSets)){
		cSCORE = myScores[geneSets[[i]],eiName]
		cECDF  = ecdf(cSCORE)
		cb     = cECDF(rb)
		L      = bb-cb
		lines(rb,L,col=geneCols[i],lwd=5)
	}
	lines(c(myRB1[reg],myRB2[reg]),c(0,0),col='black',lwd=2)
	dev.off()
}


#==================================
#  2  |  Repeat element densities
#==================================

myMirAlu=NULL
amColnames=NULL
for (ei in c('Intron')){
	for (sa in c('Sense','Antisense')){
		cTable=read.table(paste('HumanGenome/Density_Family',ei,sa,'txt',sep='.'),
						  sep='\t',stringsAsFactors=FALSE,check.names=FALSE)
		if (is.null(myMirAlu)){
			myMirAlu=cTable[,c('SINE/Alu','SINE/MIR','LINE/L1','LTR/ERVL')]
		} else {
			myMirAlu=cbind(myMirAlu,cTable[,c('SINE/Alu','SINE/MIR','LINE/L1','LTR/ERVL')])
		}
		amColnames=c(amColnames, paste('Alu',ei,sa,sep='.'),
								 paste('MIR',ei,sa,sep='.'),
								 paste('L1',ei,sa,sep='.'),
								 paste('ERVL',ei,sa,sep='.'))
	}
}
colnames(myMirAlu)=amColnames

xLIM=c(0.5,0.15,0.3,0.1)
names(xLIM)=c('Alu','MIR','L1','ERVL')

for (am in c('Alu','MIR')){
for (sa in c('Sense','Antisense')){
	
	amsa=paste(am,'Intron',sa,sep='.')
	rb=seq(0,xLIM[am],length.out=100)
	
	ei=amsa
	eiBACK = myMirAlu[tadBack,amsa]
	eiBACK = eiBACK[eiBACK>0]
	b      = ecdf(eiBACK)
	bb     = b(rb)
	
	pdf(paste('TCGA/output/CumDist/TAD',amsa,'pdf',sep='.'),width=4,height=4) ### FIGURE 2C 3D
	par(mar=c(5,5,2,2))
	cY=-0.2
	plot(rb,(-2*cY*bb)+cY,main='',type='l',lwd=2,pch=NA,col='gray',lty=2,
		 xlab=paste(am,'density'),ylab='Cumulative fraction difference\nfrom background distribution',
		 xlim=c(0,xLIM[am]),ylim=c(cY,-cY),xaxt='n')
	axis(1,seq(0,xLIM[am],length.out=6),100*seq(0,xLIM[am],length.out=6))
	axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray',col.axis='gray')
	mtext('Cumulative background distribution',4,col='gray',line=3)
	for (i in names(geneSets)){
		eiSCORE = myMirAlu[geneSets[[i]],amsa]
		eiSCORE = eiSCORE[eiSCORE>0]
		cECDF   = ecdf(eiSCORE)
		cb      = cECDF(rb)
		L       = bb-cb
		lines(rb,L,col=geneCols[i],lwd=5)
	}
	lines(c(0,xLIM[am]),c(0,0),col='black',lwd=2)
	dev.off()
}
}


