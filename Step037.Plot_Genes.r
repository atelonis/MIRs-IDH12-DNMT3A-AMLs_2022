
#=========================
# 1  |  Gene Parameters
#=========================

dir.create('output/CumDist',showWarnings=FALSE)

tadBack=read.table('output/Correlation/Step3/TADsets.Background.ENSG.txt',stringsAsFactors=FALSE,sep='\t')[,1]

fileNames=c('output/Correlation/Step3/TADsets.Negat_Overlap.ENSG.txt',
			'output/Correlation/Step3/TADsets.Negat_Intermediate.ENSG.txt',
			'output/Correlation/Step3/TADsets.Negat_Long.ENSG.txt',
			'output/Correlation/Step3/TADsets_W.Overlap.ENSG.txt')
names(fileNames)=c('Overlap','Intermediate','Long','W set')

geneSets = NULL
for (i in names(fileNames)){
	geneSets[[i]]=read.table(fileNames[i],stringsAsFactors=FALSE,sep='\t')[,1]
}
geneCols=c('#FF8674','#FF4635','#A6093D','forestgreen')
geneCols=c('#F8C1B8','#FF6655','#A6093D','forestgreen')
names(geneCols)=c('Overlap','Intermediate','Long','W set')

myScores=read.table('HumanGenome/GeneParameters.txt',header=TRUE)
rownames(myScores)=myScores[,1]
myScores=myScores[,-1]

myREG=c('GC','EvolCons')

xLAB=c('GC Content','Evolutionary conservation')
names(xLAB)=myREG

yLIM=c(0.4,0.2)
names(yLIM)=myREG

myRB1=c(0,-1.5)
names(myRB1)=myREG

myRB2=c(1,1.5)
names(myRB2)=myREG


for (reg in myREG){
	
	rb=seq(myRB1[reg],myRB2[reg],length.out=100)
	
	eiName = intersect( grep('Intron',colnames(myScores)), grep(reg,colnames(myScores)))
	eiBACK = myScores[tadBack ,eiName]
	b      = ecdf(eiBACK)
	bb     = b(rb)
	
	
	pdf(paste('output/CumDist/TAD',reg,'pdf',sep='.'),width=4,height=4) ### FIGURE 2A 3B
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
	legend('bottomright',names(geneCols),col=geneCols,pch=NA,lwd=3,cex=0.5)
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

aa=NULL
for (am in c('Alu','MIR')){
for (sa in c('Sense','Antisense')){
	
	amsa=paste(am,'Intron',sa,sep='.')
	rb=seq(0,xLIM[am],length.out=100)
	
	ei=amsa
	eiBACK = myMirAlu[tadBack,amsa]
	eiBACK = eiBACK[eiBACK>0]
	b      = ecdf(eiBACK)
	bb     = b(rb)
	
	pdf(paste('output/CumDist/TAD',amsa,'pdf',sep='.'),width=4,height=4)  ### FIGURE 2C 3D
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
	
		cKS=ks.test(eiSCORE,ecdf(eiBACK))
		aa=rbind(aa,c(i,ei,cKS$p.value))
	}
	lines(c(0,xLIM[am]),c(0,0),col='black',lwd=2)
	dev.off()
}
}




