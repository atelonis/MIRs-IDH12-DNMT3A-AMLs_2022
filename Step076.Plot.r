
filter_zero_na=function(myList){
	if (length(which(is.na(myList)))>0){
		myList=myList[-which(is.na(myList))]
	}
	if (length(which(myList==0))>0){
		myList=myList[-which(myList==0)]
	}
	return(myList)
}

dir.create('GSE60055/output/CumDist',showWarnings=FALSE)

myCombs=NULL
for (i in list.files('GSE60055/output/SAM')){
	i=strsplit(i,'[.]')[[1]]
	if (i[1]!='MEP'){
		next
	}
	myCombs=c(myCombs,paste(i[1],i[2],sep='.'))
}
myCombs=sort(unique(myCombs))

myScores=read.table('MouseGenome/GeneParameters.txt',header=TRUE)
rownames(myScores)=myScores[,1]
myScores=myScores[,-1]

for (it in myCombs){
	print(it)
	
	background=read.table(paste('GSE60055/output/SAM/',it,'.SAM.IncludedInSAM.txt',sep=''),stringsAsFactors=FALSE)[,1]
	G=NULL
	G[['Upregulated']]=read.table(paste('GSE60055/output/SAM/',it,'.SAM.PositSign.txt',sep=''),stringsAsFactors=FALSE,header=TRUE)[,3]
	G[['Downregulated']]=read.table(paste('GSE60055/output/SAM/',it,'.SAM.NegatSign.txt',sep=''),stringsAsFactors=FALSE,header=TRUE)[,3]
	
	GCols=rbind( c('gray','red'), c('gray','blue') ) # no gray should be plotted from this script, we are only looking at introns
	colnames(GCols)=c('Exon','Intron')
	rownames(GCols)=c('Upregulated','Downregulated')
	
	#==================================
	#  1  |  Gene characteristics
	#==================================
	
	myREG=c('GC','EvolCons')
	
	xLAB=c('GC Content','Evolutionary conservation')
	names(xLAB)=myREG
	
	yLIM=c(0.2,0.2)
	names(yLIM)=myREG
	
	myRB1=c(0,-3)
	names(myRB1)=myREG
	
	myRB2=c(1,5)
	names(myRB2)=myREG
	
	for (reg in myREG){
		
		rb=seq(myRB1[reg],myRB2[reg],length.out=100)
		
		L=NULL
		for (pn in c('Upregulated','Downregulated')){
			for (ei in c('Intron')){
				
				eiName  = intersect( grep(ei,colnames(myScores)), grep(reg,colnames(myScores)))
				
				eiSCORE = filter_zero_na(myScores[G[[pn]],eiName])
				eiBACK  = filter_zero_na(myScores[background ,eiName])
				
				b=ecdf(eiBACK)
				bb=b(rb)
				
				cECDF=ecdf(eiSCORE)
				cb=cECDF(rb)
				
				L[[paste(pn,ei)]]=bb-cb
			}
		}
		
		pdf(paste('GSE60055/output/CumDist/Introns.',it,'.',reg,'.pdf',sep=''),width=4,height=4) ### FIGURE 5A 5B
		par(mar=c(4,4,3,3))
		
		cY=-yLIM[reg]
		plot(rb,(-2*cY*bb)+cY,main='',type='l',lwd=2,pch=NA,col='gray60',lty=2,mgp=c(1.6,0.6,0),
			 xlab=xLAB[reg],ylab='Cumulative fraction difference\nfrom background distribution',
			 xlim=c(myRB1[reg],myRB2[reg]),ylim=c(-yLIM[reg],yLIM[reg]))
		axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray60',col.axis='gray60',mgp=c(1.6,0.6,0))
		mtext('Cumulative background distribution',4,col='gray60',line=1.6)
		
		for (pn in c('Upregulated','Downregulated')){
			for (ei in c('Intron')){
				lines(rb,L[[paste(pn,ei)]],col=GCols[pn,ei],lwd=5)
			}
		}
		lines(c(myRB1[reg],myRB2[reg]),c(0,0),col='black',lwd=2)
		dev.off()
	
	}
	
	
	#==================================
	#  2  |  Repeat element densities
	#==================================
	
	myMirAlu=NULL
	amColnames=NULL
	for (ei in c('Exon','Intron')){
		for (sa in c('Sense')){ #,'Antisense')){
			cTable=read.table(paste('MouseGenome/Density_Family',ei,sa,'txt',sep='.'),
							  sep='\t',stringsAsFactors=FALSE,check.names=FALSE)
			if (is.null(myMirAlu)){
				myMirAlu=cTable[,c('SINE/Alu','SINE/MIR')]
			} else {
				myMirAlu=cbind(myMirAlu,cTable[,c('SINE/Alu','SINE/MIR')])
			}
			amColnames=c(amColnames,c(paste('Alu',ei,sep='.'),paste('MIR',ei,sep='.')))
		}
	}
	colnames(myMirAlu)=amColnames
	
	xLIM=c(0.25,0.1)
	names(xLIM)=c('Alu','MIR')
	
	yLIM=c(0.2,0.2)
	names(yLIM)=c('Alu','MIR')
	
	
	for (am in c('Alu','MIR')){
		
		rb=seq(0,xLIM[am],length.out=100)
		
		L=NULL
		for (pn in c('Upregulated','Downregulated')){
			for (ei in c('Exon','Intron')){
				
				eiName  = intersect( grep(ei,colnames(myMirAlu)), grep(am,colnames(myMirAlu)))
				
				eiSCORE = filter_zero_na(myMirAlu[G[[pn]],eiName])
				eiBACK  = filter_zero_na(myMirAlu[background ,eiName])
				
				b=ecdf(eiBACK)
				bb=b(rb)
				
				cECDF=ecdf(eiSCORE)
				cb=cECDF(rb)
				
				L[[paste(ei,pn)]]=bb-cb
			}
		}
		
		pdf(paste('GSE60055/output/CumDist/Introns.',it,'.',am,'.pdf',sep=''),width=4,height=4) ### FIGURE 5C 5D
		par(mar=c(5,5,2,2))
		
		cY=-0.2
		plot(rb,(-2*cY*bb)+cY,main='',type='l',lwd=2,pch=NA,col='gray',lty=2,
			 xlab=paste(am,'density (%)'),ylab='Cumulative fraction difference\nfrom background distribution',
			 xlim=c(0,xLIM[am]),ylim=c(cY,-cY),xaxt='n')
		axis(1,seq(0,xLIM[am],length.out=6),100*seq(0,xLIM[am],length.out=6))
		axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray50',col.axis='gray')
		mtext('Cumulative background distribution',4,col='gray',line=3)
		
		for (pn in c('Upregulated','Downregulated')){
			for (ei in c('Intron')){
				lines(rb,L[[paste(ei,pn)]],col=GCols[pn,ei],lwd=5)
			}
		}
		lines(c(0,xLIM[am]),c(0,0),col='black',lwd=2)
		dev.off()
	}
}



