
library(car)
library(samr)
library(amap)
library(gplots)
library(ggplot2)
library(ggpubr)
library(dendextend)
library(DiscriMiner)

source('Scripts/Heatmap3.r')

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

Run_SAM=function(curda,y,myoutput){
	
	# exlude features with zero SD within at least one of the groups
	toExclude=unique( c(names(which(apply(curda[,which(y==1)],1,sd)==0)), 
					    names(which(apply(curda[,which(y==2)],1,sd)==0)) ))
	toInclude=setdiff(rownames(curda),toExclude)
	curda=curda[toInclude,]
	write.table(rownames(curda),sep='\t',file=paste(myoutput,'.SAM.IncludedInSAM.txt',sep=''),
				row.names=FALSE,col.names=FALSE,quote=FALSE)
	
	# Run R
	mylist=list(x=curda,y=y,geneid=rownames(curda),logged2=FALSE)
	mySAM=samr(mylist,resp.type="Two class unpaired",nperm=5000)
	
	# Calculate Delta
	delta.table = samr.compute.delta.table(mySAM)
	delta = -1
	for(i in 1:nrow(delta.table)){
		if(!((delta.table[i,5]>0.05) | is.na(delta.table[i,5]))){
			delta = delta.table[i,1]
			break()
		}
	}
	if (delta==-1){
		delta=delta.table[i,1]
		}
	
	# Save observed score
	write.table(mySAM$tt,file=paste(myoutput,'.SAM.Score.txt',sep=''),col.names=FALSE,quote=FALSE,sep='\t')
	
	# Save R graph
	write.table(delta,file=paste(myoutput,'.SAM.DeltaValue.txt',sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE)
	pdf(paste(myoutput,'.SAM.Graph.delta.pdf',sep=''),height=4,width=4)
	par(mar=c(5,5,2,2))
	samr.plot(mySAM,delta)
	dev.off()
	
	# Calculate and save significant features
	siggenes.table = samr.compute.siggenes.table(mySAM,delta,mylist,delta.table)
	write.table(siggenes.table$genes.up,sep='\t',file=paste(myoutput,'.SAM.PositSign.txt',sep=''),row.names=FALSE)
	write.table(siggenes.table$genes.lo,sep='\t',file=paste(myoutput,'.SAM.NegatSign.txt',sep=''),row.names=FALSE)
}

get_cells=function(gCells){
	rCells=NULL
	for (i in gCells){
		rCells=c(rCells,strsplit(i,'_')[[1]][1])
	}
	rCells=table(rCells)
	return(rCells)
}

get_FE_table=function(ctCells,D){
	p=matrix(0,nrow=length(names(D))+1,ncol=length(ctCells))
	rownames(p)=c('Background','Comparison')
	colnames(p)=sort(names(ctCells))
	p[1,]=ctCells[colnames(p)]
	p[2,]=D[colnames(p)]
	p[is.na(p)]=0
	
	FE=NULL
	PV=NULL
	for (j in colnames(p)){
		cFE = (p[i,j]/sum(p[i,])) / (p['Background',j]/sum(p['Background',]))
		FE=c(FE,cFE)
		if (cFE > 1){
			PV=c(PV,1-phyper( p[i,j]-1, p['Background',j], sum(p['Background',])-p['Background',j], sum(p[i,])) )
		} else {
			PV=c(PV,phyper( p[i,j], p['Background',j], sum(p['Background',])-p['Background',j], sum(p[i,])) )
		}
	}
	fdr=p.adjust(PV,method='fdr')
	cTable=cbind(FE,fdr)
	rownames(cTable)=colnames(p)
	
	return(cTable)
}
############################################

lookup=read.table('ID_LookupTable.txt',stringsAsFactors=FALSE)
rownames(lookup)=lookup[,2]
mirDensity=mirDensity[which(names(mirDensity) %in% rownames(lookup))]
names(mirDensity)=lookup[names(mirDensity),1]

E=read.table('Hematopoiesis/Dataset.75percMostExpr.txt',check.names=FALSE)
A=read.table('Microarrays/Dataset.50percMostExpr.txt')

coms=sort(intersect(rownames(E),rownames(A)))
coms=coms[which(coms %in% names(mirDensity))]
comDensity=mirDensity[coms]

v=NULL # The idea here was to check whether genes dense in MIRs have better/worse predictions than the 
v[['All']]     = coms # genes with the low density in MIRs. This didn't work, so I will need to 'clean' 
v[['MostMIR']] = coms[comDensity >= quantile(comDensity,0.75) ] # this code a bit, at some point.
v[['LessMIR']] = coms[comDensity <= quantile(comDensity,0.75) ]

Vdnm=matrix(0,nrow=3,ncol=2)
rownames(Vdnm)=names(v)
colnames(Vdnm)=c('HSC','Mono')
Vidh=Vdnm
Vdou=Vdnm
Pdnm=Vdnm
Pidh=Vdnm
Pdou=Vdnm

RE=apply(E,2,rank)
RA=apply(A,2,rank)

# Run PCA on ranking set
aCells=NULL
for (i in colnames(RE)){
	aCells=c(aCells,strsplit(i,'_')[[1]][1])
}

uCells=sort(unique(aCells))
cpTable=read.table('ColorsPch.txt',stringsAsFactors=FALSE)
cpTable[,1]=gsub('E9EC6B','#E9EC6B',cpTable[,1])
cpTable[,4]=gsub('E9EC6B','#E9EC6B',cpTable[,4])
cpTable[,1]=paste(col2hex(cpTable[,1]),'CC',sep='')
aColors=cpTable[aCells,1]
aPch=cpTable[aCells,2]
aLins=cpTable[aCells,3]
aColLins=cpTable[aCells,4]

mypca=prcomp(t(RE))
v1=((mypca$sdev^2)[1])/(sum(mypca$sdev^2))
v2=((mypca$sdev^2)[2])/(sum(mypca$sdev^2))
v3=((mypca$sdev^2)[3])/(sum(mypca$sdev^2))
pc1=paste("PC1 (",round(100*v1,1),"%)",sep='')
pc2=paste("PC2 (",round(100*v2,1),"%)",sep='')
pc3=paste("PC3 (",round(100*v3,1),"%)",sep='')
lims=max(abs(mypca$x[,c(1,2)]))*1.1
lims=c(-lims,lims)

### project the AML samples on this PC space
# center the data using the same centering as in the normal samples
pAML= t( RA - mypca$center )
pAML = pAML %*% mypca$rotation
amlCols=rep('orange',nrow(pAML))
amlCols[grep('IDH',rownames(pAML))]='deepskyblue'
amlCols[grep('IDH_DNMT3A',rownames(pAML))]='purple'
amlCols=paste(col2hex(amlCols),'CC',sep='')

PCnames='Lineage'
for (i in c(1:ncol(mypca$x))){
	cVar=((mypca$sdev^2)[i])/(sum(mypca$sdev^2))
	PCnames=c(PCnames,paste("PC",i," (",round(100*cVar,1),"%)",sep=''))
}
tw=cbind(aLins,mypca$x)
tw=rbind(tw, cbind(gsub('_[0-9]...','_AML',rownames(pAML)),pAML))
colnames(tw)=PCnames
rownames(tw)=c(rownames(mypca$x),rownames(pAML))
write.table(tw,file='output/PCA_Rank_Projections.txt',sep='\t')

pdf('output/PCA_Rank.pdf',height=4,width=4) ### FIGURE 6A
par(mar=c(5,5,1,1))

### Plot normal ellipses  +  AML samples
plot(0,0,pch=NA,col=NA,xlab=pc1,ylab=pc2,xlim=lims,ylim=lims*v1/v2)
myEllipses=dataEllipse(mypca$x[,1],mypca$x[,2],group.labels=NA,
					   groups=as.factor(aLins),draw=FALSE,levels=0.8)
for (e in c('Lymphoid','Myeloid','Erythroid','Progenitor','HSC')){
	eCol=cpTable[grep(e,cpTable[,3])[1],4]
	polygon(myEllipses[[e]][,1],myEllipses[[e]][,2],border=eCol,lwd=2,
			col=paste(col2hex(eCol),'CC',sep=''))
}
for (e in names(amlEllipses)){
	if (e=='IDH'){
		eCol='deepskyblue'
	} else if (e=='DNMT3A'){
		eCol='orange'
	} else {
		eCol='purple'
	}
	
	cSamples=grep(paste('^',e,'_[0-9]',sep=''),rownames(pAML),value=TRUE)
	medSam=c( median(pAML[cSamples,1]), median(pAML[cSamples,2]))
	medMatX=rbind(pAML[cSamples,1],rep(medSam[1],length(cSamples)))
	medMatY=rbind(pAML[cSamples,2],rep(medSam[2],length(cSamples)))
	
	matplot(medMatX,medMatY,type='l',col=eCol,lwd=2,lty=1,add=TRUE)
	points(pAML[cSamples,1],pAML[cSamples,2],pch=21,cex=0.5,col='black',bg=eCol) # paste(col2hex(eCol),'FF',sep=''))
	points(medSam[1],medSam[2],pch=22,col='black',bg=eCol,cex=1.25)
}

lims20=c(-20000,5000)
lineLims=c(lims20,rev(lims20),lims20[1])
lines(lineLims,rev(lineLims)*v1/v2,col='black',lty=2,lwd=1)

# Plot normal ellipses  +  AML samples  ZOOMED
plot(0,0,pch=NA,col=NA,xlab=pc1,ylab=pc2,xlim=lims20,ylim=lims20*v1/v2,xaxs='i',yaxs='i')
for (e in c('Lymphoid','Myeloid','Erythroid','Progenitor','HSC')){
	eCol=cpTable[grep(e,cpTable[,3])[1],4]
	polygon(myEllipses[[e]][,1],myEllipses[[e]][,2],border=eCol,lwd=2,
			col=paste(col2hex(eCol),'CC',sep=''))
}
for (e in names(amlEllipses)){
	if (e=='IDH'){
		eCol='deepskyblue'
	} else if (e=='DNMT3A'){
		eCol='orange'
	} else {
		eCol='purple'
	}
	
	cSamples=grep(paste('^',e,'_[0-9]',sep=''),rownames(pAML),value=TRUE)
	medSam=c( median(pAML[cSamples,1]), median(pAML[cSamples,2]))
	medMatX=rbind(pAML[cSamples,1],rep(medSam[1],length(cSamples)))
	medMatY=rbind(pAML[cSamples,2],rep(medSam[2],length(cSamples)))
	
	matplot(medMatX,medMatY,type='l',col=eCol,lwd=2,lty=1,add=TRUE)
	points(pAML[cSamples,1],pAML[cSamples,2],pch=21,cex=1.25,col='black',bg=eCol) # paste(col2hex(eCol),'FF',sep=''))
	points(medSam[1],medSam[2],pch=22,col='black',bg=eCol,cex=2)
}

# 10. Legends
plot(0,0,pch=NA,col=NA,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',legend=rownames(cpTable),col=cpTable[,1],pch=cpTable[,2],ncol=3,
	   pt.cex=1.5,cex=0.75)
plot(0,0,pch=NA,col=NA,xlab='',ylab='',xaxt='n',yaxt='n')
legend('center',legend=c('IDH','DNMT3A','Double'),col=cpTable[,1],pch=cpTable[,2],
	   pt.cex=1.5,cex=0.75)


#################################################################

# `m` is a matrix with the distances of each AML sample to each Normal sample
m=matrix(0,nrow=ncol(RA),ncol=ncol(RE))
rownames(m)=colnames(RA)
colnames(m)=colnames(RE)
for (i in colnames(RA)){
	for (j in colnames(RE)){
		m[i,j]=dist(rbind(RA[,i],RE[,j]),method='manhattan')
	}
}
rownames(m)=gsub('IDH_DNMT3A','Double',rownames(m))
idh=grep('IDH',rownames(m),value=TRUE)
dnm=grep('DNMT3A',rownames(m),value=TRUE)
dou=grep('Double',rownames(m),value=TRUE)

myCells=NULL
for (i in colnames(m)){
	myCells=c(myCells,strsplit(i,'_')[[1]][1])
}
uCells=table(myCells)

m0=m

m=m-apply(m,1,median)
write.table(t(m),file='output/MatrixDistance.txt',sep='\t')


############
### SAMs ###
m=t(m)

### Double vs DNMT3A
cM=m[,c(dnm,dou)]
y=c(rep(1,length(dnm)),rep(2,length(dou)))
Run_SAM(cM,y,'output/SAM/Distance.Double_vs_DNMT3A') ### FIGURE 6B
cPosit=read.table('output/SAM/Distance.Double_vs_DNMT3A.SAM.PositSign.txt',
				  sep='\t',stringsAsFactors=FALSE,header=TRUE)[,3]
P=get_cells(cPosit)
dnmPosit=cPosit

### Double vs IDH
cM=m[,c(idh,dou)]
y=c(rep(1,length(idh)),rep(2,length(dou)))
Run_SAM(cM,y,'output/SAM/Distance.Double_vs_IDH') ### FIGURE 6C
cPosit=read.table('output/SAM/Distance.Double_vs_IDH.SAM.PositSign.txt',
				  sep='\t',stringsAsFactors=FALSE,header=TRUE)[,3]
cNegat=read.table('output/SAM/Distance.Double_vs_IDH.SAM.NegatSign.txt',
				  sep='\t',stringsAsFactors=FALSE,header=TRUE)[,3]
PI=get_cells(cPosit)
NI=get_cells(cNegat)
idhPosit=cPosit	

FE=NULL
FE[['DNMT3A_Up']] = get_FE_table(uCells,P)$All
FE[['IDH_Up']]    = get_FE_table(uCells,PI)$All
FE[['IDH_Down']]  = get_FE_table(uCells,NI)$All

for (i in names(FE)){
	cFE=FE[[i]]
	cFE=sort(cFE[cFE[,1]>2 & cFE[,2]<0.05,1])
	pdf(paste('output/Enrichment.SAM_',i,'.pdf',sep=''),height=5,width=5) ### FIGURE 6D 6E 6F
	par(mar=c(15,10,2,2))
	barplot(cFE,col='firebrick',border=NA,las=2,xaxt='n',xlab='Fold enrichment',horiz=TRUE,xlim=c(0,5))
	axis(1)
	dev.off()
}


lins=c('Erythroid','Myeloid','Lymphoid')

doubleLins=c(0,0,0)
names(doubleLins)=lins

for (i in lins){
	doubleLins[i]=sum(NI$All[intersect(names(NI$All),rownames(cpTable)[which(cpTable[,3]==i)])])
}

backLins=table(aLins)[lins]

ht=NULL
for (i in lins){
	ht=rbind(ht,hypergeom.test( doubleLins[i], sum(doubleLins), backLins[i], sum(backLins)))
}
	
pdf(paste('output/Enrichment.SAM_IDH_Down.Lineages.pdf',sep=''),height=5,width=5) ### FIGURE 6G
par(mar=c(20,15,2,2))
barplot(log2(ht[,1]),col='firebrick',border=NA,las=2,xaxt='n',xlab='Fold enrichment (log2)',
		names.arg=lins,horiz=TRUE,xlim=c(-1,2))
axis(1)
dev.off()




