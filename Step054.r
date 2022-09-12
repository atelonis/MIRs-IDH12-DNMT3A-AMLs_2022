
library(samr)
library(amap)
library(gplots)
library(data.table)
library(VennDiagram)
source('~/Scripts/Heatmap3.r')

cutdf=function(myVector,myDelim,myField){
	myR=matrix(unlist(strsplit(myVector,myDelim)),ncol=2,byrow=TRUE)[,myField]
	return(myR)
}

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
	toInclude=unique( names(which(apply(curda[,which(y==1)],1,sd)>0)), 
					  names(which(apply(curda[,which(y==2)],1,sd)>0)) )
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
			delta = delta.table[i,1];d1=delta.table[i,]
			break()}
		}
	if (delta==-1){
		delta=delta.table[i,1]
		}
	
	# Save observed score
	write.table(mySAM$tt,file=paste(myoutput,'.SAM.Score.txt',sep=''),col.names=FALSE,quote=FALSE,sep='\t')
	
	# Save R graph
	write.table(delta,file=paste(myoutput,'.SAM.DeltaValue.txt',sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE)
	pdf(paste(myoutput,'.SAM.Graph.delta.pdf',sep=''),height=5,width=5)
	par(mar=c(5,5,2,2))
	samr.plot(mySAM,delta)
	dev.off()
	
	# Calculate and save significant features
	siggenes.table = samr.compute.siggenes.table(mySAM,delta,mylist,delta.table)
	write.table(siggenes.table$genes.up,sep='\t',file=paste(myoutput,'.SAM.PositSign.txt',sep=''),row.names=FALSE)
	write.table(siggenes.table$genes.lo,sep='\t',file=paste(myoutput,'.SAM.NegatSign.txt',sep=''),row.names=FALSE)
}


######################
### PART 1  |  DE  ###
######################

dir.create('TCGA/output/SAM',showWarnings=FALSE)

d=read.table('TCGA/RnaData.txt',check.names=FALSE)

muts=NULL
for (i in read.table('Mutations.txt',stringsAsFactors=FALSE)[,1]){
	i=strsplit(i,'[.]')[[1]]
	muts=rbind(muts,c(i[1],paste(i[1],i[2],sep='.')))
}
rownames(muts)=muts[,1]
colnames(d)=muts[colnames(d),2]

idh=grep('IDH',colnames(d),value=TRUE)
dnm=grep('DNMT3A',colnames(d),value=TRUE)
dou=grep('Double',colnames(d),value=TRUE)


y=c(rep(2,length(idh)),rep(1,length(dnm)))
Run_SAM(as.matrix(d[,c(idh,dnm)]),y,'TCGA/output/SAM/IDH_vs_DNMT3A')




#######################################
### PART 2  |  Analyze correlations ###
#######################################


dir.create('TCGA/output/Correlation/Step2',showWarnings=FALSE)
dir.create('TCGA/output/Correlation/Step3',showWarnings=FALSE)

m=as.matrix(fread('TCGA/output/Correlation/Step1/Correlations.txt'))
m=m[-which(is.na(m[,3])),]

# Compute FDR and filter for FDR<5%
fdr=p.adjust(as.numeric(m[,4]),method='fdr')

m=cbind(m[,1:4],fdr,m[,5],m[,6])
colnames(m)=c('Methylation','Gene','Rho','PValue','FDR','Annotation','Distance')

m=m[order(as.numeric(m[,3])),]
posSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])>0)
negSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])<0)


n7=as.numeric(gsub('Upstream_','',gsub('Downstream_','',gsub('GeneBody','0',m[,7]))))

aaWhite=paste(col2hex('white'),'AA',sep='')
aaGray=paste(col2hex('gray'),'AA',sep='')
aaFire=paste(col2hex('firebrick'),'AA',sep='')
aaNavy=paste(col2hex('navy'),'AA',sep='')

pdf('TCGA/output/Histogram.TadDistances.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S1B
par(mar=c(5,5,2,2))
hist(log10(n7+1),col='gray',main='Background',xlim=c(0,5),breaks=25,mgp=c(1.6,0.6,0),ylim=c(0,500000),
	 xlab='Gene-mC distance (log10 kb)')
hist(log10(n7[posSign]+1),col='#F8C1B8',main='Positive Correlations',xlim=c(0,5),breaks=25,ylim=c(0,800),mgp=c(1.6,0.6,0),
	 xlab='Gene-mC distance (log10 kb)',border='firebrick')
lines(log10(c(2,2)),c(0,800),col='black',lwd=1.5,lty=2)
lines(log10(c(500,500)),c(0,800),col='black',lwd=1.5,lty=2)
hist(log10(n7[negSign]+1),col='#C6D6E3',main='Negative Correlations',xlim=c(0,5),breaks=25,ylim=c(0,800),mgp=c(1.6,0.6,0),
	 xlab='Gene-mC distance (log10 kb)',border='navy')
lines(log10(c(2,2)),c(0,800),col='black',lwd=1.5,lty=2)
lines(log10(c(500,500)),c(0,800),col='black',lwd=1.5,lty=2)
dev.off()

# Add distance bining column
da=rep('Long',nrow(m))
da[n7<500000]='Intermediate'
da[n7<2000]='Overlap'

m=cbind(m,da)
colnames(m)[8]='DistanceBin'

v=NULL
for (i in unique(sm[,'DistanceBin'])){
	v[[gsub('Overlap','Proximal',i)]]=unique(sm[which(sm[,'DistanceBin']==i),2])
}
venn.diagram(v,fill=c('pink','palegreen','skyblue'),file='TCGA/output/TAD.GeneSetOverlap.png') ### SUPPLEMENTARY FIGURE S1C

sm=m[which(as.numeric(m[,5])<0.05),]
write.table(sm,'TCGA/output/Correlation/Step2/TadCorrelations.txt',sep='\t',row.names=FALSE,quote=FALSE)

cSpears=NULL
for (i in c('Overlap','Intermediate','Long')){
	mGenes=table(m[which(m[,8]==i),2])
	sGenes=table(sm[which(sm[,8]==i),2])
	cSpears=c(cSpears,cor(sGenes,mGenes[names(sGenes)],method='spearman'))
}
mGenes=table(m[,2])
sGenes=table(sm[,2])
cSpears=c(cSpears,cor(sGenes,mGenes[names(sGenes)],method='spearman'))
names(cSpears)=c('Proximal','Intermediate','Long','TAD-wide')
cSpears=rev(cSpears)

pdf('TCGA/output/TAD.Corr_NoCytos_NoSigns.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S1D
par(mar=c(14,8,2,5))
barplot(cSpears,xlim=c(-1,1),col='black',xlab='Spearman correlation coefficient',
		names.arg=names(cSpears),xaxt='n',las=2,horiz=TRUE,mgp=c(1.6,0.6,0))
axis(1,mgp=c(1.6,0.6,0))
dev.off()


# Incude information about DE genes and DMCs
inSAM=cutdf(read.table('TCGA/output/SAM/IDH_vs_DNMT3A.SAM.IncludedInSAM.txt',sep='\t',stringsAsFactors=FALSE)[,1],'[|]',1)
de=cutdf(c(read.table('TCGA/output/SAM/IDH_vs_DNMT3A.SAM.PositSign.txt',sep='\t',
				stringsAsFactors=FALSE,header=TRUE)[,3],
		   read.table('TCGA/output/SAM/IDH_vs_DNMT3A.SAM.NegatSign.txt',sep='\t',
				stringsAsFactors=FALSE,header=TRUE)[,3]),'[|]',1)
deColumn=rep('NotIncluded',nrow(m))
deColumn[which(cutdf(m[,2],'[|]',2) %in% inSAM)]='NotDE'
deColumn[which(cutdf(m[,2],'[|]',2) %in% de)]='DE'

withTRFs=cutdf(unique(read.table('TCGA/tRF_mRNA.LAML.txt',stringsAsFactors=FALSE)[,2]),'[|]',1)
trfColumn=rep('NotCorrelated',nrow(m))
trfColumn[which(cutdf(m[,2],'[|]',2) %in% withTRFs)]='tRF'

mdt=cbind(m,deColumn,trfColumn)
colnames(mdt)=c(colnames(m),c('DE','tRF'))
write.table(mdt,file='TCGA/output/Correlation/Step3/TAD.AnnotCorrs_DE_tRF.txt',sep='\t',quote=FALSE,row.names=FALSE)
write.table(unique(mdt[,c(2,8,10)],file='TCGA/output/Correlation/Step3/TAD.GeneCorrsTRF.txt'),sep='\t',quote=FALSE,row.names=FALSE)

### Now, let's analyze significants

bV=table(m[,8])[c('Overlap','Intermediate','Long')]
pV=table(m[posSign,8])[c('Overlap','Intermediate','Long')]
nV=table(m[negSign,8])[c('Overlap','Intermediate','Long')]

chiP=chisq.test(rbind(bV-pV,pV))
chiN=chisq.test(rbind(bV-nV,nV))

# Jaccard and Enrichment matrices
ul2=length(unique(m[,2]))
pn=list(Posit=unique(which(as.numeric(sm[,3])>0)), Negat=unique(which(as.numeric(sm[,3])<0)))

mI=matrix(0,nrow=6,ncol=6)
rownames(mI)=c(paste('Posit',names(bV),sep='_'),paste('Negat',names(bV),sep='_'))
colnames(mI)=c(paste('Posit',names(bV),sep='_'),paste('Negat',names(bV),sep='_'))
mU=mI
mF=mI
mC=mI
for (i in rownames(mI)){
	ii=strsplit(i,'_')[[1]]
	iGenes=unique( sm[ intersect(pn[[ii[1]]],unique(which(sm[,8]==ii[2]))), 2 ] )
	write.table(iGenes,file=paste('TCGA/output/Correlation/Step3/TADsets',i,'txt',sep='.'),
				col.names=FALSE,row.names=FALSE,quote=FALSE)
	for (j in colnames(mI)){
		jj=strsplit(j,'_')[[1]]
		jGenes=unique( sm[ intersect(pn[[jj[1]]],unique(which(sm[,8]==jj[2]))), 2 ] )
		mI[i,j]=length(intersect(iGenes,jGenes))
		mU[i,j]=length(union(iGenes,jGenes))
		mC[i,j]=length(iGenes)
		mF[i,j]= (mI[i,j]/length(jGenes)) / (length(iGenes)/ul2)
	}
}
diag(mF)=1

rownames(mI)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))
colnames(mI)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))
rownames(mU)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))
colnames(mU)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))
rownames(mF)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))
colnames(mF)=gsub('Overlap','Proximal',paste(rownames(mC),' (n=',mC[,1],')',sep=''))

lmf=log2(mF)
lmf[is.infinite(lmf)]=0

hc1=colorRampPalette(c('navy','white','firebrick'))(51)
hc2=colorRampPalette(c('white','darkgreen'))(51)
pdf('TCGA/output/Heatmap.TadCorrs.pdf')
cDend=as.dendrogram(hcluster(mI/mU,method='spearman'))
heatmap.3(mI/mU,col=hc2,breaks=seq(0,0.4,length.out=length(hc2)+1),margins=c(18,18),
		  Rowv=cDend,Colv=cDend,KeyValueName='Jaccard index')
dev.off()


nonSign=which(as.numeric(m[,5])>0.05)
someNons=nonSign[round(seq(1,length(nonSign),length.out=1000))]

pdf('TCGA/output/TAD.CorrsVolcanoLike.pdf',height=4,width=4) ### FIGURE 1A
par(mar=c(5,5,2,2))
plot(c(-1,1),c(-log10(0.05),-log10(0.05)),xlim=c(-1,1),ylim=c(0,10),xlab='Spearman correlation coefficient',
	 ylab='-log10(FDR)',type='l',lwd=2,lty=2,col='gray50')
lines(as.numeric(m[someNons,3]),-log10(as.numeric(m[someNons,5])),col='black',lwd=1)
points(as.numeric(m[negSign,3]),-log10(as.numeric(m[negSign,5])),col=paste(col2hex('navy'),'AA',sep=''),pch=19)
points(as.numeric(m[posSign,3]),-log10(as.numeric(m[posSign,5])),col=paste(col2hex('firebrick'),'AA',sep=''),pch=19)
text(-0.35,6,paste(length(negSign),'negative\ncorrelations'),col='navy')
text(0.3,4,paste(length(posSign),'positive\ncorrelations'),col='firebrick')
text(0,-log10(0.05),'FDR=5%\n ',col='gray50',cex=0.75)
dev.off()

pdf('TCGA/output/TAD.NumberOfCorrs.pdf',height=4,width=4) 
par(mar=c(10,10,2,2))
for (i in c('Overlap','Intermediate','Long')){
	cSM=sm[which(sm[,8]==i),]
	cGenes=tail(sort(table(cSM[,2])),10)
	names(cGenes)=matrix(unlist(strsplit(names(cGenes),'[|]')),ncol=2,byrow=TRUE)[,2]
	barplot(cGenes,horiz=TRUE,las=2,xaxt='n',xlab='Number of correlations',col='black',main=i,xlim=c(0,160))
	axis(1)
}
dev.off()

#### Overlap with enhancers
pdf('TCGA/output/TAD.EnrichEnhc.pdf',height=4,width=4) ### FIGURE 1F
par(mar=c(10,14,2,2))
for (enh in c('Enhancers')){
	enhC0=read.table(paste('BedIntersect.Cs',enh,'bed',sep='.'),stringsAsFactors=FALSE,sep='\t')
	enhCs=unique(enhC0[,4])
	HTs=NULL
	for (oil in c('Overlap','Intermediate','Long','TAD-wide')){
		if (oil!='TAD-wide'){
			cBack=unique(m[which(m[,'DistanceBin']==oil),1])
			cSign=unique(sm[which(sm[,'DistanceBin']==oil),1])
		} else {
			cBack=unique(m[,1])
			cSign=unique(sm[,1])
		}
		eBack=intersect(cBack,enhCs)
		eSign=intersect(cSign,enhCs)
		cHT=hypergeom.test(length(eSign),length(cSign),length(eBack),length(cBack))
		HTs=rbind(HTs,cHT)
	}
	
	bp=barplot(HTs[,1],col='black',las=2,yaxt='n',ylab='Fold enrichment',main=enh,
			   names.arg=c('Proximal','Intermediate','Long','TAD-wide'),ylim=c(0,3))
	starHT=rep('*',nrow(HTs))
	starHT[which(HTs[,2]>0.01)]=''
	text(bp[,1],HTs[,1]*1.2,starHT,cex=2)
	axis(2)
}
dev.off()

########################################
### Evolutionary conservation of mCs ###
########################################

ev=NULL
for (i in c(1:22,'Y')){
	ev=rbind(ev,fread(paste('TCGA/bed/BedIntersect.TadCs.EvolCons_chr',i,'.bed',sep='')))
}
write.table(ev,file='TCGA/bed/BedIntersect.TadCs.EvolCons_Cs.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
evC=as.matrix(fread('TCGA/bed/BedIntersect.TadCs.EvolCons_Cs.bed'))
rownames(evC)=evC[,4]

corTad=read.table('TCGA/output/Correlation/Step2/TadCorrelations.txt',header=TRUE,stringsAsFactors=FALSE)

evBackground=ecdf(as.numeric(evC[,8]))

Xs=seq(-10,10,length.out=101)
dBins=c('Overlap','Intermediate','Long')

binCol=c('#F7CD00','#C19117','#4A6A1D')
names(binCol)=dBins

wKS=NULL
pdf('TCGA/output/TAD.EvolCons_Cs.pdf',height=4,width=4) ### FIGURE 2E
par(mar=c(4,4,3,3))
plot(0,0,pch=NA,xlim=c(-10,10),ylim=c(-0.1,0.1),xlab='Evolutionary conservation (phyloP score)',
	 ylab='Cumulative fraction difference\nfrom background distribution',mgp=c(1.6,0.6,0))
lines(Xs,(evBackground(Xs)*0.2)-0.1,lwd=2,lty=2,col='gray60')
for (i in dBins){
	print(i)
	cC=intersect(rownames(evC),unique(corTad[which(corTad[,8]==i),1]))
	cEV=as.numeric(evC[cC,8])
	cKS=ks.test(cEV,evBackground)
	lines(Xs,evBackground(Xs)-ecdf(cEV)(Xs),col=binCol[i],lwd=5)
	print(cKS)
	wKS=rbind(wKS,c(i,cKS$p.value))
}
axis(4,c(-0.1,-0.05,0,0.05,0.1),c(0,0.25,0.5,0.75,1),col='gray60',col.axis='gray60',mgp=c(1.6,0.6,0))
mtext('Cumulative background distribution',4,col='gray60',line=1.6)
lines(c(-10,10),c(0,0),col='black',lwd=3)
dev.off()

write.table(wKS,'output/KSTest.EvolConsCs.txt',sep='\t',col.names=F,row.names=F,quote=F)





