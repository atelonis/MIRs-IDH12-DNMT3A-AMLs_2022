
library(amap)
library(gplots)
library(data.table)
library(VennDiagram)
source('Scripts/Heatmap3.r')

dir.create('output/Correlation/TAD.Plots',showWarnings=FALSE)


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

m=as.matrix(fread('output/Correlation/Step1/TadCorrelations.txt'))
m=m[-which(is.na(m[,3])),]

# Compute FDR and filter for FDR<5%
fdr=p.adjust(as.numeric(m[,4]),method='fdr')

m=cbind(m[,1:4],fdr,m[,5],m[,6])
colnames(m)=c('Methylation','Gene','Rho','PValue','FDR','Annotation','Distance')

m=m[order(as.numeric(m[,3])),]
posSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])>0)
negSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])<0)

n7=as.numeric(m[,7])

aaWhite=paste(col2hex('white'),'AA',sep='')
aaGray=paste(col2hex('gray'),'AA',sep='')
aaFire=paste(col2hex('firebrick'),'AA',sep='')
aaNavy=paste(col2hex('navy'),'AA',sep='')

pdf('output/Histogram.TadDistances.pdf',height=4,width=4) ### FIGURE 1C
par(mar=c(5,5,2,2))

hist(log10((n7/1000)+1),col='gray',main='Background',xlim=c(0,5),breaks=25,mgp=c(1.6,0.6,0),ylim=c(0,3000000),
	 xlab='Gene-mC distance (log10 kb)')

hist(log10((as.numeric(m[posSign,7])/1000)+1),col='#F8C1B8',main='Positive Correlations',xlim=c(0,5),breaks=25,ylim=c(0,60),
	 mgp=c(1.6,0.6,0),xlab='Gene-mC distance (log10 kb)',border='firebrick')
lines(log10(c(2,2)),c(0,60),col='black',lwd=1.5,lty=2)
lines(log10(c(500,500)),c(0,60),col='black',lwd=1.5,lty=2)

hist(log10((as.numeric(m[negSign,7])/1000)+1),col='#C6D6E3',main='Negative Correlations',xlim=c(0,5),breaks=25,ylim=c(0,2600),
	 mgp=c(1.6,0.6,0),xlab='Gene-mC distance (log10 kb)',border='navy')
lines(log10(c(2,2)),c(0,2500),col='black',lwd=1.5,lty=2)
lines(log10(c(500,500)),c(0,2500),col='black',lwd=1.5,lty=2)
dev.off()


# Add distance bining column
da=rep('Long',nrow(m))
da[n7<500000]='Intermediate'
da[n7<2000]='Overlap'

m=cbind(m,da)
colnames(m)[8]='DistanceBin'

sm=m[which(as.numeric(m[,5])<0.05),]
write.table(sm,'output/Correlation/Step2/TadCorrelations.txt',sep='\t',row.names=FALSE,quote=FALSE)

write.table(unique(m[,1]),file='output/Correlation/Step3/Cs.Background.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(unique(sm[,1]),file='output/Correlation/Step3/Cs.Significant.txt',sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


v=NULL
for (i in unique(sm[,'DistanceBin'])){
	v[[gsub('Overlap','Proximal',i)]]=unique(sm[which(sm[,'DistanceBin']==i),2])
}
venn.diagram(v,fill=c('pink','palegreen','skyblue'),file='output/TAD.GeneSetOverlap.png')  ### SUPPLEMENTARY FIGURE S1D


# Incude information about DE genes and DMCs
inSAM=read.table('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.IncludedInSAM.txt',sep='\t',stringsAsFactors=FALSE)[,1]
de=c(read.table('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.PositSign.txt',sep='\t',
				stringsAsFactors=FALSE,header=TRUE)[,3],
	 read.table('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.NegatSign.txt',sep='\t',
				stringsAsFactors=FALSE,header=TRUE)[,3])
deColumn=rep('NotIncluded',nrow(m))
deColumn[which(m[,2] %in% inSAM)]='NotDE'
deColumn[which(m[,2] %in% de)]='DE'

inMethylsig=read.table('output/MethylSig/Background.txt',stringsAsFactors=FALSE,sep='\t')[,1]
inMethylsig=paste(matrix(unlist(strsplit(inMethylsig,'-')),ncol=2,byrow=TRUE)[,1],'.NA0',sep='')
dm=c(read.table('output/MethylSig/Hypo.txt',sep='\t',stringsAsFactors=FALSE)[,1],
	 read.table('output/MethylSig/Hyper.txt',sep='\t',stringsAsFactors=FALSE)[,1])
dm=paste(matrix(unlist(strsplit(dm,'-')),ncol=2,byrow=TRUE)[,1],'.NA0',sep='')
dmColumn=rep('NotIncluded',nrow(m))
dmColumn[which(m[,1] %in% inMethylsig)]='NotDM'
dmColumn[which(m[,1] %in% dm)]='DM'

mdd=cbind(m,deColumn,dmColumn)
colnames(mdd)=c(colnames(m),c('DE','DM'))
write.table(mdd,file='output/Correlation/Step3/TAD.AnnotCorrs.txt',sep='\t',quote=FALSE,row.names=FALSE)
write.table(unique(mdd[,c(2,8,10)],file='output/Correlation/Step3/TAD.GeneCorrsDM.txt'),sep='\t',quote=FALSE,row.names=FALSE)
### Continued on Step015

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
	write.table(iGenes,file=paste('output/Correlation/Step3/TADsets',i,'txt',sep='.'),
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

pdf('output/Heatmap.TadCorrs.pdf') ### SUPPLEMENTARY FIGURE S1C
cDend=as.dendrogram(hcluster(mI/mU,method='spearman'))
heatmap.3(mI/mU,col=hc2,breaks=seq(0,0.4,length.out=length(hc2)+1),margins=c(18,18),
		  Rowv=cDend,Colv=cDend,KeyValueName='Jaccard index')
dev.off()


nonSign=which(as.numeric(m[,5])>0.05)
someNons=nonSign[round(seq(1,length(nonSign),length.out=1000))]

pdf('output/Correlation/TAD.VolcanoLike.pdf',height=4,width=4) ### FIGURE 1A
par(mar=c(5,5,2,2))
plot(c(-1,1),c(-log10(0.05),-log10(0.05)),xlim=c(-1,1),ylim=c(0,6),xlab='Spearman correlation coefficient',
	 ylab='-log10(FDR)',type='l',lwd=2,lty=2,col='gray50')
lines(as.numeric(m[someNons,3]),-log10(as.numeric(m[someNons,5])),col='black',lwd=1)
points(as.numeric(m[negSign,3]),-log10(as.numeric(m[negSign,5])),col=paste(col2hex('navy'),'AA',sep=''),pch=19)
points(as.numeric(m[posSign,3]),-log10(as.numeric(m[posSign,5])),col=paste(col2hex('firebrick'),'AA',sep=''),pch=19)
text(-0.35,5,paste(length(negSign),'negative\ncorrelations'),col='navy')
text(0.3,3,paste(length(posSign),'positive\ncorrelations'),col='firebrick')
text(0,-log10(0.05),'FDR=5%\n ',col='gray50',cex=0.75)
dev.off()


pdf('output/Correlation/TAD.NumberOfCorrs.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S1D
par(mar=c(10,10,2,2))
for (i in c('Overlap','Intermediate','Long')){
	cSM=sm[which(sm[,8]==i),]
	cGenes=tail(sort(table(cSM[,2])),10)
	names(cGenes)=matrix(unlist(strsplit(names(cGenes),'[|]')),ncol=2,byrow=TRUE)[,1]
	barplot(cGenes,horiz=TRUE,las=2,xaxt='n',xlab='Number of correlations',col='black',main=i,xlim=c(0,400))
	axis(1,c(0,200,400),c(0,200,400))
}
dev.off()

# Overlap with enhancers
pdf('output/TAD.EnrichEnhc.pdf',height=4,width=4) ### FIGURE 1F
par(mar=c(14,8,2,2))
enhC0=read.table('bed/BedIntersect.Cs.Enhancers.bed',stringsAsFactors=FALSE,sep='\t')
enhCs=unique(paste(enhC0[,1],enhC0[,3],'NA0',sep='.'))
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
HTs[,1]=log2(HTs[,1])
bp=barplot(HTs[,1],col='black',las=2,xaxt='n',xlab='Fold enrichment (log2)',main='Enhancers',
		   names.arg=c('Proximal','Intermediate','Long','TAD-wide'),xlim=c(-2,3),horiz=TRUE)
starHT=rep('*',nrow(HTs))
starHT[which(HTs[,2]>0.01)]=''
text(HTs[,1]+0.4,bp[,1],starHT,cex=2)
axis(1)
dev.off()



#################  PLOTS

# Read the data
d=as.matrix(fread('Datasets/CombinedData.txt'))
rownames(d)=read.table('Datasets/CombinedData.Rownames.txt',stringsAsFactors=FALSE,sep='\t')[,1]
colnames(d)=read.table('Datasets/Dataset.Colnames.txt',stringsAsFactors=FALSE,sep='\t')[,1]
d=d[,sort(colnames(d))]

fix_color=function(myColor){
	cRGB=col2rgb(myColor)/255
	return( rgb(cRGB[1],cRGB[2],cRGB[3],0.75) )
}



g=NULL
g$DNMT3A=colnames(d)[1:16]
g$IDH_DNMT3A=colnames(d)[26:36]
g$IDH=colnames(d)[17:25]

col_dnmt   = fix_color('orange')
col_double = fix_color('purple')
col_idh    = fix_color('deepskyblue')

myColors=c( rep(col_dnmt,16), rep(col_idh,9), rep(col_double,11) )

wGenes=c('DGKZ|8525','CPT1B|1375') ### FIGURE 1B
for (i in wGenes){
	
	print(i)
	
	di=d[i,]
	wMeth=sm[which(sm[,2]==i),]
	wMeth=wMeth[order(abs(as.numeric(wMeth[,3])),decreasing=TRUE),1]
	
	pdf(paste('output/Correlation/TAD.Plots/',gsub('[|]','.',i),'.pdf',sep=''),height=4,width=4)
	layout(matrix(1:4,nrow=2,ncol=2),widths=c(1,2),heights=c(2,1))
	par(mar=c(3,3,1,1))
	for (j in wMeth){
		jj=strsplit(j,'@')[[1]][1]
		dj=d[jj,]
		
		# top left
		par(mar=c(1,4,1,1))
		b3=NULL
		for (gi in names(g)){
			b3[[gi]]=d[i,g[[gi]]]
		}
		boxplot(b3,col=c(col_dnmt,col_double,col_idh),las=2,
				ylab=paste(strsplit(i,'[|]')[[1]][1],'(log2 Expression)'),xaxt='n',lwd=2)
		axis(1,c(1:3),labels=rep('',3))
		
		# bottom left
		par(mar=c(1,1,1,1))
		plot(1:10,col='white',axes=F,col.axis='white',xlab='',ylab='')
		par(srt=0)
		text(7.5,3,paste( 'rho=' , round(cor(di,dj,method='spearman'),3) ,sep=''))
		par(srt=-45)
		text(7.5,7.5,'DNMT3A', font=2,col=col_dnmt)
		text(8.5,8.5,'Double', font=2,col=col_double)
		text(9.5,9.5,'IDH',    font=2,col=col_idh)
		
		# top right
		par(mar=c(1,1,1,1))
		plot(dj,di,bg=myColors,pch=21,col.axis='white',cex=2,col=rgb(t(col2rgb(myColors)/255)))
		
		# bottom right
		par(mar=c(4,1,1,1))
		b2=NULL
		for (gi in names(g)){
			b2[[gi]]=d[jj,g[[gi]]]
		}
		boxplot(b2,col=c(col_dnmt,col_double,col_idh),
				xlab=paste(strsplit(j,'.NA0')[[1]],'(% Methylation)'),
				las=2,horizontal=TRUE,xaxt='n',yaxt='n',lwd=2)
		axis(1)
		axis(2,c(1:3),labels=rep('',3))
	}
	dev.off()
}



