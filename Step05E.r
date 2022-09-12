
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

m=as.matrix(fread('output/Correlation/Step1/NonMutantCorrs.txt'))
m=m[-which(is.na(m[,3])),]

# Compute FDR and filter for FDR<5%
fdr=p.adjust(as.numeric(m[,4]),method='fdr')

m=cbind(m[,1:4],fdr,m[,5],m[,6])
colnames(m)=c('Methylation','Gene','Rho','PValue','FDR','Annotation','Distance')


###########################################################
### Compare with the correlations of the mutant samples ### 

rownames(m)=paste(m[,1],m[,2],sep='--')

mutCorrs=read.table('output/Correlation/Step2/TadCorrelations.txt',header=TRUE,stringsAsFactors=FALSE)
rownames(mutCorrs)=paste(mutCorrs[,1],mutCorrs[,2],sep='--')
mutCorrs=mutCorrs[intersect(rownames(mutCorrs),rownames(m)),]
mnm=cbind(abs(as.numeric(mutCorrs[,3])),abs(as.numeric(m[rownames(mutCorrs),3])))

pdf('output/Barplot.CompareToWtWt.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S1H
par(mar=c(5,12,2,2))
boxplot(mnm[,1],mnm[,2],names=c('DNMT3A\nor IDH1/2','WT/WT'),col='gray',ylim=c(0,1),las=2,yaxt='n',
		mgp=c(1.6,0.6,0),ylab='Spearman correlation coefficient')
axis(2,mgp=c(1.6,0.6,0))
dev.off()


