
library(amap)
library(gplots)
library(data.table)
library(VennDiagram)
source('~/Scripts/Heatmap3.r')

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

m=as.matrix(fread('output/Correlation/Step1/CpGiCorrelations.txt'))

# Compute FDR and filter for FDR<5%
fdr=p.adjust(as.numeric(m[,4]),method='fdr')

m=cbind(m[,1:4],fdr,m[,5])
colnames(m)=c('Methylation','Gene','Rho','PValue','FDR','Annotation')

m=m[order(as.numeric(m[,3])),]
posSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])>0.5)
negSign=which(as.numeric(m[,5])<0.05 & as.numeric(m[,3])<(-0.5))

sm=m[which(as.numeric(m[,5])<0.05 & abs(as.numeric(m[,3]))>0.5),]
write.table(sm,'output/Correlation/Step2/CpGiCorrelations.txt',sep='\t',row.names=FALSE,quote=FALSE)

cpgiGenes=sort(unique(sm[,2]))


####################################
### Compare with mC correlations ###

mc=read.table('output/Correlation/Step2/TadCorrelations.txt',stringsAsFactors=FALSE,header=TRUE)
mcGenes=sort(unique(mc[,2]))


allGenes=read.table('Datasets/CombinedData.Rownames.txt',stringsAsFactors=FALSE)[,1]
allGenes=grep('[|]',allGenes,value=TRUE)

v=NULL
v$mC=mcGenes
v$CpGi=cpgiGenes

venn.diagram(v,fill=c('pink','skyblue'),file='output/Venn.Compare_mC_CpGi.png') ### SUPPLEMENTARY FIGURE S1F

### Remove log files of Venn diagrams
for (i in list.files('output/',full.names=TRUE,pattern='log')){
	ii=strsplit(i,'[.]')[[1]]
	if (ii[length(ii)]=='log'){
		file.remove(i)
	}
}





