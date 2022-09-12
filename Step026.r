
library(gplots)
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

d=read.table('TADs.GenesDeDm.txt',check.names=FALSE)
d=d[,setdiff(colnames(d),c("Satellite/acro",'Unknown','Simple_repeat','Low_complexity'))]

# Distinguish between TADs that CAN have correlations and those that CANNOT - i.e. as background
d0=d[d[,2]==0,]
d1=d[d[,2]>0,]
d10=d[d[,3]>=10 & d[,4]>=10,]

### Check enrichment of TADs in correlations, DE genes and DM cytosines
### These numbers differ from colSums as there can be overlap with more than one TADs
countDEgenes=890 # output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.*Sign.txt
countDMcytos=104730 # output/MethylSig/Hyp*txt
countBackGenes=9972 # output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.IncludedInSAM.txt
countBackCytos=1524706 # output/MethylSig/Background.txt

enrCorrs=NULL
enrDEgen=NULL
enrDMcyt=NULL
for (i in rownames(d)){
	enrCorrs=rbind(enrCorrs,hypergeom.test(d[i,'SignCorrs'],sum(d[,'SignCorrs']),d[i,'BackCorrs'],sum(d[,'BackCorrs'])))
	enrDEgen=rbind(enrDEgen,hypergeom.test(d[i,'DeGenes'],countDEgenes,d[i,'BackGenes'],countBackGenes))
	enrDMcyt=rbind(enrDMcyt,hypergeom.test(d[i,'DmCytos'],countDMcytos,d[i,'BackCytos'],countBackCytos))
}
rownames(enrCorrs)=rownames(d)
rownames(enrDEgen)=rownames(d)
rownames(enrDMcyt)=rownames(d)

eCorrs=enrCorrs[rownames(d10),]
eDEgen=enrDEgen[rownames(d10),]
eDMcyt=enrDMcyt[rownames(d10),]

eCorrs=cbind(eCorrs,p.adjust(eCorrs[,2],method='fdr'))
eDEgen=cbind(eDEgen,p.adjust(eDEgen[,2],method='fdr'))
eDMcyt=cbind(eDMcyt,p.adjust(eDMcyt[,2],method='fdr'))

colnames(eCorrs)=c('FE','PValue','FDR')
colnames(eDEgen)=c('FE','PValue','FDR')
colnames(eDMcyt)=c('FE','PValue','FDR')

write.table(eCorrs,file='output/TAD.TadHypergeom.Correlations.txt',sep='\t')
write.table(eDEgen,file='output/TAD.TadHypergeom.Genes.txt',sep='\t')
write.table(eDMcyt,file='output/TAD.TadHypergeom.Cytosines.txt',sep='\t')

myComps=list(Correlations=eCorrs,DE=eDEgen,DM=eDMcyt)

pdf('output/TAD.VolcanoTadEnrich.pdf',height=4,width=4)
par(mar=c(5,5,2,2))
for (i in names(myComps)){
	cTable=myComps[[i]]
	cSign=which(cTable[,3]<0.05 & abs(log2(cTable[,1]))>1)
	
	cColors=rep('#AAAAAAAA',nrow(cTable))
	cColors[intersect(cSign,which(cTable[,1]<1))]=paste(col2hex('navy'),'AA',sep='')
	cColors[intersect(cSign,which(cTable[,1]>1))]=paste(col2hex('firebrick'),'AA',sep='')
	
	cEnrich=log2(cTable[,1])
	cFDR=-log10(cTable[,3])
	toExclude=union(which(is.infinite(cEnrich)),which(is.infinite(cFDR)))
	if (length(toExclude)>0){
		cEnrich=cEnrich[-toExclude]
		cFDR=cFDR[-toExclude]
		cColors=cColors[-toExclude]
	}
	
	plot(cEnrich,cFDR,xlab='Fold Enrichment (log2)',ylab='FDR (-log10)',col=cColors,pch=19,main=i,ylim=c(0,40),xlim=c(-6,6))
}
dev.off()




