
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

get_fe=function(myV,myBack){
	myNames=names(myV)
	allCombs=combn(myNames,2)
	rFE=NULL
	for (i in c(1:ncol(allCombs))){
		i1=allCombs[1,i]
		i2=allCombs[2,i]
		rFE=rbind(rFE,hypergeom.test( length(intersect(myV[[i1]],myV[[i2]])), length(myV[[i1]]), length(myV[[i2]]), myBack ))
	}
	rownames(rFE)=apply(allCombs,2,paste,collapse='-')
	return(rFE)
}

#####################################################################################
### PART 1 | Examining the overlap between DE/DM genes with genes in correlations ###
#####################################################################################

v=NULL
v$DE=rbind(read.table('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.PositSign.txt',stringsAsFactors=FALSE,header=TRUE),
		   read.table('output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.NegatSign.txt',stringsAsFactors=FALSE,header=TRUE))[,3]

for (i in c('Intermediate','Long','Overlap')){
	v[[i]]=read.table(paste('output/Correlation/Step3/TADsets.Negat_',i,'.txt',sep=''),sep='\t',stringsAsFactors=FALSE)[,1]
}

dmGenes=read.table('output/Correlation/Step3/TAD.GeneCorrsDM.txt',header=TRUE,stringsAsFactors=FALSE)

TADwide=list(DM=unique(dmGenes[which(dmGenes[,3]=="DM"),1]),DE=v$DE,Correlations=unique(c(v$Overlap,v$Intermediate,v$Long)))
venn.diagram(TADwide,fill=c('pink','palegreen','skyblue'),file='output/TAD.Integrate.DeDmCorrGenes_TADwide.png')
fe=get_fe(TADwide,length(backGenes))
print(fe)

TADwide=list(DM=unique(dmGenes[which(dmGenes[,3]=="DM"),1]),DE=v$DE)
venn.diagram(TADwide,fill=c('pink','palegreen'),file='output/TAD.Integrate.DeDmGenes_TADwide.png') 

pdf('output/TAD.Integrate.DeDmCorrGenes_TADwide.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S2A
par(mar=c(14,10,2,2))
bp=barplot(fe[,1],horiz=TRUE,names.arg=rownames(fe),xlim=c(0,4),mgp=c(1.6,0.6,0),
		   xlab='Fold Enrichment over background',las=2,xaxt='n',col='black')
axis(1,mgp=c(1.6,0.6,0))
dev.off()

for (i in c('Overlap','Intermediate','Long')){
	print(i)
	cV=list(DM=unique(dmGenes[which(dmGenes[,3]=="DM" & dmGenes[,2]==i),1]),DE=v$DE,Correlations=unique(c(v[[i]])))
	venn.diagram(cV,fill=c('pink','palegreen','skyblue'),file=paste('output/TAD.Integrate.DeDmCorrGenes_',i,'.png',sep='')) ### SUPPLEMENTARY FIGURE S2A
	fe=get_fe(cV,length(backGenes))
	print(fe)
	
	pdf(paste('output/TAD.Integrate.DeDmCorrGenes_',i,'.pdf',sep=''),height=4,width=4) ### SUPPLEMENTARY FIGURE S2A
	par(mar=c(14,10,2,2))
	bp=barplot(fe[,1],horiz=TRUE,names.arg=rownames(fe),xlim=c(0,4),mgp=c(1.6,0.6,0),
			   xlab='Fold Enrichment over background',las=2,xaxt='n',col='black')
	axis(1,mgp=c(1.6,0.6,0))
	dev.off()
}


##############################################################################################
#### PART 3 | Examining whether two genes in the same TAD are correlated with the same mCs ###
##############################################################################################

m0=read.table('output/Correlation/Step3/TADsets.GenesCsInTAD.txt',header=TRUE,stringsAsFactors=FALSE)
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

pdf('output/TAD.Enrich_GeneC_pairs.pdf',height=4,width=4) ### FIGURE 1E
par(mar=c(5,5,2,2))
plot(log2(vv),-log10(sv),xlab='Fold enrichment in common mCs (log2)',
	 ylab='FDR (-log10)',col=vCols,pch=19)
dev.off()



#####################################################################################
#### PART 4 | See how the GSEA on correlations and the GSEA on DE genes correlate ###
#####################################################################################

sam=read.table('output/GSEA/Quant3.SelectKEGG/output/output.txt',stringsAsFactors=FALSE,header=TRUE,sep='\t')
ove=read.table('output/GSEA/RankTads.Overlap.SelectKEGG/output/output.txt',stringsAsFactors=FALSE,header=TRUE,sep='\t')

sam[,2]=-sam[,2]
rownames(sam)=sam[,1]
rownames(ove)=ove[,1]
coms=intersect(ove[,1],sam[,1])

m=cbind(sam[coms,2],ove[coms,2])
f=cbind(sam[coms,3],ove[coms,3])
rownames(m)=coms
rownames(f)=coms
f=f<0.1
f[which(f)]='DD'
f[which(f==FALSE)]='40'
fCols=paste('#',f[,1],f[,2],'40DD',sep='')

mLim=max(abs(m))
mLim=c(-mLim,mLim)
pdf('output/Compare_DE_IQA/TAD.Points.selectKEGG.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S2C
par(mar=c(5,5,2,2))

plot(m[,1],m[,2],xlim=mLim,ylim=mLim,xlab='Enrichment in DE',ylab='Enrichment in Correlations',
	 main='KEGG',pch=19,col=fCols)

fp=NULL
for (i in c('03010_RIBOSOME',
			'00010_GLYCOLYSIS_/_GLUCONEOGENESIS',
			'04330_NOTCH_SIGNALING_PATHWAY',
			'00051_FRUCTOSE_AND_MANNOSE_METABOLISM',
			'00561_GLYCEROLIPID_METABOLISM',
			'04310_WNT_SIGNALING_PATHWAY',
			'04012_ERBB_SIGNALING_PATHWAY',
			'04340_HEDGEHOG_SIGNALING_PATHWAY',
			'03060_PROTEIN_EXPORT')){
	iRow=which(coms==i)
	points(m[iRow,1],m[iRow,2],pch=21,col='blue',bg=fCols[iRow])
	fp=rbind(fp,(c(i,m[iRow,])))
}
dev.off()

