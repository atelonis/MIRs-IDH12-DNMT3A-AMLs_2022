
library(amap)
library(scico)
library(alluvial)
library(circlize)
library(dendextend)

source('~/Scripts/Heatmap3.r')

geneIDs=read.table('HumanGenome/GeneRanges.bed',sep='\t',stringsAsFactors=FALSE)
geneNames=matrix(unlist(strsplit(geneIDs[,4],'[|]')),byrow=TRUE,ncol=2)[,2]

toExclude=which(duplicated(geneNames))
geneIDs=geneIDs[-toExclude,]
geneNames=geneNames[-toExclude]
rownames(geneIDs)=geneNames

d0=read.table('output/TAD.EnrichPPIs.txt',stringsAsFactors=FALSE,sep='\t')
d0[,2]=gsub('Overlap','Proximal',d0[,2])
d0=cbind(d0,p.adjust(as.numeric(d0[,8]),method='fdr'))

colnames(d0)=c('Pair','Annotation','k','s','M','N','FoldEnrich','PValue','FDR')
write.table(d0,'output/TAD.EnrichPPIs.FDR.txt',sep='\t',quote=FALSE,row.names=FALSE)

n=d0[which(as.numeric(d0[,9])<0.05 & log2(d0[,7])<(-2)),]
d=d0[which(as.numeric(d0[,9])<0.05 & log2(d0[,7])>2),]
ppiPairs=sort(unique(d[,1]))


m=matrix(0,nrow=length(ppiPairs),ncol=4)
rownames(m)=ppiPairs
colnames(m)=c('Proximal','Intermediate','Long','TADwide')

for (i in c(1:nrow(d))){
	m[ d[i,1], d[i,2] ] = log2(as.numeric(d[i,7]))
}

m[is.infinite(m)]=-5

### Write files

wb=sort(unique(unlist(strsplit(d0[which(d0[,2]=='TADwide'),1],'--'))))
wb=geneIDs[wb,4]
write.table(wb,'output/PPIs.Background.ENSG.txt',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

wd=unique(unlist(strsplit(d[which(d[,2]=='TADwide'),1],'--')))
wd=geneIDs[wd,4]
write.table(wd,'output/PPIs.Enriched.ENSG.txt',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

wn=unique(unlist(strsplit(n[which(n[2]=='TADwide'),1],'--')))
wn=geneIDs[wn,4]
write.table(wn,'output/PPIs.Depleted.ENSG.txt',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


### Enriched
pdf('output/TAD.EnrichPPIs_Enriched.pdf',height=4,width=4)  ### FIGURE 3F
td=d[which(d[,2]=='TADwide'),]
td=td[order(td[,1]),]
td1=matrix(unlist(strsplit(td[,1],'--')),ncol=2,byrow=TRUE)[,1]
td2=matrix(unlist(strsplit(td[,1],'--')),ncol=2,byrow=TRUE)[,2]
dCols=colorRampPalette(c('black','red','orange','gold'))(length(unique(td1)))
dCols=scico(length(unique(td1)),palette=c('batlow'))
names(dCols)=unique(td1)
dCols=dCols[td1]
alluvial(td1,td2,freq=rep(1,nrow(td)),axis_labels=c('mC','Gene'),col=dCols,alpha=0.8,cex=1)
dev.off()


### Depleted
tt1=table(tn1)
tt2=table(tn2)
s12=c(names(tt1),paste('2__',names(tt2),sep=''))
sw=c(tt1,tt2)
names(sw)=s12

pdf('output/TAD.EnrichPPIs_Depleted.Circos.pdf',height=5,width=5) ### FIGURE 3G
circos.par('track.height'=0.2) #c( 5, rep(5,length(tt1)), 10, rep(5,length(tt2)) ))
circos.initialize(s12, xlim=c(0,1), sector.width=sw)
circos.track(s12,ylim=c(0,1))
for (i in names(tt2)){
	circos.text(0.5,0.5,gsub('2__','',i),sector.index=paste('2__',i,sep=''),niceFacing=TRUE,facing='clockwise',cex=0.5)
}
for (i in names(tt1)){
	circos.text(0.5,0.5,gsub('2__','',i),sector.index=i,niceFacing=TRUE,facing='clockwise',cex=0.5)
}
for (i in c(1:length(tn1))){
	# First
	i1=tn1[i]
	if (tt1[i1]==1){
		lim1=c(0,1)
	} else {
		noBins=1/tt1[i1]
		yetToPlot=length(which(tn1[(i+1):length(tn1)]==i1))
		lim1=c(noBins*yetToPlot,noBins*(yetToPlot+1))
	}
	# Second
	i2=tn2[i]
	if (tt2[i2]==1){
		lim2=c(0,1)
	} else {
		noBins=1/tt2[i2]
		alreadyPlot=length(which(tn2[1:(i-1)]==i2))
		yetToPlot=length(which(tn2[(i+1):length(tn2)]==i2))
		lim2=c(noBins*yetToPlot,noBins*(yetToPlot+1))
	}
	circos.link(tn1[i],lim1,paste('2__',tn2[i],sep=''),lim2,col=paste(nCols[tn1[i]],'AA',sep=''))
}
dev.off()


