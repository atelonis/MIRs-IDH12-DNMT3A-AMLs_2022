
reps=c('MIR','Alu','L1')

### Pie charts for number of DM cytosines
DM=NULL
BI=NULL
BA=NULL
pdf('output/SAM/Pies.NumberOfDM.pdf') ### FIGURE 4A
for (am in reps){
	cIncl=read.table(paste('output/SAM/Methylation',am,'IDH_vs_DNMT.SAM.IncludedInSAM.txt',sep='.'),stringsAsFactors=FALSE)[,1]
	cPos=read.table(paste('output/SAM/Methylation',am,'IDH_vs_DNMT.SAM.PositSign.txt',sep='.'),stringsAsFactors=FALSE,header=TRUE)[,3]
	pie(c(length(cIncl)-length(cPos),length(cPos)),col=c('gray','firebrick'),labels=c('Non-DM','DM'),main=am)
	
	DM[[am]]=gsub('.NA0','',cPos)
	BI[[am]]=read.table(paste('output/GEAR/SAM_',am,'.BedInters.C.GeneExt.bed',sep=''),stringsAsFactors=FALSE,sep='\t')
	BI[[am]]=cbind( paste(BI[[am]][,1],BI[[am]][,3],sep='.'), BI[[am]][,7] )
	BA[[am]]=gsub('.NA0','',cIncl)
}
dev.off()


### Pie charts for protein-coding/intergenic differences
for (am in reps){
	cTable=read.table(paste('output/GEAR/SAM_',am,'.Genes.txt',sep=''),sep='\t',stringsAsFactors=FALSE)
	pdf(paste('output/GEAR/Pie_SAM',am,'pdf',sep='.')) ### FIGURE 4B 4C 4D
	for (db in c(2,3)){
		p=NULL
		for (i in c(2,3,4)){
			p=c(p, as.numeric(strsplit(cTable[db,i],' ')[[1]][1]))
		}
		pie(p,labels=c(cTable[1,c(2,3,4)]),col=c('tan1','palegreen','plum'),main=cTable[db,1])
		#print(c(am,round(100*p/sum(p))))
	}
	dev.off()
}


### Barplots for enrichment in enhancers
enhancer=read.table('output/Enhancers/Table.txt')
EnhEnr=NULL
for (am in reps){
	EnhEnr=c(EnhEnr, (enhancer[am,4]/enhancer[am,3]) / (enhancer[am,2]/enhancer[am,1]) )
}
names(EnhEnr)=reps
pdf('output/Enhancers/Barplots.pdf',height=5,width=5) ### FIGURE 4F
par(mar=c(20,12,2,2))
barplot(log2(rev(EnhEnr)),col='firebrick',horiz=TRUE,las=2,main='Enhancers',
		xlab='Fold enrichment (log2)',xaxt='n',border=NA,xlim=c(0,3))
axis(1)
dev.off()


### Barplots for enrichment in IQA and BTM1000 genes

b=NULL
p=NULL
for (ib in c('TADsets.Negat_Overlap','TADsets_W.Overlap')){
	ibGenes=read.table(paste('output/Correlation/Step3/',ib,'.ENSG.txt',sep=''),stringsAsFactors=FALSE,sep='\t')[,1]
	for (am in c('MIR','Alu','L1')){
		cIB = unique(BI[[am]][ BI[[am]][,2] %in% ibGenes, 1])
		cDM = DM[[am]]
		cBA = BA[[am]]
		
		noSucc = length(which(cIB %in% cDM))
		noSize = length(cDM)
		inBack = length(cIB)
		noBack = length(cBA)
		b=c(b, (noSucc/noSize) / (inBack/noBack) )
		if (b[length(b)]>1){
			p=c(p,1-phyper( noSucc-1, inBack, noBack-inBack, noSize ))
		} else {
			p=c(p,phyper( noSucc, inBack, noBack-inBack, noSize ))
		}
	}
}
b=log2(rev(b))
p=rev(p)

colB=rep('gray',length(b))
colB[p<0.0102 & b>0]='firebrick'
colB[p<0.0102 & b<0]='navy'

pdf('output/GEAR/TAD.Zoomed.pdf',height=5,width=5)  ### FIGURE 4E
par(mar=c(15,12,2,2))
bp=barplot(b,col=colB,border=NA,horiz=TRUE,las=2,xlim=c(-3,3),
		names.arg=rep(c('L1','Alu','MIR'),2),
		xaxt='n',xlab='Fold enrichment (log2)')
axis(1)
s=rep('',length(b))
s[which(colB!='gray')]='*'
text(b+(sign(b)*0.25),bp,s,col=colB,cex=2)
dev.off()






