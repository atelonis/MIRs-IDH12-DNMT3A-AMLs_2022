
library(circlize)
library(data.table)
library(VennDiagram)

# Set colors
hoxCol='#FF4635'
snxCol='#1818A5'

# Coordinates of cytocines 
backc=read.table('output/ForCircos/background.txt',stringsAsFactors=FALSE,sep='\t')
hoxb5=read.table('output/ForCircos/hoxb5.txt',stringsAsFactors=FALSE,sep='\t')
snx11=read.table('output/ForCircos/snx11.txt',stringsAsFactors=FALSE,sep='\t')

# hoxb tad on chr 17
hTad=c(45880000,48560000)

# Coordinates of genes
geneCoords=read.table('HumanGenome/GeneRanges.bed',stringsAsFactors=FALSE,sep='\t')
backGenes=read.table('output/Correlation/Step3/TADsets.Background.ENSG.txt',stringsAsFactors=FALSE)[,1]

# Keep genes and cytosines within the TAD
wGenes=geneCoords[which(geneCoords[,1]=="chr17" & geneCoords[,2]>hTad[1] & geneCoords[,3]<hTad[2]),]
wGenes=wGenes[which(wGenes[,4] %in% backGenes),]
wCytos=backc[which(backc[,1]=='chr17' & backc[,2]>hTad[1] & backc[,2]<hTad[2]),]

# Start plotting
pdf('output/Correlation/HOXB_TAD.pdf')  ### FIGURE 1D
circos.par( 'track.height'=0.2,gap.degree=5)#, start.degree=85, gap.degree=10)
circos.initialize( as.factor(c("Genes","Cytos")), xlim=c(hTad[1],hTad[2]))

cLabs=pretty(hTad,3)
circos.track(as.factor(c("Genes")),ylim=c(0,1))
circos.axis(labels.cex=1,major.at=cLabs,labels=format(cLabs,scientific=FALSE))
circos.axis(labels.cex=1,major.at=-cLabs+sum(hTad),labels=format(cLabs,scientific=FALSE),sector.index="Cytos")

# Plot the genes on the "top" sector
for (i in c(1:nrow(wGenes))){
	if (wGenes[i,4] %in% c('ENSG00000120075|HOXB5','ENSG00000002919|SNX11')){
		next
	}
	circos.rect( wGenes[i,2], 0, wGenes[i,3], 1, col="#AAAAAAAA", border=FALSE,sector.index="Genes")
}
hoxbGene=wGenes[which(wGenes[,4]=='ENSG00000120075|HOXB5'),]
circos.rect(hoxbGene[2],0,hoxbGene[3],1,col=hoxCol,border=FALSE,sector.index="Genes")
snxGene=wGenes[which(wGenes[,4]=='ENSG00000002919|SNX11'),]
circos.rect(snxGene[2],0,snxGene[3],1,col=snxCol,border=FALSE,sector.index="Genes")

# Plot the cytosines on the "bottom" sector"
circos.points(-wCytos[,2]+sum(hTad),jitter(rep(0.5,nrow(wCytos)),a=0.5),pch=1,cex=0.5,col='#00000050',sector.index="Cytos")

# Link the genes with the cytosines with which they are correlated
for (i in c(1:nrow(hoxb5))){
	circos.link("Genes", mean(as.numeric(hoxbGene[2:3])), "Cytos", -hoxb5[i,2]+sum(hTad),
				col=paste(hoxCol,'AA',sep=''),border=FALSE )
}
for (i in c(1:nrow(snx11))){
	circos.link("Genes", mean(as.numeric(snxGene[2:3])), "Cytos", -snx11[i,2]+sum(hTad),
				col=paste(snxCol,'AA',sep=''),border=FALSE )
}

# Add title and finalize
title(paste('HOXB TAD | chr17:',hTad[1],'-',hTad[2],sep=''))
circos.clear()
dev.off()


### Venn diagram
v=list(HOXB5=as.character(hoxb5[,2]),SNX11=as.character(snx11[,2]))
venn.diagram(v,fill=c(hoxCol,snxCol),file='output/TAD.HOXB5_SNX11.png')


### Remove log files
for (i in list.files('output/',full.names=TRUE,pattern='log')){
	ii=strsplit(i,'[.]')[[1]]
	if (ii[length(ii)]=='log'){
		file.remove(i)
	}
}

