
library(Hmisc)
source('~/Scripts/Heatmap3.r')

myGSEAs=list.files('TCGA/output/GSEA',pattern='Select.RankTads')
vE=NULL
vF=NULL
allPaths=NULL
for (i in myGSEAs){
	cPaths=read.table(paste('TCGA/output/GSEA',i,'TCGA/output/TCGA/output.txt',sep='/'),
					  sep='\t',header=TRUE,stringsAsFactors=FALSE)
	toExclude=unique(which(is.na(cPaths),arr.ind=T)[,1])
	if (length(toExclude)>0){
		cPaths=cPaths[-toExclude,]
	}
	Keggs=cPaths[which(cPaths[,3]<0.2),]
	allPaths=c(allPaths,cPaths[,1])
	
	# Enrichments
	vE[[strsplit(i,'[.]')[[1]][3]]]=cPaths[,1:2]
	rownames(vE[[strsplit(i,'[.]')[[1]][3]]])=cPaths[,1]
	
	# FDRs
	vF[[strsplit(i,'[.]')[[1]][3]]]=cPaths[,c(1,3)]
	rownames(vF[[strsplit(i,'[.]')[[1]][3]]])=cPaths[,1]
}

m=matrix(0,nrow=length(unique(allPaths)),ncol=length(names(vE)))
rownames(m)=unique(allPaths)
colnames(m)=names(vE)
for (i in names(vE)){
	m[rownames(vE[[i]]),i] = as.numeric(vE[[i]][,2])
	m[rownames(vF[[i]])[which(vF[[i]][,2]>0.1)],i]=0
}
m=m[-which(rowSums(abs(m))==0),]

pathNames=NULL
for (i in read.table('Scripts/GMT.selKEGG_Ensembl.gmt',sep=' ',stringsAsFactors=F)[,1]){
	pathNames=c(pathNames,strsplit(i,'\t')[[1]][1])
}
names(pathNames)=substr(pathNames,1,5)
pathNames=gsub('_',' ',gsub('[0-9][0-9][0-9][0-9][0-9]_','',pathNames))

rownames(m)=pathNames[substr(rownames(m),1,5)]

m=m[,c('Overlap','Intermediate','Long','Whole')]
colnames(m)=c('Proximal','Intermediate','Long','TAD-wide')
m=m[order(m[,1],m[,2],m[,3],m[,4],decreasing=TRUE),]
fm=m[-unique(which(m>0 & m<1.9,arr.ind=T)[,1]),]
hc=colorRampPalette(c('navy','gray90','firebrick'))(21)

pdf('TCGA/output/Heatmap.SelectGseaTads.pdf')  ### FIGURE 1G
heatmap.3(m,col=hc,margins=c(12,28),breaks=seq(-2.5,2.5,length.out=length(hc)+1),
		  dendrogram='n',Colv='n',Rowv='n',KeyValueName='Enrichment score')
dev.off()


