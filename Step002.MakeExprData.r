
library(affy)
library(hgu133plus2.db)

dir.create('Microarrays',showWarnings=FALSE)

### Get the samples that also have methylation data
wFiles=read.table('Datasets/Dataset.Colnames.txt',
				  stringsAsFactors=FALSE,sep='\t')[,1]

wIDs=NULL
for (i in wFiles){
	i=strsplit(i,'_')[[1]]
	if (length(i)==3){
		i=c('IDH_DNMT3A',i[3])
	}
	wIDs=rbind(wIDs,c(i[2],paste(i[1],i[2],sep='_')))
}

### Downloaded from GEO GSE6891
allMicros=list.files( 'CEL_Files', full.names=TRUE)
sTable=read.table("sample.tsv",header=TRUE,stringsAsFactors=FALSE,sep='\t')
rownames(sTable)=gsub("AML ","",sTable[,2])

wGSM=sTable[wIDs[,1],1]

cf=NULL
for (i in wGSM){
	cf=c(cf,grep(i,allMicros,value=TRUE))
}

wIDs=cbind(wIDs,cf)
rownames(wIDs)=gsub('CEL_Files/','',wIDs[,3])

### Read and process microarray data
d=ReadAffy(filenames=cf)
cProbes=probeNames(d)
drma=rma(d,normalize=FALSE)
data=exprs(drma)
cNames=select(hgu133plus2.db,cProbes,c('SYMBOL','ENTREZID','ENSEMBL','ENSEMBLPROT'))

cNames=cNames[-grep('AFFX',cNames[,1]),]
wProbes=cNames[which( apply(is.na(cNames[,c(2,3,4)]),1,sum)==0),]
gEE=paste(wProbes[,2],wProbes[,3],sep='|')
wProbes=cbind(wProbes,gEE)
write.table(wProbes,file='Microarrays/LookupIDs.txt',row.names=FALSE,quote=FALSE,sep='\t')

### Build matrix
newData=NULL
newRows=sort(unique(gEE))
for (i in newRows){
	i1=strsplit(i,'[|]')[[1]][1]
	i2=strsplit(i,'[|]')[[1]][2]
	geneProbes=sort(unique(wProbes[intersect(which(wProbes[,2]==i1),which(wProbes[,3]==i2)),1]))
	if (length(geneProbes)>1){
		nd=apply(data[geneProbes,],2,mean)
	} else {
		nd=data[geneProbes,]
	}
	newData=rbind(newData,nd)
}
rownames(newData)=newRows
colnames(newData)=wIDs[colnames(data),2]

### Select most expressed genes
wExp=which(rowMeans(newData)>quantile(rowMeans(newData),0.5))
write.table(newData[wExp,],file='Microarrays/Dataset.50percMostExpr.txt',sep='\t')


