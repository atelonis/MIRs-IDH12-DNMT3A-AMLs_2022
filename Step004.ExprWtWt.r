
library(affy)
library(hgu133plus2.db)

### focusing only on the IDH-WT and DNMT3A-WT samples
### Run this script after running the first half of Step005

wFiles=read.table('WtWtDatasets/Dataset.Colnames.txt',
				  stringsAsFactors=FALSE,sep='\t')[,1]
wIDs=gsub('Sample_','',wFiles)

allMicros=list.files( 'CEL_Files', full.names=TRUE)

sTable=read.table("sample.tsv",header=TRUE,stringsAsFactors=FALSE,sep='\t')
rownames(sTable)=gsub("AML ","",sTable[,2])

wGSM=sTable[wIDs,1]

cf=NULL
for (i in wGSM){
	cf=c(cf,grep(i,allMicros,value=TRUE))
}

wIDs=cbind(wIDs,cf)
rownames(wIDs)=gsub('CEL_Files/','',wIDs[,2])


d=ReadAffy(filenames=cf)
cProbes=probeNames(d)
drma=rma(d,normalize=FALSE)
data=exprs(drma)
cNames=select(hgu133plus2.db,cProbes,c('SYMBOL','ENTREZID','ENSEMBL','ENSEMBLPROT'))

cNames=cNames[-grep('AFFX',cNames[,1]),]
wProbes=cNames[which( apply(is.na(cNames[,c(2,3,4)]),1,sum)==0),]
gEE=paste(wProbes[,2],wProbes[,3],sep='|')
wProbes=cbind(wProbes,gEE)

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
colnames(newData)=paste('Sample_',wIDs[colnames(data),1],sep='')

write.table(newData,file='WtWtDatasets/Expression.txt',sep='\t')
write.table(wProbes,file='WtWtDatasets/MicroarryLookupIDs.txt',row.names=FALSE,quote=FALSE,sep='\t')




