
library(samr)
library(affy)
library(u133aaofav2cdf)

######################################
###  PART 1  |  process CEL files  ###
######################################

cf=list.files('Hematopoiesis/GSE24759_RAW/',full.names=TRUE,pattern='GSM')

d=ReadAffy(filenames=cf)
cProbes=probeNames(d)
drma=rma(d,normalize=FALSE)
data=exprs(drma)
write.table(data,'RMA_AllData.txt',sep='\t')

data=read.table('RMA_AllData.txt',sep='\t')
lookup=read.table('LookupTable.txt',stringsAsFactors=FALSE,sep='\t')
myGenes=sort(unique(lookup[,2]))

d=NULL
for (i in myGenes){
	cProbes=lookup[lookup[,2]==i,1]
	if (length(cProbes)>1){
		nd=apply(data[cProbes,],2,mean)
	} else {
		nd=data[cProbes,]
	}
	d=rbind(d,nd)
}

rownames(d)=myGenes
colnames(d)=colnames(data)

myNames=read.table('Hematopoiesis/GSM_Accessions.edited.txt',sep='\t',stringsAsFactors=FALSE)
rownames(myNames)=myNames[,1]

newCols=NULL
for (i in colnames(d)){
	cGSM=strsplit(i,'_')[[1]][1]
	newCols=c(newCols,gsub(',','',myNames[cGSM,2]))
}

colnames(d)=newCols

### Select most expressed genes
wExp=which(rowMeans(d)>quantile(rowMeans(d),0.5))
write.table(d[wExp,],file='Hematopoiesis/Dataset.50percMostExpr.txt',sep='\t')

### Select most expressed genes
wExp=which(rowMeans(d)>quantile(rowMeans(d),0.25))
write.table(d[wExp,],file='Hematopoiesis/Dataset.75percMostExpr.txt',sep='\t')

### Select the most variable genes
wVAR=which(apply(d,1,var)>quantile(apply(d,1,var),0.5))
write.table(d[wVAR,],file='Hematopoiesis/Dataset.50percMostVar.txt',sep='\t')




##################################################
### PART 2  | Differential expression analyses ###
##################################################

Run_SAM=function(curda,y,myoutput){
	
	# exlude features with zero SD within at least one of the groups
	toInclude=unique( names(which(apply(curda[,which(y==1)],1,sd)>0)), 
					  names(which(apply(curda[,which(y==2)],1,sd)>0)) )
	curda=curda[toInclude,]
	write.table(rownames(curda),sep='\t',file=paste(myoutput,'.SAM.IncludedInSAM.txt',sep=''),
				row.names=FALSE,col.names=FALSE,quote=FALSE)
	
	# Run R
	mylist=list(x=curda,y=y,geneid=rownames(curda),logged2=FALSE)
	mySAM=samr(mylist,resp.type="Two class unpaired",nperm=5000)
	
	# Calculate Delta
	delta.table = samr.compute.delta.table(mySAM)
	delta = -1
	for(i in 1:nrow(delta.table)){
		if(!((delta.table[i,5]>0.05) | is.na(delta.table[i,5]))){
			delta = delta.table[i,1];d1=delta.table[i,]
			break()}
		}
	if (delta==-1){
		delta=delta.table[i,1]
		}
	
	# Save observed score
	write.table(mySAM$tt,file=paste(myoutput,'.SAM.Score.txt',sep=''),col.names=FALSE,quote=FALSE,sep='\t')
	
	# Save R graph
	write.table(delta,file=paste(myoutput,'.SAM.DeltaValue.txt',sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE)
	pdf(paste(myoutput,'.SAM.Graph.delta.pdf',sep=''),height=5,width=5)
	par(mar=c(5,5,2,2))
	samr.plot(mySAM,delta)
	dev.off()
	
	# Calculate and save significant features
	siggenes.table = samr.compute.siggenes.table(mySAM,delta,mylist,delta.table)
	write.table(siggenes.table$genes.up,sep='\t',file=paste(myoutput,'.SAM.PositSign.txt',sep=''),row.names=FALSE)
	write.table(siggenes.table$genes.lo,sep='\t',file=paste(myoutput,'.SAM.NegatSign.txt',sep=''),row.names=FALSE)
}


dir.create('Hematopoiesis/output',showWarnings=FALSE)
dir.create('Hematopoiesis/output/SAM',showWarnings=FALSE)

d=read.table('Hematopoiesis/Dataset.75percMostExpr.txt')

myCells=NULL
for (i in colnames(d)){
	myCells=c(myCells,strsplit(i,'_')[[1]][1])
}
uCells=sort(unique(myCells))

uCols=rainbow(length(uCells))
names(uCols)=uCells

hsc=which(myCells=='HSC')
for (i in uCells){
	if (i=='HSC'){
		next
	}
	cCells=which(myCells==i)
	cD=zd[,c(hsc,cCells)]
	y=c(rep(1,length(hsc)),rep(2,length(cCells)))
	Run_SAM(as.matrix(cD),y,paste('Hematopoiesis/output/SAM/ZScored_RMA_',i,'_vs_HSC',sep=''))
}




