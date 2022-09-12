
library(samr)
library(DiscriMiner)
library(VennDiagram)

source('Scripts/Heatmap3.r')

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

dir.create('output',showWarnings=FALSE)
dir.create('output/SAM',showWarnings=FALSE)

d=read.table('Microarrays/Dataset.50percMostExpr.txt')
colnames(d)=gsub('IDH_DNMT','Double',colnames(d))

idh=grep('IDH',colnames(d),value=TRUE)
dnm=grep('DNMT3A',colnames(d),value=TRUE)
dou=grep('Double',colnames(d),value=TRUE)

# Quantile normalization
dd=d[,c(idh,dnm,dou)]

dr=apply(dd,2,rank)
ds=apply(dd,2,sort)
dsm=rowMeans(ds)

qd=matrix(0,nrow=nrow(dd),ncol=ncol(dd))
rownames(qd)=rownames(dd)
colnames(qd)=colnames(dd)
for (i in c(1:length(dsm))){
	qd[dr==i]=dsm[i]
}

y=c(rep(2,length(idh)),rep(1,length(dnm)))
Run_SAM(as.matrix(qd[,c(idh,dnm)]),y,'output/SAM/Quantile_3.IDH_vs_DNMT3A')


