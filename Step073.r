
library(samr)
library(amap)
library(gplots)
library(alluvial)
library(dendextend)
library(DiscriMiner)
library(VennDiagram)

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


dir.create('GSE60055/output',showWarnings=FALSE)
dir.create('GSE60055/output/SAM',showWarnings=FALSE)

d=read.table('GSE60055/ExprMatrix.FPKM_Filter.txt',stringsAsFactors=FALSE,sep='\t',check.names=FALSE)
d=log2(d+1)

for (mg in c('MEP')){
	# Add an "average" sample
	d=cbind(d,rowMeans(d[,grep(paste('KO3a_GFP_._',mg,sep=''),colnames(d))]))
	d=cbind(d,rowMeans(d[,grep(paste('WT3a_M140_._',mg,sep=''),colnames(d))]))
	colnames(d)=c(names(myColors),paste('KO3a_GFP_A_',mg,sep=''),paste('WT3a_M140_A_',mg,sep=''))
	
	cD=d[,grep(mg,colnames(d))]
	cD=cbind(cD,d[,grep('KO3a_M140',colnames(d))])
	
	wt  = grep('WT3a_GFP',colnames(cD),value=TRUE)
	idh = grep('WT3a_M140',colnames(cD),value=TRUE)
	dnm = grep('KO3a_GFP',colnames(cD),value=TRUE)
	dou = grep('KO3a_M140',colnames(cD),value=TRUE)
	
	# IDH vs WT
	y=c(rep(1,length(wt)),rep(2,length(idh)))
	Run_SAM(as.matrix(cD[,c(wt,idh)]),y,paste('GSE60055/output/SAM/',mg,'.M140_vs_WT',sep=''))
	
	# DNMT vs WT
	y=c(rep(1,length(wt)),rep(2,length(dnm)))
	Run_SAM(as.matrix(cD[,c(wt,dnm)]),y,paste('GSE60055/output/SAM/',mg,'.KO3a_vs_WT',sep=''))
	
	# Double vs WT
	y=c(rep(1,length(wt)),rep(2,length(dou)))
	Run_SAM(as.matrix(cD[,c(wt,dou)]),y,paste('GSE60055/output/SAM/',mg,'.Double_vs_WT',sep=''))
}


