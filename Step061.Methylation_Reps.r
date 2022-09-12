
library(samr)
library(data.table)

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


d=as.matrix(fread('Datasets/Dataset.Data.txt'))
rownames(d)=read.table('Datasets/Dataset.Rownames.txt',stringsAsFactors=FALSE,sep='\t')[,1]
colnames(d)=read.table('Datasets/Dataset.Colnames.txt',stringsAsFactors=FALSE,sep='\t')[,1]

repC=read.table('bed/BI.Cs.Reps.bed',stringsAsFactors=FALSE,sep='\t')

for (am in c('MIR','Alu','L1','ERVL')){
	print(am)
	
	if (am %in% c('MIR','Alu')){
		next
		cC=repC[grep(paste('SINE',am,sep='.'),repC[,7]),]
	} else if (am=='L1') {
		next
		cC=repC[grep(paste('LINE',am,sep='.'),repC[,7]),]
	} else { # ERVL
		cC=repC[grep('LTR/ERVL[|]',repC[,7]),]
	}
	
	cC=paste(cC[,1],cC[,3],'NA0',sep='.')
	cC=intersect(rownames(d),cC)
	
	dnmt = grep('IDH',grep('DNMT3A',colnames(d),value=TRUE),value=TRUE,invert=TRUE)
	idh  = grep('DNMT3A',grep('IDH',colnames(d),value=TRUE),value=TRUE,invert=TRUE)
	
	y=c(rep(1,length(dnmt)),rep(2,length(idh)))
	Run_SAM(as.matrix(d[cC,c(dnmt,idh)]),y,paste('output/SAM/Methylation',am,'IDH_vs_DNMT',sep='.'))
}


