
library(amap)
library(gplots)
library(data.table)
library(VennDiagram)
source('~/Scripts/Heatmap3.r')

hypergeom.test=function(k,s,M,N){
	# k: number of successes
	# s: sample size
	# M: Number of successes in the population
	# N: population size
	myFE = (k/s) / (M/N)
	if (is.na(myFE)){
		p=NA
	} else if (myFE>1){
		p=1-phyper( k-1, M, N-M, s )
	} else {
		p=phyper( k, M, N-M, s )
	}
	return(c(myFE,p))
}

##############################################################
### PART 1 | Extract the enrichments/depletions of each TF ###
##############################################################

dir.create('output/Cs_Peaks',showWarnings=FALSE)

bC=read.table('output/Correlation/Step3/Cs.Background.txt' ,sep='\t',stringsAsFactors=FALSE)[,1] 
sC=read.table('output/Correlation/Step3/Cs.Significant.txt',sep='\t',stringsAsFactors=FALSE)[,1] 


### Beck
beck=list.files('Beck',pattern='BI.Cs.')
beck=matrix(unlist(strsplit(beck,'[.]')),ncol=4,byrow=TRUE)[,3]

beckHT=NULL
for (i in beck){
	print(i)
	
	cBeck=read.table(paste('Beck/BI.Cs',i,'bed',sep='.'),sep='\t',stringsAsFactors=FALSE)
	cBeck_Cs=paste(unique(gsub(' ','',apply(cBeck[,c(1,2)],1,paste,collapse='.'))),'.NA0',sep='')
	cBeck_peaks=unique(gsub(' ','',apply(cBeck[,4:6],1,paste,collapse='.')))
	
	cHT=hypergeom.test( length(intersect(sC,cBeck_Cs)), length(sC), length(intersect(bC,cBeck_Cs)), length(bC))
	
	beckHT=rbind(beckHT,c(i,cHT))
}
beckHT=cbind(beckHT,p.adjust(as.numeric(beckHT[,3]),method='fdr'))
colnames(beckHT)=c('TF','FE','PValue','FDR')

write.table(beckHT,'output/Cs_Peaks/Beck.txt',sep='\t',quote=FALSE,row.names=FALSE)


### Encode
encode=list.files('Encode',pattern='BI.Cs.')
encode=matrix(unlist(strsplit(encode,'[.]')),ncol=4,byrow=TRUE)[,3]

encodeHT=NULL
for (i in encode){
	print(i)
	
	cEncode=try(read.table(paste('Encode/BI.Cs',i,'bed',sep='.'),sep='\t',stringsAsFactors=FALSE),silent=TRUE)
	if (is(cEncode)[1]=='try-error'){
		next
	}
	cEncode_Cs=paste(unique(gsub(' ','',apply(cEncode[,c(1,2)],1,paste,collapse='.'))),'.NA0',sep='')
	cEncode_peaks=unique(gsub(' ','',apply(cEncode[,4:6],1,paste,collapse='.')))
	
	cHT=hypergeom.test( length(intersect(sC,cEncode_Cs)), length(sC), length(intersect(bC,cEncode_Cs)), length(bC))
	
	encodeHT=rbind(encodeHT,c(i,cHT))
}
encodeHT=encodeHT[-which(is.na(as.numeric(encodeHT[,2]))),]
encodeHT=encodeHT[-which(as.numeric(encodeHT[,2])==0),]
encodeHT=cbind(encodeHT,p.adjust(as.numeric(encodeHT[,3]),method='fdr'))
colnames(encodeHT)=c('DBP','FE','PValue','FDR')

write.table(encodeHT,'output/Cs_Peaks/Encode.txt',sep='\t',quote=FALSE,row.names=FALSE)


#################################
### PART 2 | Plot as barplots ###
#################################


myFiles=gsub('.txt','',list.files('output/Cs_Peaks',pattern='txt'))

pdfMars=c(44,43,2)
names(pdfMars)=c('Beck','Encode')

xlims=list( Beck=c(-3,1), Encode=c(-6,2) )

for (i in myFiles){
		cE=read.table(paste('output/Cs_Peaks/',i,'.txt',sep=''),stringsAsFactors=FALSE,sep='\t',header=TRUE)
		cE[,2]=log2(cE[,2])
		cE=cE[which(!is.infinite(cE[,2])),]
		cE=cE[order(cE[,2]),]
		
		sStar=rep('',nrow(cE))
		wStar=which(cE[,'FDR']<0.05 & abs(cE[,2])>1)
		if (length(wStar)>0){
			sStar[wStar]='*'
		}
		
		cColors=rep('gray',nrow(cE))
		cColors[sStar=='*']='navy'
		
		pdf(paste('output/Cs_Peaks/',i,'.pdf',sep=''),height=10,width=4) ### FIGURE 3D 3E
		par(mar=c(pdfMars[i],7,2,2))
		if (i=='Encode'){
			cNames=cE[,1]
			cNames[sStar=='']=''
			bp=barplot(cE[,2],main=i,horiz=TRUE,las=2,xaxt='n',xlab='Fold Enrichment (log2)',col=cColors,xlim=xlims[[i]],
					   names.arg=cNames,border=NA,cex.names=0.2)
		} else {
			bp=barplot(cE[,2],main=i,horiz=TRUE,las=2,xaxt='n',xlab='Fold Enrichment (log2)',col=cColors,xlim=xlims[[i]],
					   names.arg=cE[,1],border=NA)
		}
		axis(1)
		dev.off()
}



