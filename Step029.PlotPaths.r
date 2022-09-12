
library(Hmisc)
source('Scripts/Heatmap3.r')

###### DAVID
myDAVIDs=rev(grep('DE.TAD',list.files('output/DAVID',pattern='TADsets'),invert=TRUE,value=TRUE))

goDE=NULL
for (i in myDAVIDs){
	if (i=='TADsets.Negat_Whole'){
		next
	}
	cDE=read.table(paste('output/DAVID/DE.',i,'/Output.FDRfilter.txt',sep=''),
				   sep='\t',header=TRUE,stringsAsFactors=FALSE)
	goDE[[strsplit(i,'[.]')[[1]][2]]]=cDE[which(cDE[,1]=='GOTERM_BP_FAT'),c('Term','Fold.Enrichment')]
}

pdf('output/Barplot.DAVID_DeBackgr.pdf',height=4,width=5)
par(mar=c(4,15,2,2))
for (i in names(goDE)){
	print(i)
	myBars=goDE[[i]]
	myBars=myBars[order(myBars[,2]),]
	
	myCols=rep('black',nrow(myBars))
	myCols[which(myBars[,1] %in% go[[i]][,1])]='goldenrod'
	
	myBars[,1]=matrix(unlist(strsplit(myBars[,1],'~')),ncol=2,byrow=TRUE)[,2]
	myBars[,1]=paste(toupper(substr(myBars[,1],1,1)),substr(myBars[,1],2,1000),sep='')
	
	barplot(myBars[,2],col=myCols,border=NA,horiz=TRUE,xlab='Fold enrichment',las=2,xaxt='n',mgp=c(1.6,0.6,0),
			names.arg=myBars[,1],main=gsub('Negat_','',gsub('Overlap','Proximal',i)),cex.names=0.4)
	axis(1,mgp=c(1.6,0.6,0))
}
dev.off()


##############
###  GSEA  ###
##############

KEGGs=NULL ### This is to rename the pathways with lower-cases
for (i in list.files('HumanGenome/',pattern='PathwayLists')){
	cKegg=read.table(paste('HumanGenome/',i,sep=''),sep='\t',
					 stringsAsFactors=FALSE,colClasses="character")
	KEGGs=rbind(KEGGs,cKegg)
}
rownames(KEGGs)=KEGGs[,2]

### Read the GSEA output files
myGSEAs=list.files('output/GSEA',pattern='RankTads')
v=NULL
for (i in myGSEAs){
	cPaths=read.table(paste('output/GSEA',i,'output/output.txt',sep='/'),
					  sep='\t',header=TRUE,stringsAsFactors=FALSE)
	toExclude=unique(which(is.na(cPaths),arr.ind=T)[,1])
	if (length(toExclude)>0){
		cPaths=cPaths[-toExclude,]
	}
	Keggs=cPaths[which(cPaths[,3]<0.1),]
	v[[strsplit(i,'[.]')[[1]][2]]]=KEGGs[substr(Keggs[,1],1,5),3]
}

### Put into a matrix
m=matrix(0,nrow=length(unique(unlist(v))),ncol=length(names(v)))
rownames(m)=unique(unlist(v))
colnames(m)=gsub('Negat_','',names(v))
for (i in names(v)){
	m[rownames(v1[[i]]),i] = as.numeric(v1[[i]][,2])
}

### Select those to plot
metabolism = c(grep('metabolism',rownames(m),value=TRUE),
			 "Glycolysis / Gluconeogenesis",
			 "Pentose phosphate pathway",
			 "Lysosome",
			 "Autophagy - animal",
			 "Other types of O-glycan biosynthesis",
			 "Phagosome")
metabolism = sort(metabolism)

signaling = sort(grep('signaling',rownames(m),value=TRUE))
negative  = sort(names(which(m[,'Overlap']<0)))
# adhesion  = sort(setdiff(rownames(m),c(metabolism,signaling,negative)))

m=m[,c('Overlap','Intermediate','Long','Whole')]
colnames(m)=c('Overlap','Intermediate','Long','TAD-wide')
m=m[order(m[,1],m[,2],m[,3],m[,4],decreasing=TRUE),]
m=m[c(metabolism,signaling,negative),]
hc=colorRampPalette(c('navy','gray90','firebrick'))(21)

### Make heatmap
pdf('output/Heatmap.GseaTads.pdf',height=10)  ### FIGURE 1G
heatmap.3(m,col=hc,margins=c(22,28),breaks=seq(-2.5,2.5,length.out=length(hc)+1),
		  dendrogram='n',Colv='n',Rowv='n',KeyValueName='Enrichment score')
dev.off()


