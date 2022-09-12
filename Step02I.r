
library(Hmisc)
library(alluvial)
library(data.table)
library(VennDiagram)
source('~/Scripts/Heatmap3.r')


##################
### CpGi to mC ###

d=read.table('output/CpGi_mC_Correspond.txt',stringsAsFactors=FALSE)

d23=table(paste(d[,2],d[,3]))
t23=cbind( matrix(unlist(strsplit(names(d23),' ')),byrow=TRUE,ncol=2), d23 )

t23=gsub('NotCorr','NS',t23)
t23=gsub('Corr','Sign',t23)

pdf('output/Alluvial.CpGi_mC_correspond.pdf',height=4,width=3) ### SUPPLEMENTARY FIGURE S1E
alluvial(t23[,c(1,2)],freq=as.numeric(t23[,3]),axis_labels=c('CpGi','mC'))
dev.off()


###### GSEAs

KEGGs=NULL ### This is to rename the pathways with lower-cases
for (i in list.files('HumanGenome/',pattern='PathwayLists')){
	cKegg=read.table(paste('HumanGenome/',i,sep=''),sep='\t',
					 stringsAsFactors=FALSE,colClasses="character")
	KEGGs=rbind(KEGGs,cKegg)
}
rownames(KEGGs)=KEGGs[,2]

myGSEAs=list.files('output/GSEA/CpGi')
v=NULL
v1=NULL
for (i in rev(myGSEAs)){
	cPaths=read.table(paste('output/GSEA/CpGi',i,'output/output.txt',sep='/'),
					  sep='\t',header=TRUE,stringsAsFactors=FALSE)
	toExclude=unique(which(is.na(cPaths),arr.ind=T)[,1])
	if (length(toExclude)>0){
		cPaths=cPaths[-toExclude,]
	}
	Keggs=cPaths[which(cPaths[,3]<0.1),]
	v[[strsplit(i,'[.]')[[1]][2]]]=KEGGs[substr(Keggs[,1],1,5),3]
	
	v1[[strsplit(i,'[.]')[[1]][2]]]=Keggs[,1:2]
	rownames(v1[[strsplit(i,'[.]')[[1]][2]]])=KEGGs[substr(Keggs[,1],1,5),3]
}

m=matrix(0,nrow=length(unique(unlist(v))),ncol=length(names(v)))
rownames(m)=unique(unlist(v))
colnames(m)=gsub('RankTads.','',names(v))
for (i in names(v)){
	m[rownames(v1[[i]]),i] = as.numeric(v1[[i]][,2])
}

m=m[,c('Proximal','Intermediate','Long','Whole')]
colnames(m)=c('Proximal','Intermediate','Long','TAD-wide')
m=m[order(m[,1],m[,2],m[,3],m[,4],decreasing=TRUE),]

w=c("Cytokine-cytokine receptor interaction",
	"Phosphatidylinositol signaling system",
	"Necroptosis",
	"Biosynthesis of unsaturated fatty acids",
	"Endocytosis",
	"Fructose and mannose metabolism",
	"Inositol phosphate metabolism",
	"Lysosome",
	"N-Glycan biosynthesis",
	"Peroxisome",
	"Phagosome",
	"Sphingolipid metabolism",
	"Steroid biosynthesis",
	"Ribosome",
	"Spliceosome")
m=m[w,]

hc=colorRampPalette(c('navy','gray90','firebrick'))(21)

pdf('output/Heatmap.CpGiGseaTads.pdf',height=10) ### SUPPLEMENTARY FIGURE S1G
heatmap.3(m,col=hc,margins=c(32,28),breaks=seq(-2.5,2.5,length.out=length(hc)+1),
		  dendrogram='n',Colv='n',Rowv='n',KeyValueName='Enrichment score')
dev.off()



