
source('~/Scripts/Heatmap3.r')

aml=gsub('.txt','',list.files('Hematopoiesis/output/AmlAnalyses/'))

m=NULL
for (i in aml){
	cG=read.table(paste('Hematopoiesis/output/AmlAnalyses/',i,'.txt',sep=''),header=TRUE)
	cFE=cG[,'FE']
	cFE[cG[,'FDR']>0.01]=0
	cFE[abs(cFE)<1]=0
	m=cbind(m,sign(cFE))
}

rownames(m)=cG[,1]
colnames(m)=aml

m=m[,c('IDH_vs_DNMT3A.Posit','IDH_vs_DNMT3A.Negat')]

red=NULL
for (i in rownames(m)){
	
	i=strsplit(i,'[.]')[[1]]	
	newRed=NULL
	
	for (EvGc in c('GC_Intron','IntronEvolCons')){
		cK=read.table(paste('output/ECO/ZScored_RMA_',i[1],'_vs_HSC.',i[2],'.KS_Test.',EvGc,'.txt',sep=''))
		if (as.numeric(cK[,2])>0.01){
			newRed=c(newRed,'white')
		} else if (cK[,3]=='Higher_Median'){
			newRed=c(newRed,'gold')
		} else {
			newRed=c(newRed,'forestgreen')
		}
	}
	
	cE=read.table(paste('output/RED/ZScored_RMA_',i[1],'_vs_HSC.',i[2],'.Hypergeometric.FoldEnrichment.txt',sep=''))
	cK=read.table(paste('output/RED/ZScored_RMA_',i[1],'_vs_HSC.',i[2],'.KStest.Dens.txt',sep=''))
	for (am in c('SINE/Alu','SINE/MIR')){
		if (cK[am,'Intron.Sense']>0.01){
			newRed=c(newRed,'white')
		} else if (cE[am,'Intron.Sense']>1){
			newRed=c(newRed,'gold')
		} else {
			newRed=c(newRed,'forestgreen')
		}
	}
	
	red=rbind(red,newRed)
}
colnames(red)=c('GC content','Evol. Cons.','Alu','MIR')

mPos=m[grep('Posit',rownames(m)),]
mNeg=m[grep('Negat',rownames(m)),]
rownames(mPos)=gsub('.Posit','',rownames(mPos))
rownames(mNeg)=gsub('.Negat','',rownames(mNeg))

redPos=red[grep('Posit',rownames(m)),]
redNeg=red[grep('Negat',rownames(m)),]
rownames(m)=gsub('.Posit',' - Up',gsub('.Negat',' - Down',rownames(m)))
rownames(red)=rownames(m)

wOrder=c('CMP','MEP','GMP','Erythroid','Mega','Eosin','Baso','Gran','Mono','Dendritic','BCell','NK','NKT','TCell')
m=m[c(paste(wOrder,'- Up'),paste(wOrder,'- Down')),]
red=red[rownames(m),]

wwOrder=names(which(rowSums(abs(m))!=0))

hc=c('navy','gray90','firebrick')
pdf('Hematopoiesis/output/Heatmap.AmlAnalyses.pdf') ### FIGURE 6H 6I
heatmap.3(m,col=hc,margins=c(10,25),cexCol=1,cexRow=1.2,dendrogram='n',Colv='n',Rowv='n',labRow=c(wOrder,wOrder),
		  RowSideColors=t(red),RowSideColorsSize=3.5,colsep=c(2,4,6),rowsep=length(wOrder))
heatmap.3(m[wwOrder,],col=hc,margins=c(15,25),cexCol=1,cexRow=1.2,
		  dendrogram='n',Colv='n',Rowv='n',labRow=gsub(' - Up','',gsub(' - Down','',wwOrder)))
dev.off()


