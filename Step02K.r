

muts=read.table('output/Correlation/Step2/TadCorrelations.txt',header=TRUE,stringsAsFactors=FALSE)
wtwt=read.table('output/Correlation/Step2/WtWtTadCorrelations.txt',stringsAsFactors=FALSE)

wtwt=wtwt[-which(is.na(wtwt[,3])),]

v=NULL
v[['Observed']]=abs(muts[,3])
v[['WT/WT']]=abs(wtwt[,3])

pdf('output/Barplot.CompareToWtWt.pdf',height=4,width=4)  ### SUPPLEMENTARY FIGURE S1H
par(mar=c(5,12,2,2),mgp=c(1.6,0.6,0))
boxplot(v,names=c('DNMT3A\nor IDH1/2','WT/WT'),col='gray',ylim=c(0,1),las=2,yaxt='n',
		mgp=c(1.6,0.6,0),ylab='Spearman correlation coefficient')
axis(2,mgp=c(1.6,0.6,0))
dev.off()



