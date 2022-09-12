
library(methylSig)

dir.create('output/MethylSig',showWarnings=FALSE)

fileList=c(
# IDH (n=9)
'Datasets/GSM2310848_methylcall.Sample_2268.mincov10.txt',
'Datasets/GSM2310867_methylcall.Sample_4335.mincov10.txt',
'Datasets/GSM2310873_methylcall.Sample_5289.mincov10.txt',
'Datasets/GSM2310888_methylcall.Sample_6887.mincov10.txt',
'Datasets/GSM2310894_methylcall.Sample_7059.mincov10.txt',
'Datasets/GSM2310916_methylcall.Sample_7180.mincov10.txt',
'Datasets/GSM2310919_methylcall.Sample_7301.mincov10.txt',
'Datasets/GSM2310934_methylcall.Sample_7418.mincov10.txt',
'Datasets/GSM2310935_methylcall.Sample_7420.mincov10.txt',
# DNMT3A (n=16)
'Datasets/GSM2310823_methylcall.Sample_2185.mincov10.txt',
'Datasets/GSM2310827_methylcall.Sample_2195.mincov10.txt',
'Datasets/GSM2310831_methylcall.Sample_2206.mincov10.txt',
'Datasets/GSM2310835_methylcall.Sample_2234.mincov10.txt',
'Datasets/GSM2310865_methylcall.Sample_3331.mincov10.txt',
'Datasets/GSM2310866_methylcall.Sample_4334.mincov10.txt',
'Datasets/GSM2310871_methylcall.Sample_5286.mincov10.txt',
'Datasets/GSM2310878_methylcall.Sample_6239.mincov10.txt',
'Datasets/GSM2310881_methylcall.Sample_6374.mincov10.txt',
'Datasets/GSM2310906_methylcall.Sample_7131.mincov10.txt',
'Datasets/GSM2310917_methylcall.Sample_7188.mincov10.txt',
'Datasets/GSM2310922_methylcall.Sample_7313.mincov10.txt',
'Datasets/GSM2310924_methylcall.Sample_7318.mincov10.txt',
'Datasets/GSM2310926_methylcall.Sample_7322.mincov10.txt',
'Datasets/GSM2310930_methylcall.Sample_7407.mincov10.txt',
'Datasets/GSM2310932_methylcall.Sample_7411_repeat.mincov10.txt')

sample.id=c( paste('IDH',c(1:9),sep='_'), paste('DNMT3A',c(1:16),sep='_') )
treatment=c( rep(0,9), rep(1,16) )

meth = methylSigReadData(fileList, sample.ids=sample.id, assembly="hg19", treatment=treatment,
						 context="CpG", destranded=TRUE, minCount=10, num.cores=50)

myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=c(16,9))
myBackground = paste(as.character(myDiffSigboth[,"chr"]),'.',
					 as.character(myDiffSigboth[,"start"]),'-',
					 as.character(myDiffSigboth[,"end"]),sep='')

myHyper = myDiffSigboth[myDiffSigboth[,"qvalue"]<0.05 & myDiffSigboth[,"meth.diff"]>25,]
myHypo  = myDiffSigboth[myDiffSigboth[,"qvalue"]<0.05 & myDiffSigboth[,"meth.diff"]<(-25),]

myHyper = paste(as.character(myHyper[,"chr"]), '.', 
				as.character(myHyper[,"start"]), '-', 
				as.character(myHyper[,"end"]), sep='')

myHypo  = paste(as.character(myHypo[,"chr"]),  '.', 
				as.character(myHypo[,"start"]),  '-', 
				as.character(myHypo[,"end"]),  sep='')

write.table(myBackground,file='output/MethylSig/Background.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(myHyper,file='output/MethylSig/Hyper.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(myHypo ,file='output/MethylSig/Hypo.txt', quote=FALSE,row.names=FALSE,col.names=FALSE)
# Hyper/Hypo in dnmt3a mutants as compared to idh


