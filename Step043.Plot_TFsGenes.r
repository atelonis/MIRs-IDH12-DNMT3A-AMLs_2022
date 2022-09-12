
library(amap)
library(dendextend)

source('Scripts/Heatmap3.r')

hc=colorRampPalette(c('navy','gray90','firebrick'))(51)

margs=c(25,22,10)
names(margs)=c('Beck','Encode')

setNames=c('Proximal','Intermediate','Long','W')
names(setNames)=c('Negat_Overlap','Negat_Intermediate','Negat_Long','W_Overlap')

gataType  = sort(c('GATAD2A','GATAD2B','MTA2','MTA1'))
heLoHe    = sort(c('MXI1','MAX','ARNT','BHLHE40','MITF','SREBF1'))
rnaBind   = sort(c('FUS','JUN','POLR2A','U2AF1','XRCC5','AGO1','FIP1L1','HNRNPK','HNRNPL','HDAC2','PCBP1','SAFB','UBTF'))
forEncode = c(gataType,heLoHe,rnaBind) ### See (*) at end of this script

allM=NULL
for (be in c('Beck','Encode')){
	cTable=read.table(paste('output/TAD.',be,'TFs.txt',sep=''),sep='\t',stringsAsFactors=FALSE,header=TRUE)
	cTable=cTable[-which(cTable[,1]=='W_Long' | cTable[,1]=='W_Intermediate'),]
	m=matrix(0,nrow=length(unique(cTable[,1])),ncol=length(unique(cTable[,2])))
	rownames(m)=sort(unique(cTable[,1]))
	colnames(m)=sort(unique(cTable[,2]))
	
	for (i in rownames(m)){
		for (j in colnames(m))
			m[i,j]=as.numeric(cTable[which(cTable[,1]==i & cTable[,2]==j),'ZScore'])
	}
	m=t(m)
	colnames(m)=setNames[colnames(m)]
	
	pdf(paste('output/TAD.',be,'TFs.pdf',sep='')) ### FIGURE 3A 3B
	for (cMeth in c('euclidean','manhattan','pearson','spearman','kendall')){
		dRow=as.dendrogram(hcluster(m,method=cMeth))
		dCol=as.dendrogram(hcluster(t(m),method=cMeth))
		heatmap.3(m,col=hc,breaks=seq(-20,20,length.out=length(hc)+1),margins=c(margs[be],25),
				  Colv=dCol,Rowv=dRow,main=cMeth)
		if (be=='Encode'){
			lr=rep('',nrow(m))
			w=which(m[,1]>10 & m[,2]>10 & m[,3]<10 & m[,4]<10)
			lr[w]=rownames(m)[w]
			
			s=matrix('',nrow=length(forEncode),ncol=ncol(m))
			s[which(abs(m[forEncode,])>10)]='*'
			dCol=as.dendrogram(hcluster(t(m[forEncode,]),method=cMeth))
			heatmap.3(m[forEncode,],col=hc,breaks=seq(-20,20,length.out=length(hc)+1),margins=c(margs[be]+5,25),
					  Colv=dCol,Rowv='n',dendrogram='column',main=cMeth,cellnote=s,notecol='white',notecex=1.5)
		}
	}
	allM[[be]]=m
	dev.off()
}

### m is for ENCODE
m[m<10]=0
m[m>0]=1

toPlot=colSums(m)[c('W','Long','Intermediate','Proximal')]

pdf('output/Barplot.TadEncode.pdf',height=4,width=4) ### FIGURE 3B
par(mar=c(8,14,2,2))
barplot(rev(toPlot),las=2,col='black',yaxt='n',ylab='Number of TFs',ylim=c(0,150),border=NA)
axis(2,c(0,50,100,150),c(0,50,100,150))

dev.off()

### SELECTED (*)
### These clusters were extracted from DAVID after analyzing the ENSGs from the 'w' list:
###
###GATA-type zinc finger
###ENSG00000167491	GATA zinc finger domain containing 2A(GATAD2A)	RG	Homo sapiens
###ENSG00000143614	GATA zinc finger domain containing 2B(GATAD2B)	RG	Homo sapiens
###ENSG00000149480	metastasis associated 1 family member 2(MTA2)	RG	Homo sapiens
###ENSG00000182979	metastasis associated 1(MTA1)	RG	Homo sapiens
###
###Helix-loop-helix
###ENSG00000119950	MAX interactor 1, dimerization protein(MXI1)	RG	Homo sapiens
###ENSG00000125952	MYC associated factor X(MAX)	RG	Homo sapiens
###ENSG00000143437	aryl hydrocarbon receptor nuclear translocator(ARNT)	RG	Homo sapiens
###ENSG00000134107	basic helix-loop-helix family member e40(BHLHE40)	RG	Homo sapiens
###ENSG00000187098	melanogenesis associated transcription factor(MITF)	RG	Homo sapiens
###ENSG00000072310	sterol regulatory element binding transcription factor 1(SREBF1)	RG	Homo sapiens
###
###RNA binding
###ENSG00000089280	FUS RNA binding protein(FUS)	RG	Homo sapiens
###ENSG00000177606	Jun proto-oncogene, AP-1 transcription factor subunit(JUN)	RG	Homo sapiens
###ENSG00000181222	RNA polymerase II subunit A(POLR2A)	RG	Homo sapiens
###ENSG00000160201	U2 small nuclear RNA auxiliary factor 1(U2AF1)	RG	Homo sapiens
###ENSG00000079246	X-ray repair cross complementing 5(XRCC5)	RG	Homo sapiens
###ENSG00000092847	argonaute 1, RISC catalytic component(AGO1)	RG	Homo sapiens
###ENSG00000145216	factor interacting with PAPOLA and CPSF1(FIP1L1)	RG	Homo sapiens
###ENSG00000165119	heterogeneous nuclear ribonucleoprotein K(HNRNPK)	RG	Homo sapiens
###ENSG00000104824	heterogeneous nuclear ribonucleoprotein L(HNRNPL)	RG	Homo sapiens
###ENSG00000196591	histone deacetylase 2(HDAC2)	RG	Homo sapiens
###ENSG00000169564	poly(rC) binding protein 1(PCBP1)	RG	Homo sapiens
###ENSG00000160633	scaffold attachment factor B(SAFB)	RG	Homo sapiens
###ENSG00000108312	upstream binding transcription factor, RNA polymerase I(UBTF)	RG	Homo sapiens
###


