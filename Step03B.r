
library(data.table)

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

back=read.table('output/Correlation/Step4/Reps.Cs.Background.txt',sep='\t',stringsAsFactors=FALSE)
corrs=read.table('output/Correlation/Step4/Reps.Cs.Significant.txt',sep='\t',stringsAsFactors=FALSE)

mirBack=unique(grep('SINE.MIR',back[,1],value=TRUE))
mirCorrs=unique(grep('SINE.MIR',corrs[,1],value=TRUE))


#######################
### B-Boxes in MIRs ###

bbox=read.table('HumanGenome/MIRs_with_BBox.txt',stringsAsFactors=FALSE,sep='\t')[,1]

k=length(intersect(bbox,mirCorrs))
s=length(mirCorrs)
M=length(intersect(bbox,mirBack))
N=length(mirBack)

bHT=hypergeom.test(k,s,M,N)
print('BBox')
print(bHT)

pdf('output/MIR_EnrichBBox.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S2B
par(mar=c(5,15,2,2))
barplot(bHT[1],names.arg='B-Box',ylim=c(0,2),col='gray',las=2,ylab='Fold enrichment',main='B-Box',mgp=c(1.6,0.6,0),
		yaxt='n',border='NA')
axis(2,mgp=c(1.6,0.6,0))
dev.off()


################################
### Distance to TAD boundary ###

back =log2(back[grep('SINE.MIR',back[,1]),3])
corrs=log2(corrs[grep('SINE.MIR',corrs[,1]),3])

Xs=seq(0,max(back)*1.2,length.out=50)

ecdfBack=ecdf(back)
bb=ecdfBack(Xs)

cY=-0.2

pdf('output/MIR_TADbound_distance.pdf',width=4,height=4) ### SUPPLEMENTARY FIGURE S2C
par(mar=c(4,4,3,3))

plot(Xs,(-2*cY*bb)+cY,type='l',lwd=2,pch=NA,col='gray60',lty=2,
	 xlab='Distance to closest TAD boundary (log2 bp)',
	 ylab='Cumulative fraction difference\nfrom background distribution',
	 xlim=c(0,max(Xs)),ylim=c(cY,-cY),mgp=c(1.6,0.6,0))
axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray60',col.axis='gray60',mgp=c(1.6,0.6,0))
mtext('Cumulative background distribution',4,col='gray60',line=1.6)

L= bb - ecdf(corrs)(Xs)
lines(Xs,L,col='goldenrod',lwd=5)

lines(c(0,max(Xs)),c(0,0),col='black',lwd=2)
dev.off()


########################
### Divergence score ###

diver=as.matrix(fread('HumanGenome/MIR.Divergence.txt',sep='\t'))
rownames(diver)=diver[,1]
diver=diver[mirBack,]


pdf('output/MIR_mC_Divergence.pdf',width=4,height=4) ### SUPPLEMENTARY FIGURE S2E
par(mar=c(4,4,3,3))

Xs=seq(0,100)

ecdfBack=ecdf( as.numeric(diver[,2]) + as.numeric(diver[,3]) + as.numeric(diver[,4]) )
corrDiver= as.numeric(diver[mirCorrs,2]) + as.numeric(diver[mirCorrs,3]) + as.numeric(diver[mirCorrs,4])
bb=ecdfBack(Xs)

print(ks.test(corrDiver,ecdfBack))

cY=-0.2
plot(Xs,(-2*cY*bb)+cY,type='l',lwd=2,pch=NA,col='gray60',lty=2,mgp=c(1.6,0.6,0),
	 xlab='Divergence', ylab='Cumulative fraction difference\nfrom background distribution',
	 xlim=c(0,max(Xs)),ylim=c(cY,-cY))
axis(4,c(cY,cY/2,0,-cY/2,-cY),c(0,0.25,0.5,0.75,1),col='gray60',col.axis='gray60',mgp=c(1.6,0.6,0))
mtext('Cumulative background distribution',4,col='gray60',line=1.6)

L= bb - ecdf(corrDiver)(Xs)
lines(Xs,L,col='goldenrod',lwd=5)
lines(c(0,max(Xs)),c(0,0),col='black',lwd=2)

dev.off()


#######################
### MIR subfamilies ###

d=read.table('output/MIR_Subfamilies.txt',sep='\t',stringsAsFactors=FALSE)

d=d[grep('Posit',d[,1],invert=TRUE),]
d=d[grep('TADsets_W',d[,1],invert=TRUE),]
d[,1]=gsub('TADsets.','',gsub('Overlap','Proximal',d[,1]))
d=d[nrow(d):1,]

pdf('output/MIR_Subfamilies.pdf',height=5,width=5) ### SUPPLEMENTARY FIGURE S2F
par(mar=c(5,8,2,2))
barplot(as.numeric(d[,7]),col='gray',border=NA,mgp=c(1.6,0.6,0),horiz=TRUE,las=2,xaxt='n',
		xlab='Fold enrichment',xlim=c(0,3),names.arg=d[,2])
axis(1,mgp=c(1.6,0.6,0))
lines(c(0.5,0.5),c(0,30),col='red',lty=2,lwd=1.5)
lines(c(2,2),c(0,30),col='red',lty=2,lwd=1.5)
dev.off()


##################
### Insulators ###

d=read.table('output/MIR_EnrichInsulators.txt',stringsAsFactors=FALSE)

d=d[grep('Posit',d[,1],invert=TRUE),]
d=d[grep('W',d[,1],invert=TRUE),]
d[,1]=gsub('_','.',d[,1])
d[,1]=gsub('Negat.','',gsub('TADsets.','',gsub('Overlap','Proximal',d[,1])))

pdf('output/MIR_EnrichInsulators.pdf',height=4,width=4) ### SUPPLEMENTARY FIGURE S2D
par(mar=c(8,14,2,2))
#bp=barplot(log2(d[,6]),col='black',xlim=c(-6,1),horiz=TRUE,mgp=c(1.6,0.6,0),xlab='Fold enrichment (log2)',
bp=barplot(p.adjust(d[,7],method='fdr'),col='gray',ylim=c(0,0.6),mgp=c(1.6,0.6,0),ylab='FDR',
		   names.arg=d[,1],las=2,yaxt='n')
axis(2,mgp=c(1.6,0.6,0))
dev.off()


#####################
### Histone marks ###

d=read.table('output/MIR_EnrichHistones.txt',stringsAsFactors=FALSE) ### SUPPLEMENTARY FIGURE S3

d=d[c(grep('Negat',d[,2]),which(d[,2]=='mCs')),]
d[,2]=gsub('TADsets.Negat_','',gsub('Overlap','Proximal',d[,2]))
d=d[c( which(d[,2]=='Proximal'), which(d[,2]=='Intermediate'), which(d[,2]=='Long'), which(d[,2]=='mCs')), ]
d=d[nrow(d):1,]


pdf('output/MIR_EnrichHistones.pdf',height=4,width=4)
par(mar=c(5,10,2,2))
bp=barplot(d[,7],horiz=TRUE,las=2,xaxt='n',names.arg=d[,1],col='gray',xlim=c(0,3),
		   xlab='Fold enrichment',mgp=c(1.6,0.6,0),border=NA)
lines(c(0.5,0.5),c(0,20),col='red',lty=2,lwd=1.5)
lines(c(2,2),c(0,20),col='red',lty=2,lwd=1.5)
axis(1,mgp=c(1.6,0.6,0))
dev.off()


##############
### Encode ###

d=read.table('output/MIR_EnrichEncode.txt',stringsAsFactors=FALSE)
d=cbind(d,p.adjust(d[,8],method='fdr'))
d[,2]=gsub('TADsets.','',gsub('_','.',d[,2]))
d[,2]=gsub('Overlap','Proximal',d[,2])

d=d[-which(d[,2]=='W.Whole'),]
d=d[-grep('Posit',d[,2]), ]

m=matrix(0,nrow=length(unique(d[,1])),ncol=length(unique(d[,2])))
rownames(m)=sort(unique(d[,1]))
colnames(m)=sort(unique(d[,2]))

f=m

for (i in c(1:nrow(d))){
	m[ d[i,1], d[i,2] ] = log2(d[i,7])
	f[ d[i,1], d[i,2] ] = d[i,9]
}
toExclude=-unique(which(is.na(m),arr.ind=T)[,1])

f=f[toExclude,]
m=m[toExclude,]

inf2zero=which(is.infinite(m) & f<0.05)
inf2minus=which(is.infinite(m))

m[inf2zero]=0
m[inf2minus]=-4


w=NULL
### genes
for (i in grep('Negat',colnames(m),value=TRUE)){
#	m=m[order(m[,i]),]
	iCols=rep('gray',nrow(m))
	iCols[which(m[,i]>1 & f[,i]<0.05 & m[,gsub('Negat','W',i)]<1)]='firebrick'
	forW=sort(rownames(m)[which(iCols=='firebrick')])
	w=rbind(w,cbind(rep( i, length(forW)), forW, m[forW,i] ))
}

### mCs
i='mCs'
iCols=rep('gray',nrow(m))
iCols[which(m[,i]>1 & f[,i]<0.05)]='firebrick'
forW=sort(rownames(m)[which(iCols=='firebrick')])
w=rbind(w,cbind(rep( i, length(forW)), forW, m[forW,i] ))

w2=NULL
for (i in c('Negat.Proximal','Negat.Intermediate','Negat.Long','mCs')){
	w2=rbind(w2,w[which(w[,1]==i),])
}
w2=w2[nrow(w2):1,]

pdf('output/MIR_EnrichEncode.pdf',height=4,width=4)
par(mar=c(4,5,1,1))
bp=barplot(as.numeric(w2[,3]),names.arg=w2[,2],xlim=c(0,2.5),horiz=TRUE,mgp=c(1.6,0.6,0),las=2,xaxt='n',
		xlab='Fold enrichment (log2)',col='black',cex.names=0.8)
axis(1,mgp=c(1.6,0.6,0))
text(as.numeric(w2[,3])+0.1,bp,'*',cex=1.5)
dev.off()




