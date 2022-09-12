
library(data.table)

source('Scripts/Heatmap3.r')

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
	return(c(log2(myFE),p))
}

############################
### Overlap with repeats ###
############################

d=read.table('output/TAD.C_in_Repeats.txt')
d=d[which(rowSums(d[,2:4])>10),]
d=d[setdiff(rownames(d),c('Low_complexity','Simple_repeat')),]

fe=NULL
for (i in c(1:4)){
	fe=cbind(fe,d[,i]/(d[1,i]+d[2,i]))
}

for (i in c(2,3,4)){
	fe[,i]=fe[,i]/fe[,1]
}
rownames(fe)=rownames(d)
colnames(fe)=colnames(d)

vO=NULL
vI=NULL
vL=NULL
for (i in c(3:nrow(d))){
	vO=rbind(vO,c(rownames(d)[i],hypergeom.test(d[i,4],d[1,4]+d[2,4],d[i,1],d[1,1]+d[2,1])))
	vI=rbind(vI,c(rownames(d)[i],hypergeom.test(d[i,2],d[1,2]+d[2,2],d[i,1],d[1,1]+d[2,1])))
	vL=rbind(vL,c(rownames(d)[i],hypergeom.test(d[i,3],d[1,3]+d[2,3],d[i,1],d[1,1]+d[2,1])))
}

w=rbind( cbind(rep('Proximal',nrow(vO)),vO),
		 cbind(rep('Intermediate',nrow(vI)),vI),
		 cbind(rep('Long',nrow(vL)),vL) )

write.table(w,file='output/TAD.C_in_Repeats.FePVals.txt',sep='\t',row.names=FALSE,
			col.names=c('Gene set','Repeatitive element','Fold enrichment','P-value'))

vFE=cbind(as.numeric(vO[,2]),as.numeric(vI[,2]),as.numeric(vL[,2]))
vPV=cbind(as.numeric(vO[,3]),as.numeric(vI[,3]),as.numeric(vL[,3]))
rownames(vFE)=vO[,1]
colnames(vFE)=c('Overlap','Intermediate','Long')

vFE[vPV>0.015]=0
vFE=vFE[which(rowSums(vFE!=0)>0),]

hc=colorRampPalette(c('navy','gray90','firebrick'))(51)
pdf('output/TAD.C_in_Repeats.pdf') ### FIGURE 2F
heatmap.3(vFE,col=hc,breaks=seq(-1.5,1.5,length.out=52),margins=c(23,28),Colv='n',dendrogram='row',
		  KeyValueName='Fold enrichment')
dev.off()


#################################
### Evolutionary conservation ###
#################################

ev=NULL
for (i in c(1:22,'Y')){
	print(i)
	ev=rbind(ev,fread(paste('bed/BedIntersect.TadCs.EvolCons_chr',i,'.bed',sep='')))
}
w=which( (ev[,3]-101==ev[,9])[,1] )
evC=ev[w,]


corTad=read.table('output/Correlation/Step2/TadCorrelations.txt',header=TRUE,stringsAsFactors=FALSE)
corTad[,1]=gsub('.NA0','',corTad[,1])

evBackground=ecdf(as.numeric(evC[,10]))


Xs=seq(-10,10,length.out=101)
dBins=c('Overlap','Intermediate','Long')

binCol=c('#F8C1B8','#FF6655','#A6093D')
names(binCol)=dBins

pdf('output/TAD.EvolCons_Cs.pdf',height=4,width=4) ### FIGURE 2E
par(mar=c(5,5,2,2))

# Per distance
plot(0,0,pch=NA,xlim=c(-10,10),ylim=c(-0.1,0.1),xlab='Evolutionary conservation of mC',
	 ylab='Cumulative fraction difference\nfrom background distribution')
lines(Xs,(evBackground(Xs)*0.2)-0.1,lwd=2,lty=2,col='gray')
for (i in dBins){
	print(i)
	cC=intersect(rownames(evC),unique(corTad[which(corTad[,8]==i),1]))
	cEV=as.numeric(evC[cC,10])
	cKS=ks.test(cEV,evBackground)
	lines(Xs,evBackground(Xs)-ecdf(cEV)(Xs),col=binCol[i],lwd=5)
	print(cKS)
	axis(4,c(-0.1,-0.05,0,0.05,0.1),c(0,0.25,0.5,0.75,1),col='gray',col.axis='gray')
}
lines(c(-10,10),c(0,0),col='black',lwd=3)
legend('bottomright',names(binCol),col=binCol,lwd=3,pch=NA,cex=0.75)

dev.off()






