

### Gene symbols lookup table
lookup=read.table('Microarrays/LookupIDs.txt',stringsAsFactors=FALSE,
				  sep='\t',header=TRUE)
lookup=unique(paste(lookup[,2],lookup[,3],lookup[,4],sep='|'))
lookup=matrix(unlist(strsplit(lookup,'[|]')),ncol=3,byrow=TRUE)
rownames(lookup)=lookup[,2]

### Correlations
signCorrs=read.table('output/Correlation/Step2/TadCorrelations.txt',sep='\t',stringsAsFactors=FALSE,
					 header=TRUE)
signCorrs[,8]=gsub('Overlap','Proximal',signCorrs[,8])

### Gene characteristics
cha=read.table('output/GC_Ccount/TAD.GeneTable.txt',sep='\t',check.names=FALSE)
mir=cha[,'SINE/MIR']
names(mir)=rownames(cha)

cpgd = 1000 * cha[,'C_Count'] / ( (2**cha[,'Exon_Length']) + (2**cha[,'Intron_Length']) )
names(cpgd)=rownames(cha)

s=NULL
zd=NULL
for (i in c('Overlap','Intermediate','Long','Whole')){
	cRank=read.table(paste('output/Correlation/Step3/RankTads',i,'txt',sep='.'),stringsAsFactors=FALSE)
	cGenes=lookup[as.character(cRank[,1]),1]
	rownames(cRank)=cGenes
	cGenes=intersect(names(mir),cGenes)
	
	cMIR=mir[cGenes]
	cZero=names(which( cMIR == 0 ))
	mDens=names(which( cMIR >= quantile(as.numeric(cMIR),0.9) ))
	
	i=gsub('Overlap','Proximal',i)
	zd[[ paste( i, 'Dense', sep='.') ]] = cRank[mDens,2]
	zd[[ paste( i, 'Zero', sep='.') ]] = cRank[cZero,2]
	
	cU=wilcox.test(cRank[cZero,2],cRank[mDens,2])
	if (cU$p.value<0.01){
		s=c(s,'','*')
	} else {
		s=c(s,'','ns')
	}
	
	oZD=length(names(zd))
	while( length(names(zd)) <= oZD+2 ){
		takeSample=sample(cZero,length(mDens))
		cpgU=wilcox.test(cpgd[cZero],cpgd[takeSample])
		if (cpgU$p.value > 0.1){
			zd[[ paste( i, '.Iter', length(names(zd))-oZD+1, sep='' ) ]] = cRank[takeSample,2]
		
			cU=wilcox.test(cRank[mDens,2],cRank[takeSample,2])
			if (cU$p.value<0.01){
				s=c(s,'*')
			} else {
				s=c(s,'ns')
			}
		}
	}
	
	if (i=='Whole'){
		cCorrs=matrix(unlist(strsplit(signCorrs[,2],'[|]')),ncol=2,byrow=TRUE)
	} else {
		cCorrs=matrix(unlist(strsplit(signCorrs[which(signCorrs[,8]==i),2],'[|]')),ncol=2,byrow=TRUE)
	}
	
	zd[[ paste( i, 'Corrs', sep='.' ) ]] = cRank[cCorrs,2]
	
	cU=wilcox.test(cRank[cCorrs,2],cRank[mDens,2])
	if (cU$p.value<0.01){
		s=c(s,'*')
	} else {
		s=c(s,'ns')
	}
}


myCols=rep(c('skyblue','palegreen','gray','gray','gray','pink'),4)

pdf('output/Boxplot.ByMIRdensity.pdf',height=4,width=6) ### SUPPLEMENTARY FIGURE S2A
par(mar=c(6,5,2,2))
boxplot(zd,col=myCols,las=2,yaxt='n',ylab='Maximum absolute coefficient',mgp=c(1.6,0.6,0),ylim=c(0,1),
		names=matrix(unlist(strsplit(names(zd),'[.]')),ncol=2,byrow=TRUE)[,2])
axis(2,mgp=c(1.6,0.6,0))
text( c(1:length(names(zd))), 0.95, s,col='red')
dev.off()



