
library(VennDiagram)

d=read.table('ForComparingGlassTCGA/ForVenn.txt',sep='\t')

v=NULL
v$Glass=c(1:d[1,1])
v$TCGA=c(1:d[3,1],(d[1,1]+1):(d[1,1]+d[5,1]))

venn.diagram(v,fill=c('blue','red'),file='ForComparingGlassTCGA/Venn.png') ### SUPPLEMENTARY FIGURE S1A


### Remove log files
for (i in list.files('ForComparingGlassTCGA',full.names=TRUE,pattern='log')){
	ii=strsplit(i,'[.]')[[1]]
	if (ii[length(ii)]=='log'){
		file.remove(i)
	}
}


