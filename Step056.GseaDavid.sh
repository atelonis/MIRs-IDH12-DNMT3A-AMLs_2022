#!/bin/bash 



######################
###====  GSEA  ====###
######################

[ ! -d TCGA/output/GSEA ] && mkdir TCGA/output/GSEA
tads=`ls TCGA/output/Correlation/Step3/ | grep RankTads | sed 's/.txt//'`

for t in $tads
do
	echo $t
	[ ! -d TCGA/output/GSEA/Select.$t ] && mkdir TCGA/output/GSEA/Select.$t
	Scripts/GSEA.sh selKEGG_Ensembl TCGA/output/GSEA/$t/Genelist.rnk TCGA/output/GSEA/Select.$t/TCGA/output > TCGA/output/GSEA/Select.$t/Log.txt
done



#######################
###====  DAVID  ====###
#######################

mkdir TCGA/output/DAVID
tadSets=`ls TCGA/output/Correlation/Step3/ | grep TADsets | sed 's/.txt//' | grep Negat`
cat TCGA/output/Correlation/Step1/Correlations.txt | awk '$3!="NA"' | cut -f2 | sort | uniq > TCGA/output/DAVID/Background.TADs_WholeNames.txt
cat TCGA/output/DAVID/Background.TADs_WholeNames.txt | cut -d '|' -f1 | sort | uniq > TCGA/output/DAVID/Background.TADs.txt

for t in $tadSets
do
	echo $t
	mkdir TCGA/output/DAVID/$t
	cat TCGA/output/Correlation/Step3/$t.txt | cut -d '|' -f1 | sort | uniq > TCGA/output/DAVID/$t/Genes.txt
	
	Genes="output/DAVID/$t/Genes.txt"
	Background="output/DAVID/Background.TADs.txt"
	Idtype="ENSEMBL_GENE_ID"
	Output="output/DAVID/$t/Output.txt"
	LogFile="output/DAVID/$t/Log.txt"
	perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR
done

### All genes
echo Whole
mkdir TCGA/output/DAVID/TADsets.Negat_Whole
cat TCGA/output/DAVID/TADsets.Negat_[OIL]*/Genes.txt | sort | uniq > TCGA/output/DAVID/TADsets.Negat_Whole/Genes.txt
	
Genes="output/DAVID/TADsets.Negat_Whole/Genes.txt"
Background="output/DAVID/Background.TADs.txt"
Idtype="ENSEMBL_GENE_ID"
Output="output/DAVID/TADsets.Negat_Whole/Output.txt"
LogFile="output/DAVID/TADsets.Negat_Whole/Log.txt"
perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile

OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
cat $Output | head -n1 > $OutputFDR
cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR



##################
###   W   sets ###
##################

wSets=`ls TCGA/output/GSEA/ | grep RankTads | cut -d '.' -f2`
for t in $wSets
do
	echo $t
	mkdir TCGA/output/DAVID/TAD_W.$t
	cat TCGA/output/GSEA/RankTads.$t/Genelist.rnk | sort -k2,2n | head -n1000 | cut -f1 | sort | uniq > TCGA/output/DAVID/TAD_W.$t/Genes.txt
	cat TCGA/output/GSEA/RankTads.$t/Genelist.rnk | cut -f1 | sort | uniq > TCGA/output/DAVID/TAD_W.$t/Background.txt
	
	Genes="output/DAVID/TAD_W.$t/Genes.txt"
	Background="output/DAVID/TAD_W.$t/Background.txt"
	Idtype="ENSEMBL_GENE_ID"
	Output="output/DAVID/TAD_W.$t/Output.txt"
	LogFile="output/DAVID/TAD_W.$t/Log.txt"
	perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR
done

