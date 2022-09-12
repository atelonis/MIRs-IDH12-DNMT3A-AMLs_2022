#!/bin/bash 

################################
### For plotting a bit later ###
################################

mkdir output/ForCircos
cat output/Correlation/Step1/TadCorrelations.txt | cut -f1 | sort | uniq | tr '.' '\t' | awk '{print $1,$2-1,$2}' | tr ' ' '\t' > output/ForCircos/background.txt
cat output/Correlation/Step2/TadCorrelations.txt | grep SNX11 | cut -f1 | sort | uniq | tr '.' '\t' | awk '{print $1,$2-1,$2}' | tr ' ' '\t' > output/ForCircos/snx11.txt
cat output/Correlation/Step2/TadCorrelations.txt | grep HOXB5 | cut -f1 | sort | uniq | tr '.' '\t' | awk '{print $1,$2-1,$2}' | tr ' ' '\t' > output/ForCircos/hoxb5.txt


##############
###  GSEA  ###
##############

[ ! -d output/GSEA ] && mkdir output/GSEA
tads=`ls output/Correlation/Step3/ | grep RankTads | sed 's/.txt//'`

for t in $tads
do
	echo $t
	
	[ ! -d output/GSEA/$t ] && mkdir output/GSEA/$t
	Scripts/GSEA.sh selKEGG_Entrez output/GSEA/$t/Genelist.rnk output/GSEA/$t/output > output/GSEA/$t/Log.txt
done


###############
###  DAVID  ###
###############

tadSets=`ls output/Correlation/Step3/ | grep TADsets | sed 's/.txt//' | grep Negat`
cat output/Correlation/Step1/TadCorrelations.txt | cut -f2 | sort | uniq > output/DAVID/Background.TADs_WholeNames.txt
cat output/Correlation/Step1/TadCorrelations.txt | cut -f2 | cut -d '|' -f2 | sort | uniq > output/DAVID/Background.TADs.txt

for t in $tadSets
do
	echo $t
	mkdir output/DAVID/$t
	cat output/Correlation/Step3/$t.txt | cut -d '|' -f2 | sort | uniq > output/DAVID/$t/Genes.txt
	
	Genes="output/DAVID/$t/Genes.txt"
	Background="output/DAVID/Background.TADs.txt"
	Idtype="ENTREZ_GENE_ID"
	Output="output/DAVID/$t/Output.txt"
	LogFile="output/DAVID/$t/Log.txt"
	perl Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR
done

### All genes
mkdir output/DAVID/TADsets.Negat_Whole
cat output/DAVID/TADsets.Negat_[OIL]*/Genes.txt | sort | uniq > output/DAVID/TADsets.Negat_Whole/Genes.txt
	
Genes="output/DAVID/TADsets.Negat_Whole/Genes.txt"
Background="output/DAVID/Background.TADs.txt"
Idtype="ENTREZ_GENE_ID"
Output="output/DAVID/TADsets.Negat_Whole/Output.txt"
LogFile="output/DAVID/TADsets.Negat_Whole/Log.txt"
perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile

OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
cat $Output | head -n1 > $OutputFDR
cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR



#########################
###  DAVID   W   sets ###
#########################

wSets=`ls output/GSEA/ | grep RankTads | cut -d '.' -f2`
for t in $wSets
do
	echo $t
	mkdir output/DAVID/TAD_W.$t
	cat output/GSEA/RankTads.$t/Genelist.rnk | sort -k2,2n | head -n1000 | cut -f1 | sort | uniq > output/DAVID/TAD_W.$t/Genes.txt
	cat output/GSEA/RankTads.$t/Genelist.rnk | cut -f1 | sort | uniq > output/DAVID/TAD_W.$t/Background.txt
	
	Genes="output/DAVID/TAD_W.$t/Genes.txt"
	Background="output/DAVID/TAD_W.$t/Background.txt"
	Idtype="ENTREZ_GENE_ID"
	Output="output/DAVID/TAD_W.$t/Output.txt"
	LogFile="output/DAVID/TAD_W.$t/Log.txt"
	perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR
done

############################
### Background: DE Genes ###
############################

tadSets=`ls output/Correlation/Step3/ | grep TADsets | sed 's/.txt//' | grep Negat | grep -v Names | grep -v ENSG`
cat output/DAVID/Background.TADs.txt | grep -w -F -f <(cat output/SAM/Quantile_3.IDH_vs_DNMT3A.SAM.*Sign.txt | cut -f3 | grep -v Gene | sed 's/"//g' | cut -d '|' -f2 | sort | uniq) > output/DAVID/Background.DE.txt

for t in $tadSets
do
	echo DE.$t
	mkdir output/DAVID/DE.$t
	cat output/Correlation/Step3/$t.txt | cut -d '|' -f2 | sort | uniq | grep -w -F -f output/DAVID/Background.DE.txt > output/DAVID/DE.$t/Genes.txt
	
	Genes="output/DAVID/DE.$t/Genes.txt"
	Background="output/DAVID/Background.DE.txt"
	Idtype="ENTREZ_GENE_ID"
	Output="output/DAVID/DE.$t/Output.txt"
	LogFile="output/DAVID/DE.$t/Log.txt"
	perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<10' | sed 's/__^__/ /g' >> $OutputFDR
done







