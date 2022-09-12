#!/bin/bash 

[ ! -d output/GSEA/CpGi ] && mkdir output/GSEA/CpGi
tads=`ls output/Correlation/Step3/CpGi/ | grep RankTads | sed 's/.txt//'`

for t in $tads
do
	echo $t
	
	[ ! -d output/GSEA/CpGi/$t ] && mkdir output/GSEA/CpGi/$t
	cat output/Correlation/Step3/CpGi/$t.txt | awk '{print $1,$2*100}' | tr ' ' '\t' > output/GSEA/CpGi/$t/Genelist.rnk
	Scripts/GSEA.sh selKEGG_Entrez output/GSEA/CpGi/$t/Genelist.rnk output/GSEA/CpGi/$t/output > output/GSEA/CpGi/$t/Log.txt
done


