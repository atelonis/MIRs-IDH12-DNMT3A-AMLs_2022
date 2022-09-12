#!/bin/bash 

[ ! -d GSE60055/output/DAVID ] && mkdir GSE60055/output/DAVID

for mg in MEP
do
for i in Double KO3a M140
do
for pn in Posit Negat
do
	echo $mg $i $pn
	[ ! -d GSE60055/output/DAVID/$mg.$i.$pn ] && mkdir GSE60055/output/DAVID/$mg.$i.$pn
	
	cat GSE60055/output/SAM/$mg.$i"_vs_WT.SAM.IncludedInSAM.txt" | cut -d '|' -f1 | sort | uniq > GSE60055/output/DAVID/$mg.$i.$pn/Background.txt
	cat GSE60055/output/SAM/$mg.$i"_vs_WT.SAM."$pn"Sign.txt" | cut -f 3 | sed 's/"//g' | cut -d '|' -f1 | grep ENSMUSG > GSE60055/output/DAVID/$mg.$i.$pn/Genes.txt
	
	Genes="GSE60055/output/DAVID/$mg.$i.$pn/Genes.txt"
	Background="GSE60055/output/DAVID/$mg.$i.$pn/Background.txt"
	Idtype="ENSEMBL_GENE_ID"
	Output="GSE60055/output/DAVID/$mg.$i.$pn/Output.txt"
	LogFile="GSE60055/output/DAVID/$mg.$i.$pn/Log.txt"
	perl ~/Scripts/David.pl $Genes $Background $Idtype $Output > $LogFile
	
	OutputFDR=`echo $Output | sed 's/Output.txt/Output.FDRfilter.txt/'`
	cat $Output | head -n1 > $OutputFDR
	cat $Output | sed 's/ /__^__/g' | awk '$13<5' | sed 's/__^__/ /g' >> $OutputFDR
done
done
done

