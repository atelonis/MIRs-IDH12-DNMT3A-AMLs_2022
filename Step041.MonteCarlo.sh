#!/bin/bash

####### Preparing
preps="Negat_Overlap Negat_Intermediate Negat_Long"
for i in $preps
do
	cat HumanGenome/GeneRanges.bed | grep -F -f output/Correlation/Step3/TADsets.$i.ENSG.txt | awk '{print $1,$2-2000,$3+2000,$4}' | tr ' ' '\t' > output/Correlation/Step3/TADsets.$i.ENSG.2kbExt.bed
	
	w=`echo $i | cut -d '_' -f2`
	cat HumanGenome/GeneRanges.bed | grep -F -f output/Correlation/Step3/TADsets_W.$w.ENSG.txt | awk '{print $1,$2-2000,$3+2000,$4}' | tr ' ' '\t' > output/Correlation/Step3/TADsets.W_$w.ENSG.2kbExt.bed
done

tadPaths="Negat_Overlap Negat_Intermediate Negat_Long W_Overlap W_Intermediate W_Long"


######## Beck et al.
[ ! -d Beck ] && mkdir Beck
[ ! -d output/BeckMonteCarlo ] && mkdir output/BeckMonteCarlo
[ ! -d output/BeckMonteCarlo/Shuffling ] && mkdir output/BeckMonteCarlo/Shuffling

backGenes="output/Correlation/Step3/TADsets.Background.ENSG.txt"
cat HumanGenome/GeneRanges.bed | grep -F -f $backGenes | awk '{print $1,$2-2000,$3+2000,$4,$5,$6}' | tr ' ' '\t' | sort -k1,1 -k2,2n > output/BeckMonteCarlo/Background.bed
bedtools merge -i output/BeckMonteCarlo/Background.bed > output/BeckMonteCarlo/MergedBackground.bed

beckDir="HumanGenome/Beck.PMID23974199"
Beck=`ls $beckDir | sed 's/_Peaks.txt//'`

for b in $Beck
do
 	echo $b
	bedtools intersect -a Datasets/Dataset.Rownames.bed -b $beckDir/$b"_Peaks.txt" -wa -wb > Beck/BI.Cs.$b.bed
	
	### Consider only the peaks in the background | For shuffling within the background
	bedtools intersect -a $beckDir/$b"_Peaks.txt" -b output/BeckMonteCarlo/MergedBackground.bed > output/BeckMonteCarlo/BI.Background.$b.bed
	
	for i in {1..1000}
	do
		bedtools shuffle -i output/BeckMonteCarlo/BI.Background.$b.bed -incl output/BeckMonteCarlo/Background.bed -g HumanGenome/ChromosomeLengths.txt > output/BeckMonteCarlo/Shuffling/Genespace.$b.$i.bed &
	done
	wait
done

for ci in $tadPaths
do
	for b in $Beck
	do
	echo $ci $b
	bedtools intersect -a output/Correlation/Step3/TADsets.$ci.ENSG.2kbExt.bed -b $beckDir/$b"_Peaks.txt" -wa -wb > Beck/$ci.$b.bed
	for i in {1..1000}
		do
			bedtools intersect -a output/Correlation/Step3/TADsets.$ci.ENSG.2kbExt.bed -b output/BeckMonteCarlo/Shuffling/Genespace.$b.$i.bed -wa -wb > output/BeckMonteCarlo/Shuffling/BI.Genespace.$ci.$b.$i.bed &
		done
	done
	wait
done

###### ENCODE

[ ! -d Encode ] && mkdir Encode
[ ! -d output/EncodeMonteCarlo ] && mkdir output/EncodeMonteCarlo
[ ! -d output/EncodeMonteCarlo/Shuffling ] && mkdir output/EncodeMonteCarlo/Shuffling

backGenes="output/Correlation/Step3/TADsets.Background.ENSG.txt"
cat HumanGenome/GeneRanges.bed | grep -F -f $backGenes | awk '{print $1,$2-2000,$3+2000,$4,$5,$6}' | tr ' ' '\t' | sort -k1,1 -k2,2n > output/EncodeMonteCarlo/Background.bed
bedtools merge -i output/EncodeMonteCarlo/Background.bed > output/EncodeMonteCarlo/MergedBackground.bed

encodeDir="HumanGenome/Encode"
Encode=`ls $encodeDir | sed 's/.bed//'`

for e in $Encode
do
 	echo $e
	bedtools intersect -a Datasets/Dataset.Rownames.bed -b $encodeDir/$e.bed -wa -wb > Encode/BI.Cs.$e.bed
	
	### Consider only the peaks in the background | For shuffling within the background
	bedtools intersect -a $encodeDir/$e.bed -b output/EncodeMonteCarlo/MergedBackground.bed > output/EncodeMonteCarlo/BI.Background.$e.bed
	
	for i in {1..1000}
	do
		bedtools shuffle -i output/EncodeMonteCarlo/BI.Background.$e.bed -incl output/EncodeMonteCarlo/Background.bed -g HumanGenome/ChromosomeLengths.txt > output/EncodeMonteCarlo/Shuffling/Genespace.$e.$i.bed &
	done
	wait
done


for ci in $tadPaths
do
	echo $ci
	for e in $Encode
	do
		bedtools intersect -a output/Correlation/Step3/TADsets.$ci.ENSG.2kbExt.bed -b $encodeDir/$e.bed -wa -wb > Encode/$ci.$e.bed
		for i in {1..1000}
		do
			echo $ci $e $i
			bedtools intersect -a output/Correlation/Step3/TADsets.$ci.ENSG.2kbExt.bed -b output/EncodeMonteCarlo/Shuffling/Genespace.$e.$i.bed -wa -wb > output/EncodeMonteCarlo/Shuffling/BI.Genespace.$ci.$e.$i.bed &
		done
		wait
	done
done


