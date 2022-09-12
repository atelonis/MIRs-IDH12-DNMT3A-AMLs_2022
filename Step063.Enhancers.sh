#!/bin/bash 

[ ! -d output/Enhancers ] && mkdir output/Enhancers

echo AllBack allDM enhBack enhDM > output/Enhancers/Table.txt 
for am in MIR Alu L1
do
	for db in DiffMeth Background
	do
		bedtools intersect -a <(cat Human/Genome/Yng_*Enh_peaks.bed) -b <(cat output/GEAR/SAM_$am.$db.txt | tr '.' '\t' | tr '-' '\t') -wa -wb > output/Enhancers/BI.$am.$db.bed
	done
	allBack=`cat output/GEAR/SAM_$am.Background.txt | wc -l`
	allDM=`cat output/GEAR/SAM_$am.DiffMeth.txt | wc -l`
	enhBack=`cat output/Enhancers/BI.$am.Background.bed | cut -f 4-6 | sort | uniq | wc -l`
	enhDM=`cat output/Enhancers/BI.$am.DiffMeth.bed | cut -f 4-6 | sort | uniq | wc -l`
	echo $am $allBack $allDM $enhBack $enhDM >> output/Enhancers/Table.txt
done

