#!/bin/bash 

[ ! -d TCGA/bed ] && mkdir TCGA/bed
[ ! -d TCGA/output ] && mkdir TCGA/output
[ ! -d TCGA/output/Correlation ] && mkdir TCGA/output/Correlation

cat TCGA/Data/LAML.Methylation.Level_3/jhu-usc.edu_LAML.HumanMethylation450.Level_3.2.3.0/jhu-usc.edu_LAML.HumanMethylation450.2.lvl-3.TCGA-AB-2971-03A-01D-0741-05.txt | awk 'NF==5' | grep -v Composite | grep -v Hybridization | awk '{print "chr"$4,$5-1,$5,$1}' | tr ' ' '\t' | grep -v chrNA > cg.bed

bedtools intersect -a cg.bed -b <(cat TADs.bed TADs.Not.bed) -wa -wb > TCGA/bed/BedInters.cg.TADs.bed
bedtools intersect -a cg.bed -b <(cat HumanGenome/Yng*Enh_peaks.bed) -wa -wb > TCGA/bed/BedIntersect.Cs.Enhancers.bed

### For evolutionary conservation
cat cg.bed | grep -w -F -f <(cat output/TAD.pairs.txt | cut -d '.' -f1 | sort | uniq) > TCGA/bed/ForTadCEvol.bed

### evol cons
chros=`ls HumanGenome/EvolCons/ | grep bed`
for c in $chros
do
	echo $c
	bedtools intersect -a TCGA/bed/ForTadCEvol.bed -b HumanGenome/EvolCons/$c -wa -wb > TCGA/bed/BedIntersect.TadCs.EvolCons_$c
done


