#!/bin/bash 

for i in bed output
do
	[ ! -d $i ] && mkdir $i
done

# Merge overlaping TADs
bedtools merge -i <(cat HumanGenome/FromQin_TADs.bed | sort -k1,1 -k2,2n) | awk '{print $1,$2,$3,"TAD_"$1"_"$2"_"$3}' | tr ' ' '\t' > TADs.bed
# Areas of the genome that are not covered by TADs, consider them as pseudo-TADs
bedtools subtract -a <(cat HumanGenome/ChromosomeLengths.txt | awk '{print $1,1,$2}' | tr ' ' '\t') -b TADs.bed | awk '{print $1,$2,$3,"TADnot_"$1"_"$2"_"$3}' | tr ' ' '\t' > TADs.Not.bed

bedtools intersect -a TADs.bed -b Datasets/Dataset.Rownames.bed -wa -wb > output/TAD.Cs.bed
bedtools intersect -a TADs.bed -b HumanGenome/GeneRanges.bed -wa -wb > output/TAD.genes.bed

bedtools intersect -a TADs.Not.bed -b Datasets/Dataset.Rownames.bed -wa -wb > output/TADnot.Cs.bed
bedtools intersect -a TADs.Not.bed -b HumanGenome/GeneRanges.bed -wa -wb > output/TADnot.genes.bed

bedtools intersect -a <(cat TADs.bed TADs.Not.bed) -b HumanGenome/RepeatMasker.bed -wa -wb > TADs.Repeats.bed

### Below, merge repeats of the same family
allTADs=`cat TADs.bed TADs.Not.bed | cut -f 4`
printf '' > TADs.MergedRepeats.bed
for t in $allTADs
do
	echo $t
	cat TADs.Repeats.bed | grep $t | cut -d '|' -f1 | sed 's/\//__/' > temp.bed
	allReps=`cat temp.bed | cut -f 8 | sort | uniq | grep -v "?"`
	for r in $allReps
	do
		bedtools merge -i <(cat temp.bed | awk -v r=$r '$8==r' | awk '{print $5,$6,$6}' | tr ' ' '\t' | sort -k1,1 -k2,2n) | awk -v r=$r -v t=$t '{print $1,$2,$3,t,r}' | tr ' ' '\t' | sed 's/__/\//' >> TADs.MergedRepeats.bed
	done
done

### Enhancers from Adelman et al., PMID31085557
bedtools intersect -a Datasets/Dataset.Rownames.bed -b <(cat HumanGenome/Yng*Enh_peaks.bed) -wa -wb > bed/BedIntersect.Cs.Enhancers.bed

### CpG islands
bedtools intersect -a Datasets/Dataset.Rownames.bed -b bed/cpgi.hg19.bed.txt -wa -wb > bed/BedIntersect.Cs.CpGi.bed

### Repeats
bedtools intersect -a Datasets/Dataset.Rownames.bed -b HumanGenome/RepeatMasker.bed -wa -wb > bed/BI.Cs.Reps.bed


