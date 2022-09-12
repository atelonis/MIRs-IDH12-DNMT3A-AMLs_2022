#!/bin/bash 

########################################################################
##### PART 1 | Preparing the genome for gene architecture analyses #####
########################################################################

for hm in Human Mouse
do
	### Find overlaps with Repeats
	for ei in Exon Intron
	do
		echo $ei
		bedtools intersect -a $hm"Genome"/$ei"Ranges.bed" -b $hm"Genome/RepeatMasker.bed" -wa -wb > bed/$hm.BedInters.$ei.Reps.bed
	done
	
	### Find overlaps with evolutionary conservation files
	for ei in Exon Intron
	do
		echo $ei
		bedtools intersect -a $hm"Genome"/$ei"Ranges.bed" -b EvolCons/$hm/chr*.bed -wa -wb > bed/$hm.BedInters.$ei.EvolCons.bed
	done
done

################################################
##### PART 2 | Properties of the cytosines #####
################################################

### For evolutionary conservation
cat output/TAD.pairs.txt | cut -f1 | sort | uniq | tr '.' '\t' | awk '{print $1,$2-102,$2+100,$1"."$2,1000,"+"}' | tr ' ' '\t' | awk '$2>0' > bed/ForTadCEvol.bed
cat output/TAD.pairs.txt | cut -f1 | sort | uniq | tr '.' '\t' | awk '{print $1,$2-3,$2,$1"."$2,1000,"+"}' | tr ' ' '\t' > bed/ForTadSeqExtr.bed
python ExtractSequence.py Human.GRCh37 bed/ForTadSeqExtr.bed tadCytosines.fa

### Repeats
bedtools intersect -a ForTadCEvol.bed -b HumanGenome/RepeatMasker.bed -wa -wb > bed/BedIntersect.TadCs.Repeats.bed 

### evol cons
chros=`ls ~/Genome/Human.GRCh37/EvolCons/ | grep bed`
for c in $chros
do
	echo $c
	bedtools intersect -a ForTadCEvol.bed -b EvolCons/Human/$c -wa -wb > bed/BedIntersect.TadCs.EvolCons_$c
done

