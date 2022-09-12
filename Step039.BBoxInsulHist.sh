#!/bin/bash 

#############
### B-Box ###
#############

cat HumanGenome/RepeatMasker.bed | grep MIR | awk '{print $1,$2-150,$3+150,$4,$5,$6}' | tr ' ' '\t' > HumanGenome/MIRs.Ext.bed

python Scripts/ExtractSequence.py Human.GRCh37 HumanGenome/MIRs.Ext.bed HumanGneome/MIRs.Ext.fa


##################
### Insulators ###
##################

bedtools intersect -a <(cat PMID26216945.MirInsulators/liftover_fromSupp_to_hg19.bed | awk '{print $1,$2-100,$3+100}' | tr ' ' '\t') -b HumanGenome/RepeatMasker.bed -wa -wb > PMID26216945.MirInsulators/BI.Insul.Reps.bed


#####################
### Histone marks ###
#####################

for h in H3K27ac H3K27me3 H3K4me1 H3K4me3
do
	bedtools intersect -a PMID31085557.Adelman/$h.bed -b HumanGenome/RepeatMasker.bed -wa -wb > PMID31085557.Adelman/BI.$h.Repeats.bed &
done


#########################################
### For TFs binding the specific MIRs ###
#########################################

TFs=`ls HumanGenome/Encode/ | sed 's/.bed//'`
repMask="HumanGenome/RepeatMasker.bed"

for tf in $TFs
do
	echo $tf
	bedtools intersect -a HumanGenome/Encode/$tf.bed -b $repMask -wa -wb > HumanGenome/BI.$tf.RM.bed
done


