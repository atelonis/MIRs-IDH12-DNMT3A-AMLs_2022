#!/bin/bash 

mkdir GSE60055/LogFiles
mkdir GSE60055/Datasets
mkdir GSE60055/StarMapping

srr_column=`cat GSE60055/SraRunTable.txt | head -n1 | tr '\t' '\n' | cat -n | awk '$2=="Run"' | awk '{print $1}'`
samples=`cat GSE60055/SraRunTable.txt | tr ' ' '_' | cut -f $srr_column | grep SRR`
for i in $samples
do
	echo $i
	/home/axt5207/Scripts/sratoolkit.2.9.6-1-centos_linux64/bin/prefetch $i > GSE60055/LogFiles/$i.txt
done

mv ~/ncbi/public/sra/* GSE60055/Datasets/

# paired-end RNA-seq
srr_column=`cat GSE60055/SraRunTable.txt | head -n1 | tr '\t' '\n' | cat -n | awk '$2=="Run"' | awk '{print $1}'`
samples=`cat GSE60055/SraRunTable.txt | tr ' ' '_' | cut -f $srr_column | grep SRR`

cd GSE60055/Datasets
for i in $samples
do
	echo $i
	/home/axt5207/Scripts/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump -I --split-files $i.sra
done
cd ../../

for s in $samples
do
	cGSM=`cat GSE60055/SraRunTable.txt | grep $s | cut -f 17`
	cName=`cat GSE60055/GSM_Accession.txt | grep $cGSM | cut -f2`
	echo $s $cName
	STAR --runThreadN 30 --genomeDir MouseGenome/for_star/ --readFilesIn GSE60055/Datasets/$s"_1.fastq" GSE60055/Datasets/$s"_2.fastq" --outFileNamePrefix GSE60055/StarMapping/$cName. --quantMode TranscriptomeSAM &> GSE60055/LogFiles/STAR.$cName.txt
	rsem-calculate-expression -p 30 --paired-end --bam GSE60055/StarMapping/$cName.Aligned.toTranscriptome.out.bam MouseGenome/for_rsem/mouse StarMapping/$cName &> GSE60055/LogFiles/RSEM.$cName.txt
done


