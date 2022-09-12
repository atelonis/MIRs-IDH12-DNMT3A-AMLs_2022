#!/bin/bash 

#Juicer step1: from mapping to get interation
bash /opt/juicer/scripts/juicer.sh -g hg19 -d /opt/juicer/work/HiC/HiC1_rmchr/ -y /opt/juicer/restriction_sites/hg19_GATC_GANTC_rmchr.txt -z /opt/juicer/references/Homo_sapiens_assembly19.fasta -p /opt/juicer/references/hg19_rmchr.chrom.sizes -s Arima -t 15 &

zcat /opt/juicer/work/HiC/HiC1_rmchr/fastq/HiC_1_R1.fastq.gz | split -a 3 -l 90000000 -d - /opt/juicer/work/HiC/HiC1_rmchr/splits/HiC_1_R1.fastq &

zcat /opt/juicer/work/HiC/HiC1_rmchr/fastq/HiC_1_R2.fastq.gz | split -a 3 -l 90000000 -d - /opt/juicer/work/HiC/HiC1_rmchr/splits/HiC_1_R2.fastq &

bash /opt/juicer/scripts/juicer.sh -g hg19 -d /opt/juicer/work/HiC/HiC1_rmchr/ -y /opt/juicer/restriction_sites/hg19_GATC_GANTC_rmchr.txt -z /opt/juicer/references/Homo_sapiens_assembly19.fasta -p /opt/juicer/references/hg19_rmchr.chrom.sizes -s Arima -t 15 -S dedup &

bash after_dedup.sh

bash /opt/juicer/scripts/juicer.sh -g hg19 -d /opt/juicer/work/HiC/HiC1_rmchr/ -y /opt/juicer/restriction_sites/hg19_GATC_GANTC_rmchr.txt -z /opt/juicer/references/Homo_sapiens_assembly19.fasta -p /opt/juicer/references/hg19_rmchr.chrom.sizes -s Arima -t 15 -S final &

#Juicer step2: get TAD regions using arrowhead
java -jar /opt/juicer/scripts/juicer_tools.7.0.jar arrowhead -r 40000 -k KR --threads 10 inter_Qin.hic inter_Qin_domains_list



