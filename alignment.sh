#!/bin/bash

ref=/home/bsb/PublicData/ucsc.hg19.ref/ucsc.hg19.fasta
agilent_v4_71m_bed=/home/bsb/PublicData/Agilent_v4_71m_reduced.bed
mills_indel=/home/bsb/Public_data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
one_k_indel=/home/bsb/Public_data/1000G_phase1.indels.hg19.sites.vcf
dbsnp=/home/bsb/Public_data/dbsnp_138.hg19.vcf

project=${1}
read1=${2}
read2=${3}
sample=$( basename ${read1} .R1.fastq.gz )
dir=/home/bsb/projects/${project}/${sample}
samle_prefix=${dir}/${sample}
sam=${samle_prefix}.sam
bam=${samle_prefix}.view.bam
sorted_bam=${samle_prefix}.sorted.bam
arrg_bam=${samle_prefix}.arrg.bam
md_bam=${samle_prefix}.md.mab
metrix_txt=${samle_prefix}.metrix.txt
RTC_intervals=${samle_prefix}.RTC.intervals
IR_bam=${samle_prefix}.IR.bam
BR_table=${samle_prefix}.BR_table.txt
BQSR_txt=${samle_prefix}.BQSR_table.txt
BQSR_bam=${samle_prefix}.BQSR.bam

mkdir -p ${dir}

/home/bsb/soft/bwa/bwa mem -M -t 8 ${ref} ${read1} ${read2} > ${sam}

/home/bsb/soft/samtools/samtools view -bT ${ref} ${sam} > ${bam}

/home/bsb/soft/samtools/samtools sort -l 9 -O bam ${bam} > ${sortedbam}

java -jar /home/bsb/soft/picard_tools/picard.jar AddOrReplaceReadGroups \
  INPUT=${sorted_bam} \
  OUTPUT=${arrg_bam} \
  SORT_ORDER=coordinate \
  RGID=${sample} \
  RGLB=${sample} \
  RGPL=ILLUMINA \
  RGPU=SureSelectV4 \
  RGSM=${sample} \
  RGCN=NLA \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_RECORDS_IN_RAM=1000000

java -jar /home/bsb/soft/picard_tools/picard.jar MarkDuplicates \
  INPUT=${arrg_bam} \
  OUTPUT=${md_bam} \
  METRICS_FILE=${metrix_txt} \
  ASSUME_SORTED=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R ${ref} \
  -I ${md_bam} \
  -L ${agilent_v4_71m_bed} \
  -known ${mills_indel} \
  -known ${one_k_indel} \
  -o ${RTC_intervals}

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R ${ref} \
  -I ${md_bam} \
  -L ${agilent_v4_71m_bed} \
  -targetIntervals ${RTC_intervals} \
  -known ${mills_indel} \
  -known ${one_k_indel} \
  -o ${IR_bam}

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${ref} \
  -I ${IR_bam} \
  -knownSites ${dbsnp} \
  -knownSites ${mills_indel} \
  -knownSites ${one_k_indel} \
  -o ${BR_table}

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T BaseRecalibrator \
  -R ${ref} \
  -I ${IR_bam} \
  -knownSites ${dbsnp} \
  -knownSites ${mills_indel} \
  -knownSites ${one_k_indel} \
  -BQSR ${BR_table} \
  -o ${BQSR_table}

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T PrintReads \
  -R ${ref} \
  -I ${IR_bam} \
  -BQSR ${BQSR_table} \
  -o ${BQSR_bam}
