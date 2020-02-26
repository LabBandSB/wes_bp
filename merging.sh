#!/bin/bash

ref=/home/bsb/PublicData/ucsc.hg19.ref/ucsc.hg19.fasta
agilent_v4_71m_bed=/home/bsb/PublicData/Agilent_v4_71m_reduced.bed
mills_indel=/home/bsb/Public_data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
one_k_indel=/home/bsb/Public_data/1000G_phase1.indels.hg19.sites.vcf
dbsnp=/home/bsb/Public_data/dbsnp_138.hg19.vcf
hapmap=/home/bsb/Public_data/hapmap_3.3.hg19.sites.vcf
omni=/home/bsb/Public_data/1000G_omni2.5.hg19.sites.vcf
one_k_snp=/home/bsb/Public_data/1000G_phase1.snps.high_confidence.hg19.sites.vcf

project=${1}

cd ${project}
mkdir -p /home/bsb/projects/${project}/MS_vcf
project_suffix=/home/bsb/projects/${project}/MS_vcf/${project}
ms_vcf=${project_suffix}.vcf

input_bams=''
for i in $( ls * -d ); do input_bams="${input_bams} -I ${i}/*.BQSR.bam "; done

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T HaplotypeCaller \
  -R ${ref} \
  ${input_bams} \
  -D ${dbsnp} \
  --genotyping_mode DISCOVERY \
  -stand_call_conf 30 \
  -o ${ms_vcf} \
  -nct 8

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${ref} \
  -input ${ms_vcf} \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
  -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${omni} \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
  -an DP \
  -an QD \
  -an FS \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode SNP \
  -tranche 100.0 \
  -tranche 99.9 \
  -tranche 99.0 \
  -tranche 90.0 \
  -recalFile ${project_suffix}.MSHC.VR_SNP.recal \
  -tranchesFile ${project_suffix}.MSHC.VR_SNP.tranches \
  -rscriptFile ${project_suffix}.MSHC.VR_SNP_plots.R

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -R ${ref} \
  -input ${ms_vcf} \
  -mode SNP \
  --ts_filter_level 99.0 \
  -recalFile ${project_suffix}.MSHC.VR_SNP.recal \
  -tranchesFile ${project_suffix}.MSHC.VR_SNP.tranches \
  -o ${project_suffix}.MSHC.VQSR_AR_SNP.vcf

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ${ref} \
-input ${project_suffix}.MSHC.VQSR_AR_SNP.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 ${mills_indel} \
-an DP \
-an FS \
-an MQRankSum \
-an ReadPosRankSum \
-mode INDEL \
-tranche 100.0 \
-tranche 99.9 \
-tranche 99.0 \
-tranche 90.0 \
--maxGaussians 4 \
-recalFile ${project_suffix}.MSHC.VR_INDEL.recal \
-tranchesFile ${project_suffix}.MSHC.VR_INDEL.tranches \
-rscriptFile ${project_suffix}.MSHC.VR_INDEL_plots.R

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R ${ref} \
-input ${project_suffix}.MSHC.VQSR_AR_SNP.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile ${project_suffix}.MSHC.VR_INDEL.recal \
-tranchesFile ${project_suffix}.MSHC.VR_INDEL.tranches \
-o ${project_suffix}.MSHC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref} \
-V ${project_suffix}.MSHC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf \
-selectType SNP \
-o ${project_suffix}.MSHC.VQSR.raw_snp.vcf

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${ref} \
-V ${project_suffix}.MSHC.VQSR.raw_snp.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "SNP_FAIL" \
-o ${project_suffix}.MSHC.VQSR.fil_snp.vcf

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T SelectVariants \
-R ${ref} \
-V ${project_suffix}.MSHC.VQSR_AR_SNP.VQSR_AR_INDEL.vcf \
-selectType INDEL \
-o ${project_suffix}.MSHC.VQSR.raw_indel.vcf

java -jar /home/bsb/soft/gatk/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${ref} \
-V ${project_suffix}.MSHC.VQSR.raw_indel.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "INDEL_FAIL" \
-o ${project_suffix}.MSHC.VQSR.fil_indel.vcf

/home/bsb/soft/vcftools/bin/vcf-concat ${project_suffix}.MSHC.VQSR.fil_snp.vcf ${project_suffix}.MSHC.VQSR.fil_indel.vcf > ${project_suffix}.MSHC.VQSR.fil.concat.vcf

cat ${project_suffix}.MSHC.VQSR.fil.concat.vcf | /home/bsb/soft/vcftools/bin/vcf-sort > ${project_suffix}.MSHC.VQSR.fil.concat.sorted.vcf
