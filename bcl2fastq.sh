#!/bin/bash
#
# Bcl conversion and demultiplexing to .fastq.gz files
# /home/bsb/run_dir/run_YYMMDD/ is an Illumina HiSeq2000 run output dir
# samplesheet.csv is a samplesheet with WES project sample list
perl  /usr/local/bin/configureBclToFastq.pl  \
  --input-dir  /home/bsb/run_dir/run_YYMMDD/Data/Intensities/BaseCalls/  \
  --output-dir  /home/bsb/run_dir/run_YYMMDD/Unaligned_fastq_gz/  \
  --sample-sheet /home/bsb/run_dir/run_YYMMDD/samplesheet.csv \
  --no-eamss  \
  --use-bases-mask  y101,i6n,y101  \
  --ignore-missing-bcl  \
  --ignore-missing-stats  \
  --ignore-missing-control

# merge several .fastq.gz files to single .fastq.gz file
mkdir -p /home/bsb/run_dir/run_YYMMDD/Unaligned_fastq_gz/fastq_gz/
cd /home/bsb/run_dir/run_YYMMDD/Unaligned_fastq_gz/Project_WES
for i in $( ls -d WE* ); do
  cat ${i}/*R1*.fastq.gz > ../fastq_gz/${i}.R1.fastq.gz;
  cat ${i}/*R2*.fastq.gz > ../fastq_gz/${i}.R2.fastq.gz;
done
