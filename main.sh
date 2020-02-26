#!/bin/bash

project=KAZWES

bash bcl2fastq.sh

cd /home/bsb/run_dir/run_YYMMDD/Unaligned_fastq_gz/fastq_gz/

for read1 in $( ls *R1* ); do
  sample=$( basename ${read1} .R1.fastq.gz ) ;
  read2=${sample}.R2.fastq.gz ;
  bash alignment.sh ${project} ${read1} ${read2} ;
done

bash merging.sh ${project}
