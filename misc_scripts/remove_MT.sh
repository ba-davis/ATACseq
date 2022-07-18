#!/bin/bash

#conda activate dnaseq_2022

dir=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/sorted_bams
outdir=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/noMT

for i in $dir/*.sorted.bam; do samtools idxstats $i | cut -f 1 | grep -v MT | xargs samtools view -b $i > $i.noMT.bam; done
