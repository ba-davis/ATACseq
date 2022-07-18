#!/bin/bash

# conda activate dnaseq_2022

dir=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/noMT_dedup

# remove multimapped reads (i.e. those with mapq < 30)
for i in $dir/*.bam; do samtools view -h -q 30 $i > $i.mapq30.bam; done
