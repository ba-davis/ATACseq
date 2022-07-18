#!/bin/bash


dir=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/noMT_dedup_mapq30

# retain properly paired reads only (defined by samtools 0x02 flag)
for i in $dir/*.bam; do samtools view -b -f 0x02 $i > $i.proper.bam; done
