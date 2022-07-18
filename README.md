# ATACseq

Scripts and Functions related to ATACseq data analysis

## General Outline
FastQC
Trimmomatic
Bowtie2 (-X 2000 parameter)

### Filter Alignments
Remove alignments that:
  - map to the MT chromosome
  - are duplicate as found by Picard Tools markduplicates
  - have a mapq score < 30
  - are non properly paired as defined by samtools 0x02 flag

Call peaks with MACS2
Differential analysis with DiffBind
