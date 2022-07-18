#!/bin/bash

#SBATCH -p exacloud                # partition (queue)
#SBATCH -N 1                       # number of nodes
#SBATCH -n 12                      # number of cores
#SBATCH --mem 99000                # memory pool for all cores
#SBATCH -t 0-24:00                 # time (D-HH:MM)
#SBATCH -o rmdup_%A_%a.out         # Standard output
#SBATCH -e rmdup_%A_%a.err         # Standard error
#SBATCH --array=1-9                # sets number of jobs in array


source activate dnaseq_2022

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

dir=$(pwd)

# set infile directory path
inpath=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/noMT

# set outdir path
OUTDIR=/home/groups/hoolock2/u0/bd/Projects/ECP60/data/bowtie2/noMT_dedup

# path to picard.jar file
picard_path=/home/groups/hoolock2/u0/bd/miniconda3/envs/dnaseq_2022/share/picard-2.18.29-0

# set path to arguments file
args_file=$dir/args_file_dedup

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`
# set variables based on the output of the line ARGS
IFS=: read INFILE BASENAME <<< $ARGS

# run picard mark duplicates on each pair
java -jar $picard_path/picard.jar MarkDuplicates \
I=$inpath/$INFILE \
O=$OUTDIR/$BASENAME.dedup.bam \
M=$OUTDIR/$BASENAME.dup_metrics.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
