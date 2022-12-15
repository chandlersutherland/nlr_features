#!/bin/bash
#SBATCH --job-name=methylation_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load cutadapt
module load fastqc 

#set the trim directory, output directory, and input directory, containing the untrimmed fastq files
TRIM_DIR=/global/home/users/chandlersutherland/programs/TrimGalore-0.6.6
#trim_output=/global/scratch/users/chandlersutherland/e14/trim_williams/pooled
#trim_input=/global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled

#run trim galore in default mode 
for f in $trim_input/*
do 
   BASENAME=$(basename $f .fastq)
   $TRIM_DIR/trim_galore -o $trim_output --fastqc --illumina $f
   echo "finished trimming ${BASENAME}"
done 
