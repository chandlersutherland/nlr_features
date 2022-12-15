#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load cutadapt
module load fastqc 

TRIM_DIR=/global/home/users/chandlersutherland/programs/TrimGalore-0.6.6
OUT_DIR=/global/scratch/users/chandlersutherland/e14/trim_williams/rna/
rna='SRR17281236 SRR17281235 SRR17281234 SRR17281233'

cd /global/scratch/users/chandlersutherland/e14/rna_fastq_files/williams/

for f in $rna
do 
   $TRIM_DIR/trim_galore -o $OUT_DIR --fastqc --illumina "${f}".fastq 
done 
