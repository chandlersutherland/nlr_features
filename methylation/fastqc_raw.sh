#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load fastqc 

SCRATCH_DIR=/global/scratch/users/chandlersutherland/e14/
OUTPUT_DIR=$SCRATCH_DIR/fastqc_raw/
cd $SCRATCH_DIR
FILES=$(find . -type f -print)

fastqc -o $OUTPUT_DIR $FILES 
