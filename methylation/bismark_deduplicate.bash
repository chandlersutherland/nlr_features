#!/bin/bash
#SBATCH --job-name=bismark_deduplicate
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load bowtie2
module load samtools

BISMARK=/global/home/users/chandlersutherland/programs/Bismark-0.23.0
ARAPORT11=/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly
OUTPUT_DIR=/global/scratch/users/chandlersutherland/e14/bismark
BISULFITE='SRR17281088 SRR17281087 SRR17281086 SRR17281085'

cd $BISMARK 

for f in $BISULFITE
do 
	./deduplicate_bismark -p \
		--output_dir  $SCRATCH/e14/bismark/deduplicate_bismark/ \
		/global/scratch/users/chandlersutherland/e14/bismark/bam_files/${f}_1_val_1_bismark_bt2_pe.bam
done 
