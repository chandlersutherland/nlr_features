#!/bin/bash
#SBATCH --job-name=bismark_deduplicate_n_extract
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
module load python 
source activate e14

mkdir -p $deduplicate_input/deduplicate

deduplicate_bismark -p \
	--output_dir  $deduplicate_input/deduplicate \
	$dedpulicate_input/${sample}_1_val_1_bismark_bt2_pe.bam
echo "finished deduplication of ${sample}" 
	
bismark_methylation_extractor -p \
	--output $extraction_output \
	--ignore_r2 2 \
	--comprehensive \
	--bedGraph \
	--CX \
	--parallel $SLURM_NTASKS \
	$deduplicate_input/deduplicate/${sample}*.bam
		
echo "finished extraction of ${sample}" 
