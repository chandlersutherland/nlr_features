#!/bin/bash
#SBATCH --job-name=methylation_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load bowtie2
source activate e14
module load samtools

#need to export three variables, trim_output directory, which is our input directory, bismark_output directory, and the genome file 
cd $trim_output
echo "beginning bismark on sample ${sample}"

bismark --genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	--non_directional \
	-1 "${sample}"_1_val_1.fq -2 "${sample}"_2_val_2.fq
