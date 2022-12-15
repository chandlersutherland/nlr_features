#!/bin/bash
#SBATCH --job-name=methylation_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=06:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
module load bowtie2
source activate e14
module load samtools

#need to export three variables, trim_output directory, which is our input directory, bismark_output directory, and the genome file 
cd trim_output

#for each sample, write the mates to a csv file 
samples=$(basename -s _1_trimmed.fq  *_1*.fq)
touch mates_1.csv 
touch mates_2.csv 

for sample in $samples 
do 
	echo -n $trim_output/${sample}_1_trimmed.fq, >> mates_1.csv 
	echo -n $trim_output/${sample}_2_trimmed.fq, >> mates_2.csv
done 

MATES1=$(cat mates_1.csv)
MATES2=$(cat mates_2.csv)

bismark --genome $genome \
	--temp_dir $SCRATCH \
	--output_dir $bismark_output \
	-p 4 \
	--non_directional \
	-1 $MATES1 -2 $MATES2