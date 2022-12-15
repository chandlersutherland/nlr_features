#!/bin/bash
#SBATCH --job-name=bismark_extractor
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load bowtie2
module load samtools

BISMARK=/global/home/users/chandlersutherland/programs/Bismark-0.23.0
OUTPUT_DIR=/global/scratch/users/chandlersutherland/e14/bismark/extraction

#already ran SRR17281088 in interactive to test, so don't need it here 
BISULFITE='SRR17281087 SRR17281086 SRR17281085'

cd $BISMARK 

#williams time 
for f in $BISULFITE
do 
	./bismark_methylation_extractor -p \
		--output  $OUTPUT_DIR \
		--ignore_r2 2 \
		--comprehensive \
		--parallel $SLURM_NTASKS \
		/global/scratch/users/chandlersutherland/e14/bismark/deduplicate_bismark/${f}_1_val_1_bismark_bt2_pe.deduplicated.bam
	echo "finished" $f 
done 

#ecker hehe. Keeping comprehensive for now  
./bismark_methylation_extractor -s \
		--output  $OUTPUT_DIR \
		--comprehensive \
		--parallel $SLURM_NTASKS \
		/global/scratch/users/chandlersutherland/e14/bismark/deduplicate_bismark/SRR771698_bismark_bt2.deduplicated.bam
