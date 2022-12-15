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

#pool two biological reps to increase genomic coverage, leaving two biological reps 
cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281085_1.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281086_1.fastq  > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep1_1.fastq
cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281085_2.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281086_2.fastq > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep1_2.fastq
cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281087_1.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281088_1.fastq  > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep2_1.fastq
cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281087_2.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281088_2.fastq > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep2_2.fastq

#trim the two pooled biological reps 
trim_input="/global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled"
trim_output="/global/scratch/users/chandlersutherland/e14/trim_williams/pooled"

sbatch --job-name=bio_reps.trim_galore --export=trim_input=$trim_input,trim_output=$trim_output /global/home/users/chandlersutherland/nlr_features/methylation/trim_galore.sh
