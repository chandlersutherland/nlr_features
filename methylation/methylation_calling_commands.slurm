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

#pool two biological reps to increase genomic coverage, leaving two biological reps 
#cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281085_1.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281086_1.fastq  > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep1_1.fastq
#cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281085_2.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281086_2.fastq > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep1_2.fastq
#cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281087_1.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281088_1.fastq  > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep2_1.fastq
#cat /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281087_2.fastq /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams/SRR17281088_2.fastq > /global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled/biorep2_2.fastq

#define samples 
samples=("biorep1" "biorep2")

#trim the two pooled biological reps 
trim_input="/global/scratch/users/chandlersutherland/e14/bs_fastq_files/williams_pooled"
trim_output="/global/scratch/users/chandlersutherland/e14/trim_williams/pooled"

#sbatch --job-name=bio_reps.trim_galore --export=trim_input=$trim_input,trim_output=$trim_output,samples=$samples -A co_minium /global/home/users/chandlersutherland/nlr_features/methylation/trim_galore.sh

#run bismark alignment 
genome="/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly"
bismark_output="/global/scratch/users/chandlersutherland/e14/bismark/bam_files/pooled"

sbatch --job-name=bio_reps.bismark --export=trim_output=$trim_output,genome=$genome,bismark_output=$bismark_output,samples=$samples -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/bismark.bash

#check read coverage of NLRs, generating a clean_coverage.tsv file for each .sam file in the coverage_input directory  
#untested 12/21 
coverage_input=$bismark_output
sbatch --job-name=bio_reps.nlr_coverage --export=coverage_input=$coverage_input -A co_minium /global/home/users/chandlersutherland/nlr_features/methylation/samtools_coverage_runner.bash

