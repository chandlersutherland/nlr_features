#!/bin/bash
#SBATCH --job-name=vcf_thin
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=15:00:00 
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load vcftools
cd /global/scratch/users/chandlersutherland/e14/popgen/vcf_1001_full

#goal: subset the giant all sites VCF by the NLRome accessions and by the CDS sequences
#create one file per chromosome to run in parallel 

bed='/global/scratch/users/chandlersutherland/Athaliana/Atha_cds.bed'
gvcf_input='1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz'
vcftools --gzvcf $gvcf_input --chr 1 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_1_cds 
vcftools --gzvcf $gvcf_input --chr 2 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_2_cds
vcftools --gzvcf $gvcf_input --chr 3 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_3_cds
vcftools --gzvcf $gvcf_input --chr 4 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_4_cds
vcftools --gzvcf $gvcf_input --chr 5 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_5_cds

