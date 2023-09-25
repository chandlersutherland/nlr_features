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

#vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --keep nlrome_IDs.txt --recode --recode-INFO-all --out nlrome_partial 
bed='/global/scratch/users/chandlersutherland/Athaliana/Atha_cds.bed'
gvcf_input='1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz'
#vcftools --gzvcf $gvcf_input --chr 1 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_1_cds 
#vcftools --gzvcf $gvcf_input --chr 2 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_2_cds
vcftools --gzvcf $gvcf_input --chr 3 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_3_cds
vcftools --gzvcf $gvcf_input --chr 4 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_4_cds
vcftools --gzvcf $gvcf_input --chr 5 --keep nlrome_IDs.txt --recode --bed $bed --out nlrome_invar_5_cds

#repeat with gene bed file
#bed_gene='/global/scratch/users/chandlersutherland/Athaliana/Atha_genes.recode.bed'
#vcftools --vcf nlrome_partial.recode.vcf --chr 1 --recode --bed $bed_gene --out nlrome_partial_1_gene
#vcftools --vcf nlrome_partial.recode.vcf --chr 2 --recode --bed $bed_gene --out nlrome_partial_2_gene
#vcftools --vcf nlrome_partial.recode.vcf --chr 3 --recode --bed $bed_gene --out nlrome_partial_3_gene
#vcftools --vcf nlrome_partial.recode.vcf --chr 4 --recode --bed $bed_gene --out nlrome_partial_4_gene
#vcftools --vcf nlrome_partial.recode.vcf --chr 5 --recode --bed $bed_gene --out nlrome_partial_5_gene
