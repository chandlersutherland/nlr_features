#!/bin/bash
#SBATCH --job-name=vcf_thin
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00 
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load vcftools
cd /global/scratch/users/chandlersutherland/e14/popgen/vcf_1001

#vcftools --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf --keep nlrome_IDs.txt --recode --recode-INFO-all --out nlrome_partial 
bed='/global/scratch/users/chandlersutherland/Athaliana/Atha_cds.bed'
#vcftools --vcf nlrome_partial.recode.vcf --chr 1 --recode --bed $bed nlrome_partial_1_cds 

#vcftools --vcf nlrome_partial.recode.vcf --chr 2 --recode --bed $bed --out nlrome_partial_2_cds
#vcftools --vcf nlrome_partial.recode.vcf --chr 3 --recode --bed $bed --out nlrome_partial_3_cds
#vcftools --vcf nlrome_partial.recode.vcf --chr 4 --recode --bed $bed --out nlrome_partial_4_cds
#vcftools --vcf nlrome_partial.recode.vcf --chr 5 --recode --bed $bed --out nlrome_partial_5_cds

#repeat with gene bed file
bed_gene='/global/scratch/users/chandlersutherland/Athaliana/Atha_genes.recode.bed'
vcftools --vcf nlrome_partial.recode.vcf --chr 1 --recode --bed $bed_gene --out nlrome_partial_1_gene
vcftools --vcf nlrome_partial.recode.vcf --chr 2 --recode --bed $bed_gene --out nlrome_partial_2_gene
vcftools --vcf nlrome_partial.recode.vcf --chr 3 --recode --bed $bed_gene --out nlrome_partial_3_gene
vcftools --vcf nlrome_partial.recode.vcf --chr 4 --recode --bed $bed_gene --out nlrome_partial_4_gene
vcftools --vcf nlrome_partial.recode.vcf --chr 5 --recode --bed $bed_gene --out nlrome_partial_5_gene
