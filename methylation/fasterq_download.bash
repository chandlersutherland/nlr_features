#!/bin/bash
#SBATCH --job-name=fastq_download
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=10:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

## download the sra files for epigenomes and rnaseq 

module load python
module load sra-tools

#download 1001 epigenomes methylc and rna 
SCRATCH_DIR=/global/scratch/users/chandlersutherland/e14/

#echo "beginning ecker" 

#fasterq-dump --threads $SLURM_NTASKS -O $SCRATCH_DIR/bs_fastq_files/ecker -t $SCRATCH_DIR -p SRR771698
#fasterq-dump --threads $SLURM_NTASKS -O $SCRATCH_DIR/rna_fastq_files/ecker -t $SCRATCH_DIR -p SRR3465232

echo "ecker finished, beginning williams" 

#download williams data 
bisulfite='SRR17281088 SRR17281087 SRR17281086 SRR17281085'
rna='SRR17281236 SRR17281235 SRR17281234 SRR17281233'

for i in $bisulfite
	do 
	fasterq-dump --threads $SLURM_NTASKS -O $SCRATCH_DIR/bs_fastq_files/williams -t $SCRATCH_DIR -p $i
done 

echo "finished williams bisulfite, beginning RNA"

for i in $rna
	do 
	fasterq-dump --threads $SLURM_NTASKS -O $SCRATCH_DIR/bs_fastq_files/williams -t $SCRATCH_DIR -p $i
done
