#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:10:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
source activate e14

#The sjdb overhang depends on the read length. Here we have 50 reads, so set to 49 
GENOME_DIR=/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly/STAR_genome_williams
FASTA_FILE=/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly/Athaliana_447_TAIR10.fa
GTF_FILE=/global/scratch/users/chandlersutherland/Athaliana/Araport11_GTF_genes_transposons.current.gtf

STAR --runThreadN $SLURM_NTASKS \
	 --runMode genomeGenerate \
	 --genomeDir $GENOME_DIR_WILLIAMS \
	 --genomeFastaFiles $FASTA_FILE \
	 --genomeSAindexNbases 12 \
	 --sjdbGTFfile $GTF_FILE  \
	 --sjdbOverhang 49

echo 'genome_dir_williams made' 