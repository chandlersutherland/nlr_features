#!/bin/bash
#SBATCH --job-name=process_polyester
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:40:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#First, align to 

#define directory variables
INPUT=/global/scratch/users/chandlersutherland/e14/polyester/primary_simulated_reads_1006
SCRATCH_DIR=/global/scratch/users/chandlersutherland/e14/

cd $INPUT
FILES=$(for f in *.fasta; do echo $f; done)

#step 1: STAR 
GENOME_DIR_WILLIAMS=/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly/STAR_genome_williams
STAR_OUTPUT=$INPUT/STAR
mkdir $STAR_OUTPUT

#define star function, which takes in a file and writes the STAR output to the STAR_OUTPUT directory 
STAR_RUN (){ 
    GENOME_DIR_WILLIAMS=/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly/STAR_genome_williams
    STAR_OUTPUT=/global/scratch/users/chandlersutherland/e14/polyester/primary_simulated_reads_1006/STAR
    BASENAME=$(basename $1 .fasta)
    STAR --runThreadN $SLURM_NTASKS --genomeDir $GENOME_DIR_WILLIAMS --outFileNamePrefix "${STAR_OUTPUT}"/"${BASENAME}"_ --readFilesIn $1
    echo "finished" $BASENAME
}

#HTCOUNT_RUN, which takes in a .sam file and produces a NLR-only htseq count file and then a all genome htseq-count file 
COUNTS_OUTPUT=$INPUT/htseq_counts


HTCOUNT_RUN () {
    COUNTS_OUTPUT=/global/scratch/users/chandlersutherland/e14/polyester/primary_simulated_reads_1006/htseq_counts
    BASENAME=$(basename $1 .sam)
    htseq-count -r pos -s yes -c $COUNTS_OUTPUT/${BASENAME}_NLRs.tsv $1 /global/scratch/users/chandlersutherland/Athaliana/GTFs/all_NLRs.gtf
    htseq-count -r pos -s yes -c $COUNTS_OUTPUT/${BASENAME}.tsv $1 /global/scratch/users/chandlersutherland/Athaliana/GTFs/Araport11_GTF_genes_transposons.current.gtf
    echo 'finished' ${BASENAME}
}

#apply STAR_RUN to the input files
export -f STAR_RUN
export -f HTCOUNT_RUN

#parallel STAR_RUN ::: $FILES 
#echo "finished STAR, moving on to htseq_counts" 

parallel HTCOUNT_RUN ::: $STAR_OUTPUT/*.sam
echo "finished HTCOUNT, ready for processing" 
