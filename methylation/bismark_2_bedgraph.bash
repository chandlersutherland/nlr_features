#!/bin/bash
#SBATCH --job-name=bismark_extractor
#SBATCH --partition=savio2
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=01:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python 
source activate e14
module load samtools

#could not get this to work, no errors but would not save the file! 
BISMARK_BEDGRAPH () {
	basename = $(basename $1)
	bismark2bedGraph --output $basename.bed \
	--dir $2 \
	--cutoff 5 \
	--CX \
	$1
    echo 'finished' $1
}

export -f BISMARK_BEDGRAPH

CHG_in=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/whole_genome/CHG*; do echo $f; done)
CHG_out=/global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHG_highcov

CHH_in=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/whole_genome/CHH*; do echo $f; done)
CHH_out=/global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHH_highcov

parallel BISMARK_BEDGRAPH ::: $CHG_in ::: $CHG_out
parallel BISMARK_BEDGRAPH ::: $CHH_in ::: $CHH_out

BISMARK_BEDGRAPH_CpG () {
	basename = $(basename $1)
	bismark2bedGraph --output $basename.bed \
	--dir $2 \
	--cutoff 5 \
	$1
    echo 'finished' $1
}

export -f BISMARK_BEDGRAPH_CpG
CpG_in=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/whole_genome/CpG*; do echo $f; done)
CpG_out=/global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CpG_highcov

parallel BISMARK_BEDGRAPH_CpG ::: $CpG_in ::: $CpG_out