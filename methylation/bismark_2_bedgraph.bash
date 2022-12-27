#!/bin/bash
#SBATCH --job-name=bismark_extractor
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
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

#function takes two positional arguments: $1 being the .txt file to extract from, and $2 being the output directory
#cutoff set to 5 means that only cytosines with >5 reads are outputted, reducing file size and simplifying future processing 
BISMARK_BEDGRAPH () {
	basename=$(basename $1 _1_val_1_bismark_bt2_pe.deduplicated.txt)
	bismark2bedGraph --output $basename.bed \
	--dir $2 \
	--cutoff 5 \
	--CX \
	$1
    echo 'finished' $1
}

export -f BISMARK_BEDGRAPH
mkdir -p $extraction_output/bedGraph_highcov

CHG_in=$extraction_output/CHG_context_${sample}_1_val_1_bismark_bt2_pe.deduplicated.txt
CHG_out=$extraction_output/bedGraph_highcov

CHH_in=$extraction_output/CHH_context_${sample}_1_val_1_bismark_bt2_pe.deduplicated.txt
CHH_out=$extraction_output/bedGraph_highcov

BISMARK_BEDGRAPH $CHG_in $CHG_out 
BISMARK_BEDGRAPH $CHH_in $CHH_out

#need a separate function for CpG since it cannot take the --CX tag 
BISMARK_BEDGRAPH_CpG () {
	basename=$(basename $1 _1_val_1_bismark_bt2_pe.deduplicated.txt)
	bismark2bedGraph --output $basename.bed \
	--dir $2 \
	--cutoff 5 \
	$1
    echo 'finished' $1
}

export -f BISMARK_BEDGRAPH_CpG

CpG_in=$extraction_output/CpG_context_${sample}_1_val_1_bismark_bt2_pe.deduplicated.txt
CpG_out=$extraction_output/bedGraph_highcov 

BISMARK_BEDGRAPH_CpG $CpG_in $CpG_out 

cd $extraction_output/bedGraph_highcov
gunzip *.cov.gz 

echo 'finished generating context specific coverage files for ${sample}'
