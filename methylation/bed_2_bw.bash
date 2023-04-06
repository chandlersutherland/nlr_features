#!/bin/bash
#SBATCH --job-name=deeptools 
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --time=00:20:00
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load parallel 
module load python/3.8.8

#I've installed the deeptools package using conda
source activate e14_deeptools 

#input is a bedgraph file
bed_to_big() {
	
	BASENAME=$(basename $1)
	CHROM_SIZE=$SCRATCH/e14/deeptools/Athaliana_447_TAIR10.genome.sizes
	DIR=$(dirname $(readlink -f ${1}))
	BW_DIR=${DIR}/bigwig; mkdir -p ${BW_DIR}
	
	#first, need to remove the track line from the bedgraph files, and save them as *.fixed.bed 
	grep -v track ${1} |
	sort -k1,1 -k2,2n > ${1}_sorted.bed
	
	#now, convert to bigwig 
	bedGraphToBigWig ${1}_sorted.bed $CHROM_SIZE $BW_DIR/${BASENAME}.bw 
	
	echo $1 'converted to bigwig' 
	
}


export -f bed_to_big

gunzip $bigwig_input/*.gz

input=$(for f in $bigwig_input/*_context_${sample}.bed; do echo $f; done)

parallel bed_to_big ::: $input 
