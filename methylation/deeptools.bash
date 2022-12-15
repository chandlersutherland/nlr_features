#!/bin/bash
#SBATCH --job-name=deeptools 
#SBATCH --partition=savio2
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --time=00:20:00
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load parallel 


#input is a file


bed_to_big() {
	module load python/3.8.8
	source activate e14_deeptools
	
	BASENAME=$(basename $1)
	CHROM_SIZE=$SCRATCH/e14/deeptools/Athaliana_447_TAIR10.genome.sizes
	DIR=$(dirname $(readlink -f ${1}))
	BW_DIR=${DIR}/bigwig; mkdir -p ${BW_DIR}
	
	#first, need to remove the track line from the bedgraph files, and save them as *.fixed.bed 
	grep -v track ${1} |
	sort -k1,1 -k2,2n > ${1}_sorted.bed
	
	#now, convert to bigwig 
	bedGraphToBigWig ${1}_sorted.bed $CHROM_SIZE $BW_DIR/${BASENAME}.bw 
	
	echo 'converted ${1} to bigwig' 
	
	conda deactivate
	module purge 
}




#input is a bw file 
plotter(){
	#ok, now to do some profile plots 
	module load python 
	source activate e14 

	density_file= ${1}
	THREADS=$SLURM_NTASKS
	HV_BED=/global/home/users/chandlersutherland/e14/data/hv_NLR.bed
	NONHV_BED=/global/home/users/chandlersutherland/e14/data/nonhv_NLR2.bed
	DIR=$(dirname $(readlink -f ${1}))
	BASENAME=$(basename ${1})
	DEEPTOOLS_DIR=${DIR}/deeptools_out
	mkdir -p $DEEPTOOLS_DIR
	
	
	#not currently working on interactive mode, but no error messages?
	computeMatrix scale-regions -S $density_file \
		-R ${HV_BED} ${NONHV_BED} \
		--regionBodyLength 4000 \
		-o ${DEEPTOOLS_DIR}/${BASENAME}.both.mat.gz
	
	plotProfile -m ${DEEPTOOLS_DIR}/${BASENAME}.both.mat.gz \
		-out ${DEEPTOOLS_DIR}/${BASENAME}.both.pdf \
		--numPlotsPerRow 2 \
		--plotTitle "${1}" \
		--outFileNameData ${DEEPTOOLS_DIR}/${BASENAME}.both.tab 
		
	echo 'finished' ${BASENAME}
}

export -f plotter
export -f bed_to_big
echo 'parallel finished'

CHG_highcov=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHG_highcov; do echo $f; done)
CHH_out=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHH_highcov; do echo $f; done)
CpG_out=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CpG_highcov; do echo $f; done)

CHG=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHG; do echo $f; done)
CHH=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CHH; do echo $f; done)
CpG=$(for f in /global/scratch/users/chandlersutherland/e14/bismark/extraction/bedGraph/whole_genome/CpG
; do echo $f; done)
