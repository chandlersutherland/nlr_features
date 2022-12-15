#!/bin/bash
#SBATCH --job-name=coverage_py_runner
#SBATCH --partition=savio2
#SBATCH --qos=savio_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

#this script passes a bam file to the python script samtools_coverage, which outputs a file $BASENAME_clean_coverage.tsv which gives the mean depth over the NLRs
#can easily be made into a for loop to pass multiple files through 
#pass variable coverage_input, which should have a directory of sam files ready for coverage processing 

module load python
module load samtools/1.14
#very important to have this version 

cd $coverage_input
mkdir -p $coverage_input/sort_index

for file in *.sam
do 
	BASENAME=$(basename ${file} .sam)
	echo $BASENAME
	#first, convert to bam, sort and index for samtools coverage to run appropriately 
	samtools view -@ $SLURM_NTASKS -b $file |\
	samtools sort -@ 20 - -o $1/sort_index/${BASENAME}.bam
	samtools index $coverage_input/sort_index/${BASENAME}.bam
done 

cd $coverage_input/sort_index/

for file in *.bam 
do 
	#run the coverage file 
	python $HOME/e14/samtools_coverage.py ${file}
done 

#clean up working directory 
mkdir -p $coverage_input/coverage
mv *_clean_coverage.tsv $1/coverage
rm *coverage.tsv 