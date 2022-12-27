#!/bin/bash
#SBATCH --job-name=coverage_py_runner
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --time=00:03:00
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

#first, convert to bam, sort and index for samtools coverage to run appropriately 
samtools sort -@ 20 $file -o $coverage_input/sort_index/${sample}.bam ${sample}_1_val_1_bismark_bt2_pe.bam
samtools index $coverage_input/sort_index/${sample}.bam

cd $coverage_input/sort_index/

#run the coverage file 
python $HOME/e14/samtools_coverage.py ${sample}.bam 

#clean up working directory 
mkdir -p $coverage_input/coverage
mv *_clean_coverage.tsv $coverage_input/coverage
rm *coverage.tsv 
