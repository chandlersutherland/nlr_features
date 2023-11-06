#!/bin/bash
#SBATCH --job-name=methylation_processing
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:01:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

cd /global/home/users/chandlersutherland/nlr_features/methylation

#Define the input directory, with fastq files ready for trimming
trim_input="/global/scratch/users/chandlersutherland/e14/twim_williams/williams_tissues"
trim_output="/global/scratch/users/chandlersutherland/e14/trim_williams/williams_tissues"

mkdir -p $trim_output

while read sample; do sbatch --job-name=$sample.trim_galore --export=trim_input=$trim_input,trim_output=$trim_output,sample=$sample \
 -A co_minium /global/home/users/chandlersutherland/nlr_features/methylation/trim_galore.sh; done < samples.txt

#run bismark alignment 
genome="/global/scratch/users/chandlersutherland/phytozome/Athaliana/Araport11/assembly"
bismark_output="/global/scratch/users/chandlersutherland/e14/bismark/bam_files/tissues"

mkdir -p $bismark_output

while read sample; do sbatch --job-name=$sample.bismark --export=trim_output=$trim_output,genome=$genome,bismark_output=$bismark_output,sample=$sample -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/bismark.bash; done < samples.txt 

#check read coverage of NLRs, generating a clean_coverage.tsv file for each .sam file in the coverage_input directory  
coverage_input=$bismark_output
while read sample; do sbatch --job-name=$sample.nlr_coverage --export=coverage_input=$coverage_input,sample=$sample -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/samtools_coverage_runner.bash; done < samples.txt

#deduplicate files + extraction 
extraction_output="/global/scratch/users/chandlersutherland/e14/bismark/extraction/tissues"

mkdir -p $extraction_output

while read sample; do sbatch --job-name=$sample.dedup --export=deduplicate_input=$bismark_output,extraction_output=$extraction_output,sample=$sample -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/bismark_deduplicate.bash; done < samples.txt

#finish extraction to coverage files, and filter to >5 reads per cytosine 
while read sample; do sbatch --job-name=$sample.extract --export=extraction_output=$extraction_output,sample=$sample\
 -A co_minium /global/home/users/chandlersutherland/nlr_features/methylation/bismark_2_bedgraph.bash; done < samples.txt

#generate bigwig files for viewing in IGV of the high coverage extracted cytosines 
bigwig_input=$extraction_output/bedGraph_highcov

while read sample; do sbatch --job-name=$sample.bw --export=bigwig_input=$bigwig_input,sample=$sample -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/bed_2_bw.bash; done < samples.txt 

#calculate the per gene methylation percentage for each context 
cov_dir='/global/scratch/users/chandlersutherland/e14/bismark/extraction/tissues/bedGraph_highcov'
gene_positions='/global/scratch/users/chandlersutherland/Athaliana/Atha_genes.bed'


while read sample; do sbatch --job-name=$sample.per_gene_meth --export=cov_dir=$cov_dir,sample=$sample,gene_positions=$gene_positions -A co_minium \
/global/home/users/chandlersutherland/nlr_features/methylation/per_gene_methylation_runner.bash; done < samples.txt 
