#!/bin/bash
#SBATCH --job-name=pal2nal
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

base=/global/scratch/users/chandlersutherland/e14/popgen/clades

#for each clade, identify the alignment file, then run pal2nal with the transcript fasta file. 
while read clade
do 
	#name inputs, initialize output directory 
	alignment=${base}/${clade}/*.afa
	fasta=${base}/${clade}/*.transcript.fa
	output=${base}/${clade}/${clade}.pal2nal.phy
	
	#actually run pal2nal 
	cd $HOME/programs/pal2nal.v14/
	./pal2nal.pl $alignment $fasta -output paml > $output
	
	echo "finished $clade"
done < /global/scratch/users/chandlersutherland/e14/popgen/clades.txt 
