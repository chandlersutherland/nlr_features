#!/bin/bash
#SBATCH --job-name=hyphy_fubar
#SBATCH --partition=savio4_htc
#SBATCH --qos=minium_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --mail-user=chandlersutherland@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --error=/global/home/users/chandlersutherland/slurm_stderr/slurm-%j.out
#SBATCH --output=/global/home/users/chandlersutherland/slurm_stdout/slurm-%j.out

module load python
module load parallel 
source activate e14 

#goal: run FEL in parallel on all NLR clades 
base=/global/scratch/users/chandlersutherland/e14/popgen/clades
clade=$(cat /global/scratch/users/chandlersutherland/e14/popgen/clades.txt)

FEL_RUN(){
	#function that takes in clade name, finds alignment and tree file, creates log file, then runs FEL 
	alignment=${base}/${1}/popgenome/${1}.pal2nal.fas
	tree=${base}/${1}/RAxML*.out
	
	log_file=${base}/${1}/hyphy_fel.log
	echo "Running Hyphy on $clade"
	
	#actually run hyphy in multithreading mode to hopefully speed up 
	hyphy fel CPU=4 --alignment ${alignment} --tree ${tree} | tee -a $log_file
	echo "finished $1"
}

export base=$base
export -f FEL_RUN

parallel FEL_RUN ::: $clade 