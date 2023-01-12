library(polyester)
library(Biostrings)

#The polyester library simulates RNAseq reads based on an input fold change matrix and a fasta file containing transcripts to simulate from
#Generate a fold change matrix, with only "1" values since I am just interested in coverage, not differential expression  

primary_fold_changes=matrix(c(1), nrow=27654, ncol=1)
primary_fasta_file = '/global/scratch/users/chandlersutherland/Athaliana/Athaliana_447_Araport11.cds_primaryTranscriptOnly.fa'

#convert to a biostrings format 
primary_fasta=readDNAStringSet(primary_fasta_file)

# set the reads_per_transcript: baseline mean number of reads for each transcript
# since I'm just determining mappability, not crucial, but derived this number from the average number of counts per gene and how that related to input reads. 
primary_readspertx=1092

#simulate four replicates 
simulate_experiment(primary_fasta_file, reads_per_transcript=primary_readspertx, 
	num_reps=c(4), fold_changes=primary_fold_changes, 
	readlen=50, paired=FALSE, 
	outdir='/global/scratch/users/chandlersutherland/e14/polyester/primary_simulated_reads_1006')
	