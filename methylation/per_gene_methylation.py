import os 
import sys
import pandas as pd 
import glob 
import numpy as np 

#this python file does the per gene processing after bismark methylation extractor, returning a per gene average methylation
#In the case of CG context, the symmetrical CGs are averaged as they are not statistically independent  
#The output is a tsv file with the average methylation % per gene, and the number of high coverage cytosines (or CG sites for CG context) per gene 

#input is an unzipped .cov file 
cov_file=str(sys.argv[1])
methylation = pd.read_csv(cov_file, skiprows=1, sep = '\t', names=['Chrom', 'start_pos', 'end_pos', 'meth_percentage', 'count_methylated', 'count_unmethylated'])
filename=cov_file.split('/')[-1]
context=filename.split('_')[0]
sample=filename.split('_')[2].split('.')[0]
print(filename+' loaded successfully')

#initialize output directory (same as input directory)
f=len(cov_file.split('/'))
b=cov_file.split('/')[0:f-1]
out_dir='/'.join(b)
print('Output directory is '+out_dir)

#input necessary gene name files variables (bed files, chromosomes)
#second passed argument should be the bed file for the entire genome
chroms = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
pos_file=str(sys.argv[2])
all_positions = pd.read_csv(pos_file, sep = '\t',header=0, names=['Chrom', 'start', 'end', 'strand', 'gene'], index_col=False)

#step 1: filter out ChrC and ChrM, since we don't care about them 
met=methylation.loc[(methylation['Chrom']!='ChrC') & (methylation['Chrom']!='ChrM')] 
print('Filtering out ChrC and ChrM...')

#step 2: Average the symmetrical CG pairs since they are not independent 
def cg_average(df):
    cpg_cov=df.sort_values(['Chrom','start_pos'])
    
    avg_cpg = []
    
    for i in range(0,len(cpg_cov)-1):
    #find and average the symmetrical CpG pairs 
    #this relies on sorting by chromosome and then by position, searching for neighbors
        if cpg_cov.iloc[i, 1] == cpg_cov.iloc[i+1,1]-1:
            obj = {'Chrom':cpg_cov.iloc[i,0],
                   'start_pos':cpg_cov.iloc[i,1],
                   'end_pos':cpg_cov.iloc[i+1,1],
                   'meth_percentage':np.mean([cpg_cov.iloc[i,3],
                                              cpg_cov.iloc[i+1,3]])
            }
            avg_cpg.append(obj)
        else:
            continue
    avg_cpg_df = pd.DataFrame.from_dict(avg_cpg)
    return avg_cpg_df

if context=='CpG':
    print('CpG context detected, averaging symmetrical positions')
    met=cg_average(met)
    print('averaging complete')
else: 
    print('CHH or CHG context, no averaging necessary')
        

#step 3: assign cytosines to genes. Define a few functions first 

#takes in a dataframe and a list of chromosomes, and returns a list of dataframes broken by chromosome.
#Run on positions to create pos 
def break_chroms(df, chroms):
  df_chroms = []
  for i in chroms: 
    chrom_met = df[df['Chrom'] == i]
    df_chroms.append(chrom_met)
  return df_chroms

#takes in each individual chunk by chromosome, and "names" it, assigning it to the correct gene  
def gene_namer(pos, met):
  for j in range(0, len(pos)):
    start = pos.iloc[j, 1]
    end = pos.iloc[j, 2]
    gene = pos.iloc[j, 4]
    met.loc[
      (met['Pos'].between(start,end)), 'Gene'] = gene
  return met.dropna()
  
#applies namer to each broken chromosome, then re-concatenates to a final value 
def bring_together(pos_chroms, df):
  final = pd.concat(
    [gene_namer(pos_chroms[i], 
           break_chroms(df, chroms)[i]) 
    for i in range(0,5)]
    )
  return final.sort_values('Gene').reset_index()

#Applies bring_together, and aggregates by gene to return average methylation per gene and the count of high cov cytosines per gene   
def merge(df, pos_chroms):
  cov_cpg2 = df.rename(columns={"start_pos":"Pos"})
  cov_cpg3 = bring_together(pos_chroms, cov_cpg2).reset_index()

  #add in a line that does the high cov CG count per gene 
  cpg_per_gene = cov_cpg3.groupby(['Chrom', 'Gene']).agg({
      'meth_percentage': ['mean', 'count']
  }).reset_index()

  return cpg_per_gene

pos_chroms = break_chroms(all_positions, chroms)  
print('Chromosomes broken')
all_gene_count=merge(met,pos_chroms) 
print('Gene names assigned')

#step 4: save file to a tsv 
#untested
all_gene_count.to_csv(out_dir+'/'+context+'_per_gene_met_'+sample+'.tsv', sep='\t', header=True, index=False)
print('Gene counts saved to '+out_dir+'/'+context+'_per_gene_met_'+sample+'.tsv')

print('Finished with '+sample+' in the '+context+' context!!')
