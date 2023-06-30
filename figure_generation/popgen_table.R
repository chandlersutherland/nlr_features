library(tidyverse)
setwd("/global/scratch/users/chandlersutherland/e14/popgen/")  

#collect population genetics statistics into a single tsv for plotting 

#download the updated Gene table
clades <- read_delim("/global/scratch/users/chandlersutherland/e14/NLRCladeFinder/Atha_NLRome/Atha_NLRome_GeneTable.txt",delim = "\t")

#get just the Col0 genes and add the common names, write to col0_clean 
col0 <- clades %>% filter(Ecotype == 'ATHALIANA')
gene_names <- read_delim("/global/scratch/users/chandlersutherland/e14/hvnlr_clades/hvNLR/Analysis/AthaKnownGenes.txt",
                         delim='\t', col_names=c('Gene', 'name')) %>%
  mutate(name= gsub("\t","",as.character(name))) %>%
  separate(Gene, c(NA, 'Gene'), sep='_') %>%
  mutate(Gene=str_sub(Gene, end = -3))
col0_clean <- col0 %>% separate(Gene, c(NA, 'Gene', NA), sep='_') %>% merge(gene_names, on=Gene, all.x=TRUE)

#load omega values calculated by hyphy, busted_p results, and egglib pi and tajima's D 
omega <- read.csv('/global/scratch/users/chandlersutherland/e14/popgen/hyphy_w.csv')
busted_p <- read_csv('/global/scratch/users/chandlersutherland/e14/popgen/hyphy_busted_p.csv')
egglib <- read_tsv('/global/scratch/users/chandlersutherland/e14/popgen/egglib_pi.tsv') %>% rename(Clade=clade)

#combine popgen data into a single df 
col_popgen <- merge(col0_clean, omega, by='Clade', all.x=TRUE) %>% 
  merge(busted_p, all.x=TRUE)%>%
  merge(egglib, all.x=TRUE)%>% 
  mutate(busted_q = p.adjust(p, method='BH')) %>% #correct BUSTED p value for multiple testing 
  mutate(busted=case_when(busted_q < 0.05 ~ 'yes', busted_q >= 0.05 ~ 'no')) %>% #create categorical and remove old p value
  subset(select=-c(p, Clade_0, Clade_1, Clade_2, Clade_3, Ecotype, Allele, File, log_likelihood, GT, CT, CG, AT, AC, busted_q, `...1`)) %>%
  mutate(HV=recode(HV, `0` = "non-hv", `1`="hv")) %>%
  relocate("Gene", "Clade", "name", "HV")

write_tsv(col_popgen, '/global/scratch/users/chandlersutherland/e14/popgen/col0_popgen_stats.tsv')
