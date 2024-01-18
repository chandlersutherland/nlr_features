library(data.table)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(ggsignif)
library(ggpubr)
library(singscore)
library(GSEABase)
library(patchwork)

#repeat permutations.Rmd with 10,000 replicates to allow for running in the background 

zenodo_path <- "C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\Zenodo V2\\"
table <-  read.csv(paste(zenodo_path, 'all_gene_table.csv')) 

#set number of replicates throughout document
replicates <- 10000

#test 
replicates 

# write a general permutation test function that can be used across features 
# first, calculate the observed value
observe <- function(df, col_index, hv_status){
  observed <- mean(df%>%filter(HV==hv_status)%>% pull(col_index)%>%na.omit()) - 
    mean(df %>% filter(HV=='all_genes')%>%pull(col_index)%>%na.omit())
  observed
}

permute <- function(df, col_index, hv_status) {
  # build null distribution
  permutation = replicate(replicates, { 
    df <- df %>% drop_na(col_index)
    sample_small <- df[sample(nrow(df), nrow(df[df$HV == hv_status,]) , replace = FALSE), ] # sample the size of hvNLRs
    sample_large <- df[! rownames(df) %in% rownames(sample_small), ] # sample the rest
    mean(sample_small %>% pull(col_index))-mean(sample_large %>% pull(col_index)) # get expected mean
})
  permutation
}


p_calc <- function(df, col_index){
  o <- observe(df, col_index, 'hv')
  p <- permute(df, col_index, 'hv')
  pval <- mean(o > p)
  
  o2 <- observe(df, col_index, 'non-hv')
  p2 <- permute(df, col_index, 'non-hv')
  pval2 <- mean(o2 > p2)
  
  print(paste(colnames(df)[col_index]))
  print(paste('hv_p: ', pval, ' nonhv_p:', pval2))
}

#apply to TE distance
p_calc(table, 6)

observe_med <- function(df, col_index, hv_status){
  observed <- median(df%>%filter(HV==hv_status)%>%pull(col_index)%>%na.omit()) - 
    median(df %>% filter(HV=='all_genes')%>%pull(col_index)%>%na.omit())
  observed
}

permute_med <- function(df, col_index, hv_status) {
  # build null distribution
  permutation = replicate(replicates, { 
    df <- df %>% drop_na(col_index)
    sample_small <- df[sample(nrow(df), nrow(df[df$HV == hv_status,]) , replace = FALSE), ] # sample the size of hvNLRs
    sample_large <- df[! rownames(df) %in% rownames(sample_small), ] # sample the rest
    median(sample_small %>% pull(col_index))-median(sample_large %>% pull(col_index)) # get expected mean
})
  permutation
}

p_calc_med <- function(df, col_index){
  o <- observe_med(df, col_index, 'hv')
  p <- permute_med(df, col_index, 'hv')
  pval <- mean(o > p)
  
  o2 <- observe_med(df, col_index, 'non-hv')
  p2 <- permute_med(df, col_index, 'non-hv')
  pval2 <- mean(o2 > p2)
  
  print(paste(colnames(df)[col_index], ' median'))
  print(paste('hv_p: ', pval, ' nonhv_p:', pval2))
}

p_calc_med(table, 6)
p_calc(table, 5)

SRR17281085 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281085_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281086 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281086_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281087 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281087_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281088 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281088_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))

#clean up 
samples <- list(SRR17281085, SRR17281086, SRR17281087, SRR17281088)
names(samples) <- c('SRR17281085', 'SRR17281086', 'SRR17281087', 'SRR17281088')

clean <- function(df) {
  df2 <- df %>% merge(table, by='Gene', all=TRUE)
  df2$HV <- factor(df2$HV , levels=c("all_genes", "non-hv", "hv"))
  df3 <- df2 %>% subset(select=c('Gene', 'HV', 'mean_percent_methylation', 'count')) %>%
    filter(!is.na(mean_percent_methylation)) %>% 
    filter(Gene != 'AT4G16860') %>% 
    filter(Gene != 'AT1G58602')
  df3
}

clean_samples <- lapply(names(samples), function(x) clean(samples[[x]]))
names(clean_samples) <- c('SRR17281085', 'SRR17281086', 'SRR17281087', 'SRR17281088')

#create a sampled gene list weighted by counts 
density_sample <- function(df, hv_status){
  split_df <- df %>% filter(HV == hv_status)
  
  #calculate CG count density KDE function 
  d <- density(split_df$count)
  density_function <- approxfun(d)
  
  # write probabilities, and write NAs to zero as they are outliers in CG content (90+ per gene)
  probabilities <- density_function(df$count)
  probabilities[is.na(probabilities)] <- 0
  
  # sample genes, with weights according to their number of CGs, according to above function
  all_genes <- df$Gene
  sampled_genes <- data.frame(replicate(replicates, 
                                        sample(all_genes, 
                                               size = nrow(split_df), 
                                               replace = FALSE, 
                                               prob = probabilities)))
  
  return(sampled_genes)
  }


#write function that gets the difference in means between the permuted sample and the actual 
get_mean_difference <- function(gene_list, table) {
  permute_mean <- table %>% filter(Gene %in% gene_list) %>% pull(mean_percent_methylation) %>% mean()
  actual_mean <- table %>% filter(!Gene %in% gene_list)  %>% pull(mean_percent_methylation) %>% mean()
  return(actual_mean-permute_mean)
}

#calculate the pvalue 
pval <- function(sampled_df, table, hv_status, observed){
  differences <- list()
  for (x in 1:replicates) { 
    x <- get_mean_difference(sampled_df[[x]], table)
    differences <- append(differences, x)
    }
  
  p_val <- mean(differences < observed)
  
  return(p_val)
}

print('hv_nlr p value kernel density enrichment methylation')
hv_sampled_genes <- lapply(names(clean_samples), function(x) density_sample(clean_samples[[x]], 'hv'))
observed_hv <- lapply(names(clean_samples), function(x) observe(clean_samples[[x]], 3, 'hv'))

print(pval(hv_sampled_genes[[1]], clean_samples[[1]], 'hv', observed_hv[[1]]))
print(pval(hv_sampled_genes[[2]], clean_samples[[2]], 'hv', observed_hv[[2]]))
print(pval(hv_sampled_genes[[3]], clean_samples[[3]], 'hv', observed_hv[[3]]))
print(pval(hv_sampled_genes[[4]], clean_samples[[4]], 'hv', observed_hv[[4]]))

print('nonhv_nlr p value kernel density erinchment methylation') 
nhv_sampled_genes <- lapply(names(clean_samples), function(x) density_sample(clean_samples[[x]], 'non-hv'))
observed_nhv <- lapply(names(clean_samples), function(x) observe(clean_samples[[x]], 3, 'non-hv'))

print(pval(nhv_sampled_genes[[1]], clean_samples[[1]], 'non-hv', observed_nhv[[1]]))
print(pval(nhv_sampled_genes[[2]], clean_samples[[2]], 'non-hv', observed_nhv[[2]]))
print(pval(nhv_sampled_genes[[3]], clean_samples[[3]], 'non-hv', observed_nhv[[3]]))
print(pval(nhv_sampled_genes[[4]], clean_samples[[4]], 'non-hv', observed_nhv[[4]]))