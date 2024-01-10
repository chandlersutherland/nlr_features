R Notebook
================

Copyright (c) Chandler Sutherland Email:
<chandlersutherland@berkeley.edu>

Purpose: generate iTOL annotation files corresponding to the tree

``` r
library(tidyverse)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.1

    ## Warning: package 'purrr' was built under R version 4.3.1

    ## Warning: package 'dplyr' was built under R version 4.3.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.3     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(openxlsx)
library(ggsignif)
library(ggpubr)
```

First, add the itol branch length names and common gene names to the
methylation/expression/TE distance information

``` r
#read in the methylation, expression, and TE distance for all genes, filter to NLRs
zenodo_path <- "C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\Zenodo V2\\"
NLR_table <- read.csv(paste(zenodo_path, 'NLR_gene_table.csv'))

#add gene name information and the itol branch labels 
gene_names <- readxl::read_xlsx(path = "C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//gene_names_Daniil_tree.xlsx", col_names=c('itol_name', 'gene_name'))

names <- gene_names %>% separate(col=itol_name, sep=' ', into=c(NA, 'transcript_name'), remove=FALSE) %>% filter(itol_name != 'I264') %>% filter(itol_name != 'I0')
```

    ## Warning: Expected 2 pieces. Missing pieces filled with `NA` in 2 rows [135,
    ## 136].

``` r
names$Gene <- substr(names$transcript_name, 1, nchar(names$transcript_name)-2)

named_table <- merge(names, NLR_table, by='Gene', all.x=TRUE)
```

Starting with some binaries: HV/nonhv

``` r
itol_HV <- named_table[c('itol_name', 'HV')]
itol_HV <- itol_HV %>% mutate(HV_itol=recode(HV, `0`='-1'))
write.xlsx(itol_HV, 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\nlr_features\\figure_generation\\Source Data\\Figure 3\\3B\\itol_HV.xlsx')
```

Color gradient for expression and methylation

``` r
itol_TPM <- named_table[c('itol_name', 'log2_TPM')]
write.xlsx(itol_TPM, 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\nlr_features\\figure_generation\\Source Data\\Figure 3\\3B\\itol_tpm.xlsx')

itol_meth <- named_table[c('itol_name', 'meth_percentage')]
write.xlsx(itol_meth, 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\nlr_features\\figure_generation\\Source Data\\Figure 3\\3B\\itol_meth.xlsx')
```

TE could easily become a binary, which may be more informative anyway

``` r
#how about TE within 1kb
named_table['te_neighbor_500'] = named_table['te_dist'] <= 500
named_table$te_neighbor_500 <- named_table$te_neighbor_500 %>% as.integer() %>% as.character()

itol_te_500 <- named_table[c('itol_name', 'te_neighbor_500')]
itol_te_500 <- itol_te_500 %>% mutate(te_neighbor_500=recode(te_neighbor_500, '0'='-1'))
itol_te_500[is.na(itol_te_500)] <- '1'

write.xlsx(itol_te_500, 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\nlr_features\\figure_generation\\Source Data\\Figure 3\\3B\\itol_te_500.xlsx')
```

``` r
#HV color strip 
itol_HV <- itol_HV %>% mutate(HV_color=recode(HV, `1`='#F8766D', `0` = ''))
write.xlsx(itol_HV, 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\nlr_features\\figure_generation\\Source Data\\Figure 3\\3B\\itol_HV.xlsx')
```
