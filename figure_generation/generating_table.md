generating_table
================
Chandler Sutherland
2022-08-25

``` r
knitr::opts_knit$set(root.dir = 'C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14')
```

Copyright (c) Chandler Sutherland Email:
<chandlersutherland@berkeley.edu>

Purpose: combine outputs of feature investigation, population genetics,
and previously published work to create all_gene_table.csv and
NLR_gene_table.csv. These intermediate files will be published to
Zenodo.

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

``` r
#set path to output directory
zenodo_path <- "C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\Zenodo V2\\"
```

Read in the supplemental table from Lee and Chae 2020 and the HV status
from Prigozhin and Krasileva 2021

``` r
chae <- readxl::read_xlsx(path="C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\chae_2020_supplemental.xlsx")

clades <- read.table(file='C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\Atha_NLRome_GeneTable.txt', header=1) %>% 
  filter(Ecotype=='ATHALIANA') %>% 
  separate(Gene, c(NA, 'Gene'))%>%
  subset(select=c('Gene', 'Clade', 'HV')) %>% 
  distinct()
```

    ## Warning: Expected 2 pieces. Additional pieces discarded in 133 rows [1, 2, 3, 4, 5, 6,
    ## 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

``` r
#merge on Gene name 
table <- merge(chae, clades, by.x='gene', by.y='Gene', all.y=TRUE)

#Manually add information for At4G19050, which is missing from Chae analysis 
index <- which(table$gene == 'AT4G19050')
table[index, 3] = 'Chr4'
table[index, 4] = '10439983'
table[index, 5] = '10444121'
table[index, 6] = '-'
table[index, 7] = 'singleton'
table[index, 18] = 'C' 
```

QC

``` r
HV_count <- sum(table$HV == 1)
nHV_count <- sum(table$HV == 0)

paste("There are ", HV_count, "HV genes in my working table and ", nHV_count, " non-HV genes in my working table.")
```

    ## [1] "There are  35 HV genes in my working table and  97  non-HV genes in my working table."

``` r
hvNLR <- table %>% filter(HV=='hv') %>% pull('gene')
nonhvNLR <- table %>% filter(HV=='non-hv') %>% pull('gene') %>% unique()
```

Average the four rosette leaf biological replicates %CG from Williams
2022 for all genes.

``` r
#per gene average % CpG methylation 
SRR17281085 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281085_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281086 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281086_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281087 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281087_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))
SRR17281088 <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//SRR17281088_per_gene_meth_count.xlsx", skip=2, col_names = c('Index', 'Chrom', 'Gene', 'mean_percent_methylation', 'count'))

all_meth <- rbind(SRR17281085, SRR17281086, SRR17281087, SRR17281088) %>%
  group_by(Gene) %>% summarize(meth_percentage=weighted.mean(mean_percent_methylation, count/sum(count)))

#compare weighted mean to mean 
```

Average the TPM from the four biological replicates from Williams 2022
for all genes. Merge with methylation.

``` r
#TPM 
TPM <- read.table("C:/Users/chand/Box Sync/Krasileva_Lab/Research/chandler/Krasileva Lab/E14/e14_R/tpm_matrix.csv", header=TRUE, sep=",", row.names=2)%>% dplyr::select(-X)

TPM_df <- TPM %>% mutate('mean' = rowMeans(TPM[1:4])) 
all_expression <- merge(TPM_df, all_meth, by.x=0, by.y='Gene', all=TRUE)
all_expression <- all_expression %>% mutate(log2_TPM = log2(all_expression$mean+1))

#have to manually set the unmappable genes, since they just look like 0 
unmappable <- c('AT1G58807', 'AT1G58848', 'AT1G59124', 'AT1G59218')

all_expression$log2_TPM[all_expression$Row.names %in% unmappable] <- NA   
all_expression$mean[all_expression$Row.names %in% unmappable] <- NA 
```

Read in TE distance information, convert to kbp distance. Merge with
expression and methylation.

``` r
#TE distance
TE_table <- read.table(file='C:/Users/chand/Box Sync/Krasileva_Lab/Research/chandler/Krasileva Lab/E14/e14_R/Atha_genes_with_TE_dist.tsv', header=1) %>% distinct()
#converte TE distance to kb 
TE <- merge(TE_table, all_expression, by.x = 'Gene', by.y='Row.names', all=TRUE) %>%
  mutate(te_dist = te_dist/1000)
```

Retain only relevant columns and clean categorical variables.

``` r
#clean and concentrate 
merged <- merge(TE, table, by.x='Gene', by.y='gene', all=TRUE)
all_table <- merged %>% subset(select=c('Gene', 'HV', 'log2_TPM', 'meth_percentage', 'te_dist')) %>%
  mutate(HV=recode(HV, `0` = "non-hv", `1`="hv"))

all_table$HV[is.na(all_table$HV)]<-'all_genes'
all_table$HV <- factor(all_table$HV , levels=c("all_genes", "non-hv", "hv"))

HV_count <- sum(all_table$HV == 'hv')
nHV_count <- sum(all_table$HV == 'non-hv')

paste("There are ", HV_count, "HV genes in my working table and ", nHV_count, " non-HV genes in my working table.")
```

    ## [1] "There are  35 HV genes in my working table and  97  non-HV genes in my working table."

Pi and D were calculated across the genome based on the 1001 genome VCF,
and across NLRs based on the protein alignment from long read
sequencing. Combine these results, and add to the all gene table.

``` r
#popgen stats from egglib generated on NLRome coding sequences as described in egglib_nlr.ipynb 
NLR_popgen <- read_tsv("//wsl.localhost//Ubuntu//home//chandlersutherland//scratch//egglib_results.tsv") %>% 
  filter(domain=='cds')%>%
  mutate(Pi=Pi_raw/lseff) %>%
  subset(select=c('Gene', 'D','Pi'))
```

    ## New names:
    ## Rows: 518 Columns: 13
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (4): clade, domain, Gene, HV dbl (9): ...1, D, PiN, PiS, Pi_raw, dna_aln,
    ## lseff, start, stop
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
#rest of the genome from vcf as described in egglib_whole_genome.ipynb
whole_genome <- read_csv("//wsl.localhost//Ubuntu//home//chandlersutherland//scratch//cds_egglib_all.csv") %>% 
  subset(select=c('gene', 'Pi_by_lseff', 'D')) %>%
  dplyr::rename('Gene'='gene', 'Pi'='Pi_by_lseff') %>% distinct()
```

    ## New names:
    ## Rows: 25703 Columns: 8
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): gene dbl (7): ...1, D, Pi_raw, cds_length, lseff, Pi_by_dna, Pi_by_lseff
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
#replace NLRs in whole genome with more accurate long read versions 
NLRs <- NLR_popgen %>% pull('Gene') %>% unique()
no_NLRs <- whole_genome %>% filter(!Gene %in% NLRs)

popgen_table <- rbind(no_NLRs, NLR_popgen) %>% arrange(Gene)

final_table <- merge(all_table, popgen_table) %>% distinct()
HV_count <- sum(final_table$HV == 'hv')
nHV_count <- sum(final_table$HV == 'non-hv')

final_table %>% filter(HV!='all_genes') %>% distinct() %>% group_by(Gene) %>% filter(n()>1)
```

    ## # A tibble: 0 × 7
    ## # Groups:   Gene [0]
    ## # ℹ 7 variables: Gene <chr>, HV <fct>, log2_TPM <dbl>, meth_percentage <dbl>,
    ## #   te_dist <dbl>, Pi <dbl>, D <dbl>

``` r
paste("There are ", HV_count, "HV genes in my working table and ", nHV_count, " non-HV genes in my working table.")
```

    ## [1] "There are  35 HV genes in my working table and  97  non-HV genes in my working table."

Write this table to the Zenodo directory.

``` r
colnames(final_table)
```

    ## [1] "Gene"            "HV"              "log2_TPM"        "meth_percentage"
    ## [5] "te_dist"         "Pi"              "D"

``` r
write.csv(final_table, paste(zenodo_path, 'all_gene_table.csv'))
```

Create an NLR table with some more gene-wide information, things like
cluster status, distance to nearest NLR, some NLR specific popgen stats

``` r
#read in col0_popgen_stats, generated by the egglib_nlr.ipynb script 
NLR_popgen <- read_tsv("//wsl.localhost//Ubuntu//home//chandlersutherland//scratch//egglib_results.tsv") %>%
  mutate(Pi=Pi_raw/lseff)%>% #Pi normalized by effective sites 
  mutate(PiNPiS = PiN/PiS) %>%
  filter(domain=='cds') %>% 
  subset(select=c("Gene", 'HV', 'D', 'PiN', 'PiS',  'Pi', 'PiNPiS'))
```

    ## New names:
    ## Rows: 518 Columns: 13
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (4): clade, domain, Gene, HV dbl (9): ...1, D, PiN, PiS, Pi_raw, dna_aln,
    ## lseff, start, stop
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
#add in gene names 
names <- read_csv(paste(zenodo_path, "Atha_NLR_common_names.csv", sep='')) 
```

    ## New names:
    ## Rows: 132 Columns: 3
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (2): Gene, name dbl (1): ...1
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
NLR_popgen <- merge(NLR_popgen, names, all.x=TRUE)

#add in cluster information 
NLR_table <- merged %>% 
  filter(!is.na(HV)) %>% #nlr only 
  subset(select=c('Gene', 'chr', 'start.x', 'end.x', 'strand.x', 
                  'te_dist', 'meth_percentage', 'log2_TPM', 
                  'cluster_type', 'distance_front', 'distance_back', 'cluster', 'whyCluster',
                  'Nterm', 
                  'Clade')) %>% #relevant columns
  merge(NLR_popgen, by='Gene', all.x=T) %>%
  dplyr::rename('start'='start.x', 'end'='end.x', 'strand'='strand.x')

#change cluster definition slightly, we just want distance clusters not similarity based clusters 
NLR_table <- NLR_table %>% 
  mutate(cluster=case_when((distance_front > 50000 & distance_back > 50000) ~ 'singleton', 
                           .default=cluster)) %>% 
  mutate(clustered=case_when(cluster!='singleton' ~ 'clustered', 
                             cluster=='singleton' ~ 'singleton', 
                             whyCluster=='similarity' ~ 'singleton')) %>%
  unique()
colnames(NLR_table)
```

    ##  [1] "Gene"            "chr"             "start"           "end"            
    ##  [5] "strand"          "te_dist"         "meth_percentage" "log2_TPM"       
    ##  [9] "cluster_type"    "distance_front"  "distance_back"   "cluster"        
    ## [13] "whyCluster"      "Nterm"           "Clade"           "HV"             
    ## [17] "D"               "PiN"             "PiS"             "Pi"             
    ## [21] "PiNPiS"          "...1"            "name"            "clustered"

``` r
nrow(NLR_table)
```

    ## [1] 132

``` r
#add Mutation.Probability.Score from Monroe 2020
mutation <- read.csv("C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E8//Ath_somatic_mutations_gene_level_data.csv") %>% dplyr::rename('Gene'='gene') %>% subset(select=c('Gene', 'Mutation.Probability.Score'))
NLR_table <- merge(NLR_table, mutation, by='Gene', all.x=T)

HV_count <- sum(NLR_table$HV == 'hv')
nHV_count <- sum(NLR_table$HV == 'non-hv')

paste("There are ", HV_count, "HV genes in my NLR table and ", nHV_count, " non-HV genes in my NLR table.")
```

    ## [1] "There are  35 HV genes in my NLR table and  97  non-HV genes in my NLR table."

``` r
write.csv(NLR_table, file=paste(zenodo_path, 'NLR_gene_table.csv'))
```
