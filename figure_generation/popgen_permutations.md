popgen_permutations
================
Chandler Sutherland
2024-01-11

Author: Chandler Sutherland Copyright (c) Chandler Sutherland Email:
<chandlersutherland@berkeley.edu>

Goal: Use permutation tests to determine number of NLRs expected by
chance in tails of empirical distribution.

``` r
library(data.table)
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.1

``` r
library(tidyverse)
```

    ## Warning: package 'purrr' was built under R version 4.3.1

    ## Warning: package 'dplyr' was built under R version 4.3.1

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.3     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::between()     masks data.table::between()
    ## ✖ dplyr::filter()      masks stats::filter()
    ## ✖ dplyr::first()       masks data.table::first()
    ## ✖ lubridate::hour()    masks data.table::hour()
    ## ✖ lubridate::isoweek() masks data.table::isoweek()
    ## ✖ dplyr::lag()         masks stats::lag()
    ## ✖ dplyr::last()        masks data.table::last()
    ## ✖ lubridate::mday()    masks data.table::mday()
    ## ✖ lubridate::minute()  masks data.table::minute()
    ## ✖ lubridate::month()   masks data.table::month()
    ## ✖ lubridate::quarter() masks data.table::quarter()
    ## ✖ lubridate::second()  masks data.table::second()
    ## ✖ purrr::transpose()   masks data.table::transpose()
    ## ✖ lubridate::wday()    masks data.table::wday()
    ## ✖ lubridate::week()    masks data.table::week()
    ## ✖ lubridate::yday()    masks data.table::yday()
    ## ✖ lubridate::year()    masks data.table::year()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggsignif)
library(ggpubr)
library(patchwork)
```

    ## Warning: package 'patchwork' was built under R version 4.3.1

load data

``` r
#Change to your path to the zenodo download to repeat 
zenodo_path <- "C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\E14\\Zenodo V2\\"
#load in all gene data 
table <- read.csv(paste(zenodo_path, 'all_gene_table.csv'))
#set number of replicates 
replicates <- 1000
```

# D and Pi difference in means

``` r
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
  pval <- mean(o < p)
  
  o2 <- observe(df, col_index, 'non-hv')
  p2 <- permute(df, col_index, 'non-hv')
  pval2 <- mean(o2 < p2)
  
  print(paste(colnames(df)[col_index]))
  print(paste('hv_p: ', pval, ' nonhv_p:', pval2))
}
```

Pi permutation test difference in means:

``` r
#apply to Pi
p_calc(table, 7)
```

    ## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
    ## ℹ Please use `all_of()` or `any_of()` instead.
    ##   # Was:
    ##   data %>% select(col_index)
    ## 
    ##   # Now:
    ##   data %>% select(all_of(col_index))
    ## 
    ## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## [1] "Pi"
    ## [1] "hv_p:  0  nonhv_p: 0"

``` r
#apply to D
p_calc(table, 8)
```

    ## [1] "D"
    ## [1] "hv_p:  0.006  nonhv_p: 1"

# D and Pi top 5%

We are interested in the probability of observing NLRs, hvNLRs, or
non-hvNLRs in the top 5% of an empirical distribution

``` r
#function that takes in table, column name as a string, and n replicates, and outputs probability of observing NLRs, hvNLRs, and non-hvNLRs in the top 5% of the distribution by chance 
on_top <- function(table, col_name, n){
  
  q <- quantile(table %>% pull(col_name), c(.05, .95), na.rm=TRUE)
  print(paste('calculating probability of ', col_name, '. The 95th percentile is: ', q[2]))
  
  observed_nlr <- table[table['HV'] != 'all_genes' & table[col_name]>q[2],] %>% nrow()
  
  permutation_nlr = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV != 'all_genes',]) , replace = FALSE), ] # sample the size of all NLRs
  sample_small[sample_small[col_name]>q[2],] %>% nrow()  # get expected above 95% 
  })
  
  print(paste('observed NLRs in top 5%:', observed_nlr))
  print(paste('expected NLRs in top 5% by chance:', mean(permutation_nlr)))
  print(paste('p value non-hv in top 5%:', mean(observed_nlr < permutation_nlr)))
  
  observed_hv <- table[table['HV'] == 'hv' & table[col_name]>q[2],] %>% nrow()
  permutation_hv = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV == 'hv',]) , replace = FALSE), ] 
  sample_small[sample_small[col_name]>q[2],] %>% nrow()  # get expected above 95% 
  })
  
  print(paste('observed hvNLRs in top 5%:', observed_hv))
  print(paste('expected hvNLRs in top 5% by chance:', mean(permutation_hv)))
  print(paste('p value hv in top 5%:', mean(observed_hv < permutation_hv)))
  
  observed_nhv <- table[table['HV'] == 'non-hv' & table[col_name]>q[2],] %>% nrow()
  
  permutation_nhv = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV == 'non-hv',]) , replace = FALSE), ] # sample the size of non-hvNLRs
  sample_small[sample_small[col_name]>q[2],] %>% nrow()  # get expected above 95% 
  })
  
  print(paste('observed hvNLRs in top 5%:', observed_nhv))
  print(paste('expected hvNLRs in top 5% by chance:', mean(permutation_nhv)))
  print(paste('p value hv in top 5%:', mean(observed_hv < permutation_nhv)))
  
}

on_top(table, 'Pi', 10000)
```

    ## [1] "calculating probability of  Pi . The 95th percentile is:  0.0140489927540266"
    ## [1] "observed NLRs in top 5%: 63"
    ## [1] "expected NLRs in top 5% by chance: 9.4217"
    ## [1] "p value non-hv in top 5%: 0"
    ## [1] "observed hvNLRs in top 5%: 35"
    ## [1] "expected hvNLRs in top 5% by chance: 2.4933"
    ## [1] "p value hv in top 5%: 0"
    ## [1] "observed hvNLRs in top 5%: 28"
    ## [1] "expected hvNLRs in top 5% by chance: 7.0036"
    ## [1] "p value hv in top 5%: 0"

``` r
on_top(table, 'D', 10000)
```

    ## [1] "calculating probability of  D . The 95th percentile is:  1.50881426721351"
    ## [1] "observed NLRs in top 5%: 1"
    ## [1] "expected NLRs in top 5% by chance: 12.102"
    ## [1] "p value non-hv in top 5%: 1"
    ## [1] "observed hvNLRs in top 5%: 0"
    ## [1] "expected hvNLRs in top 5% by chance: 3.2033"
    ## [1] "p value hv in top 5%: 0.9629"
    ## [1] "observed hvNLRs in top 5%: 1"
    ## [1] "expected hvNLRs in top 5% by chance: 8.9656"
    ## [1] "p value hv in top 5%: 0.9998"

``` r
on_bottom <- function(table, x, n){
  
  q <- quantile(table %>% pull(x), c(.05, .95), na.rm=TRUE)
  print(paste('calculating probability of ', x, '. The 5th percentile is: ', q[1]))
  
  observed_nlr <- table[table['HV'] != 'all_genes' & table[x]<q[1],] %>% nrow()
  
  permutation_nlr = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV != 'all_genes',]) , replace = FALSE), ] # sample the size of all NLRs
  sample_small[sample_small[x]<q[1],] %>% nrow()  # get below 5% 
  })
  
  print(paste('observed NLRs in bottom 5%:', observed_nlr))
  print(paste('expected NLRs in bottom 5% by chance:', mean(permutation_nlr)))
  print(paste('p value non-hv in bottom 5%:', mean(observed_nlr < permutation_nlr)))
  
  observed_hv <- table[table['HV'] == 'hv' & table[x]<q[1],] %>% nrow()
  
  permutation_hv = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV == 'hv',]) , replace = FALSE), ] # sample the size of hvNLRs
  sample_small[sample_small[x]<q[1],] %>% nrow()  # get expected above 95% 
  })
  
  print(paste('observed hvNLRs in bottom 5%:', observed_hv))
  print(paste('expected hvNLRs in bottom 5% by chance:', mean(permutation_hv)))
  print(paste('p value hv in bottom 5%:', mean(observed_hv < permutation_hv)))
  
  observed_nhv <- table[table['HV'] == 'non-hv' & table[x]<q[1],] %>% nrow()
  
  permutation_nhv = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV == 'non-hv',]) , replace = FALSE), ] # sample the size of non-hvNLRs
  sample_small[sample_small[x]<q[1],] %>% nrow()  # get expected above 95% 
  })
  
  print(paste('observed non-hvNLRs in bottom 5%:', observed_nhv))
  print(paste('expected non-hvNLRs in bottom 5% by chance:', mean(permutation_nhv)))
  print(paste('p value non-hv in bottom 5%:', mean(observed_hv > permutation_nhv)))
  
}

on_bottom(table, 'Pi', 1000)
```

    ## [1] "calculating probability of  Pi . The 5th percentile is:  0.000316624887148824"
    ## [1] "observed NLRs in bottom 5%: 0"
    ## [1] "expected NLRs in bottom 5% by chance: 9.467"
    ## [1] "p value non-hv in bottom 5%: 1"
    ## [1] "observed hvNLRs in bottom 5%: 0"
    ## [1] "expected hvNLRs in bottom 5% by chance: 2.523"
    ## [1] "p value hv in bottom 5%: 0.916"
    ## [1] "observed non-hvNLRs in bottom 5%: 0"
    ## [1] "expected non-hvNLRs in bottom 5% by chance: 6.754"
    ## [1] "p value non-hv in bottom 5%: 0"

``` r
on_bottom(table, 'D', 1000)
```

    ## [1] "calculating probability of  D . The 5th percentile is:  -1.78922602865888"
    ## [1] "observed NLRs in bottom 5%: 46"
    ## [1] "expected NLRs in bottom 5% by chance: 12.05"
    ## [1] "p value non-hv in bottom 5%: 0"
    ## [1] "observed hvNLRs in bottom 5%: 0"
    ## [1] "expected hvNLRs in bottom 5% by chance: 3.251"
    ## [1] "p value hv in bottom 5%: 0.958"
    ## [1] "observed non-hvNLRs in bottom 5%: 46"
    ## [1] "expected non-hvNLRs in bottom 5% by chance: 8.921"
    ## [1] "p value non-hv in bottom 5%: 0"

Check hvNLRs in bottom and top of empirical distribution of D

``` r
x <- 'D'
n <- 1000
q <- quantile(table %>% pull(x), c(.05, .95), na.rm=TRUE)
observed_hv <- table %>% filter(HV == 'hv') %>% filter(D<q[1]) %>% filter(D>q[2]) %>% nrow()

permutation_hv = replicate(n, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV == 'hv',]), replace = FALSE), ] 
  sample_small %>% filter(D<q[1] | D>q[2]) %>% nrow()  # get expected above 95% 
  })

print(paste('observed hvNLRs in either tail:', observed_hv))
```

    ## [1] "observed hvNLRs in either tail: 0"

``` r
print(paste('expected hvNLRs in either tail by chance:', mean(permutation_hv)))
```

    ## [1] "expected hvNLRs in either tail by chance: 3.389"

``` r
print(paste('p value hv in either tail:', mean(observed_hv < permutation_hv)))
```

    ## [1] "p value hv in either tail: 0.971"

Supplemental figure 3, which shows the distribution of D and Pi with the
percentiles

``` r
d_distro <- ggplot(data=table)+geom_histogram(aes(x=D))+
   geom_segment(x=1.50913785453558, xend=1.50913785453558, yend=0, y=20000, aes(color='5th and 95th percentiles'))+
  geom_segment(x=-1.7870429615165, xend=-1.7870429615165, yend=0, y=20000, aes(color='5th and 95th percentiles'))+
  geom_segment(data= table %>% filter(HV=='hv'), aes(
      x= D,
        y = 0,
        xend = D,
        yend = -250, color='hvNLR'), alpha = 0.75, linewidth=0.5) +
    labs(colour = c('hvNLR', 'nonhvNLR')) +
geom_segment(data= table %>% filter(HV=='non-hv'), aes(
      x= D,
        y = -250,
        xend = D,
        yend = -500, color='non-hvNLR'), alpha = 0.75, linewidth=0.5) +
    labs(colour = c('hvNLR', 'nonhvNLR')) +
  theme_classic() +
  theme(legend.position='right', text=element_text(size=10), legend.title=element_blank())+
 # labs(title='Empirical Distribution of D values')+
  xlab("Tajima's D")+
  ylab('count')+
  scale_color_manual(values=c('red', '#F8766D', '#00BFC4'))

pi_distro <- ggplot(data=table)+geom_histogram(aes(x=Pi), bins=50)+
  geom_segment(x=0.000316624887148824, xend=0.000316624887148824, yend=0, y=20000, aes(color='5th and 95th percentiles'))+
  geom_segment(x=0.0140259207477425, xend=0.0140259207477425, yend=0, y=20000, aes(color='5th and 95th percentiles'))+
  geom_segment(data= table %>% filter(HV=='hv'), aes(
      x= Pi,
        y = 0,
        xend = Pi,
        yend = -500, color='hvNLR'), alpha = 0.75, linewidth=0.5) +
    labs(colour = c('hvNLR', 'nonhvNLR')) +
geom_segment(data= table %>% filter(HV=='non-hv'), aes(
      x= Pi,
        y = -500,
        xend = Pi,
        yend = -1000, color='non-hvNLR'), alpha = 0.75, linewidth=0.5) +
    labs(colour = c('hvNLR', 'nonhvNLR')) +
  theme_classic() +
  xlim(0, 0.10)+
  theme(legend.position='right', text=element_text(size=10), legend.title=element_blank())+
  #labs(title='Empirical Distribution of Pi Values')+
  xlab(expression(pi))+
  ylab('count')+
  scale_color_manual(values=c('red', '#F8766D', '#00BFC4'))

both <- d_distro+pi_distro+plot_layout(guides = 'collect')
both
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 1134 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 590 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 2 rows containing missing values (`geom_bar()`).

![](popgen_permutations_files/figure-gfm/EVfig3-1.png)<!-- -->

``` r
ggsave(filename='C:\\Users\\chand\\Box Sync\\Krasileva_Lab\\Research\\chandler\\Krasileva Lab\\Outputs\\NLR Features Paper\\EMBO Submission\\Figure Panels\\EV_fig4a.png', plot=both, dpi=1000, width=180, height=70, unit='mm')
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

    ## Warning: Removed 1134 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 590 rows containing non-finite values (`stat_bin()`).

    ## Warning: Removed 2 rows containing missing values (`geom_bar()`).

``` r
#write source data 
write.csv(table %>% subset(select=c('Gene', 'HV', 'Pi', 'D')), file="./Source Data/EV Figure 3/EV 3AB/pi_d_genome.csv")
```

Repeat for the density of both distributions

``` r
#double double toil and trouble, two distros! 
q_D <- quantile(table$D, c(.05, .95), na.rm = TRUE)
q_Pi <- quantile(table$Pi, c(.05, .95), na.rm=TRUE)

observed_nlr <- table %>% filter(HV != 'all_genes') %>% filter(D > q_D[2]) %>% filter(Pi > q_Pi[2]) %>% nrow()

permutation_nlr = replicate(10000, {
  sample_small <- table[sample(nrow(table), nrow(table[table$HV != 'all_genes',]) , replace = FALSE), ] # sample the size of non-hvNLRs
  #sample_large <- meth[! rownames(meth) %in% rownames(sample_small), ] # sample the rest
  sample_small%>% filter(D > q_D[2]) %>% filter(Pi > q_Pi[2]) %>% nrow()  # get expected above 95%
})

print(paste('observed NLRs in top 5% of both distros:', observed_nlr))
```

    ## [1] "observed NLRs in top 5% of both distros: 1"

``` r
print(paste('expected NLRs in top 5% of both distros by chance:', mean(permutation_nlr)))
```

    ## [1] "expected NLRs in top 5% of both distros by chance: 1.7119"

``` r
print(paste('p value NLRs in top of both distros:', mean(observed_nlr < permutation_nlr)))
```

    ## [1] "p value NLRs in top of both distros: 0.5107"
