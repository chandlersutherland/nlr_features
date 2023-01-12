setwd("/Users/prigozhin/Dropbox/Shared_Files/Daniil_William_Share")
getwd()
require(tidyverse)

gff <- read_delim("/Users/prigozhin/Library/CloudStorage/Box-Box/hvNLRs-nt/GFF/TAIR10_GFF3_genes.gff", delim = "\t", col_names = c("chr","source","type","start","end","score","strand","phase","attributes"))
genes <- gff %>% filter(type == "gene")
genes <- genes %>% mutate(Gene = attributes %>% str_remove(".*Name="))

### Add distance to the nearest transposon--------
te_data <- read_delim("/Users/prigozhin/Library/CloudStorage/Box-Box/hvNLRs-nt/GFF/TAIR10_Transposable_Elements.txt", delim = "\t", col_names = T)
te_data <- te_data %>% mutate(chr = paste0("Chr",Transposon_Name %>% str_remove("^AT")%>% str_remove("TE.*")))%>%
  mutate(source = "TAIR10", type = "TE", start =Transposon_min_Start, end =Transposon_max_End, strand = ifelse(orientation_is_5prime, "+","-"))%>%
  mutate(name = Transposon_Name, family = Transposon_Family, superfamily = Transposon_Super_Family)%>%
  select(chr, source, type, start, end, strand ,name  ,     family,     superfamily)
genes
te_data

find_te_distance <- function(name){
  line <- genes %>% filter(Gene==name)  
  if (nrow(line) != 1){stop("Gene name not found, or too many hits in gff")}
  gene_start <- line$start %>% unlist()
  gene_end <- line$end %>% unlist()
  gene_chr <- line$chr %>% unlist()
  candidates <- te_data %>% filter (chr == gene_chr)
  overlaps <- candidates %>% filter(end >= gene_start & start <= gene_end)
  if (nrow(overlaps)>0){return(0)}else{
    return(min(abs(gene_start-candidates$start),
               abs(gene_start-candidates$end),
               abs(gene_end-candidates$start),
               abs(gene_end-candidates$end)    ))
  }  
}
name<-"AT3G07040"
find_te_distance("AT3G07040")

genes <- genes %>% mutate(te_dist = mapply(genes$Gene,FUN = find_te_distance))
write_delim(genes,"Atha_genes_with_TE_dist.tsv",delim = "\t")
