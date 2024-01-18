require(entropy)
require(msa)
require(tidyverse)
require(parallel)
install.packages("parallel")
getwd()

Atha_Common <- read_delim("~/Dropbox/NLRomes/Atha_SnkMk_NLRome/Atha_SnkMk_GeneTable.txt",delim = "\t")
atha_clades <- Common %>% filter(!is.na(File))%>%group_by(Clade,File,HV) %>%count%>%ungroup
clade <- "Int10172_376_423_L_90_128_L_41"
get_atha_ent <- function(clade){
  ali <- paste0("~/Dropbox/NLRomes/Atha_SnkMk_NLRome/", (atha_clades %>% filter(Clade == clade))$File)
  maa <- readAAMultipleAlignment(ali)
  seqs <- tibble(Names = maa@unmasked@ranges@NAMES)
  entlist <- list()
  for(jj in seq_along(seqs$Names)){
    seq <- seqs$Names[[jj]]
    ent <- get_ent(maa,seq) %>% unlist()
    entlist[[jj]] <- tibble(Name = seq %>%str_remove(" "),Ent = ent,Clade = clade,Pos = 1:length(ent))
  }
  return(bind_rows(entlist))
}
get_atha_ent(clade)

atha_clades <- atha_clades %>% filter(!is.na(File))%>%ungroup


atha_clades_ent <- mclapply(atha_clades$Clade,get_atha_ent,mc.cores = 8)
Atha_Ent <- bind_rows(atha_clades_ent)
Atha_Ent
write_delim(Atha_Ent,"Athaliana_NLR_Entropy.tsv",delim = "\t")
Atha_Ent <- read_delim("Athaliana_NLR_Entropy.tsv")

Atha_Common <- Common
Atha_per_seq <- Atha_Ent %>% group_by(Name) %>% summarise(mean = mean(Ent), 
                                                          median = median(Ent),
                                                          top_10 = quantile(Ent,probs = 0.9,na.rm=T),
                                                          top_5 = quantile(Ent,probs = 0.95,na.rm=T),
                                                          top_2 = quantile(Ent,probs = 0.98,na.rm=T),
                                                          top_1 = quantile(Ent,probs = 0.99,na.rm=T)
                                                          
)
Atha_per_seq <- left_join(Atha_per_seq,Atha_Common %>% transmute(Name = Gene, HV = HV))

Atha_per_seq %>% ggplot(aes(x=mean,color = as.factor(HV)))+geom_histogram()
Atha_per_seq %>% ggplot(aes(x=median,color = as.factor(HV)))+geom_histogram()
Atha_per_seq %>% ggplot(aes(x=top_10,color = as.factor(HV)))+geom_histogram()
Atha_per_seq %>% ggplot(aes(x=top_5,color = as.factor(HV)))+geom_histogram()
Atha_per_seq %>% ggplot(aes(x=top_2,color = as.factor(HV)))+geom_histogram()
Atha_per_seq %>% ggplot(aes(x=top_1,color = as.factor(HV)))+geom_histogram()

HighEntnonHV <- Atha_per_seq %>%filter(mean>0.3,HV==0)
Atha_Common %>% filter(Gene %in% HighEntnonHV$Name) %>%count(by=Clade)

atha_pos_ent <- Atha_Ent %>% filter(!is.na(Ent))%>%group_by(Pos) %>% summarise(mean = mean(Ent))
ggplot(atha_pos_ent,aes(x=Pos,y=mean))+geom_point()+ylim(0,0.5)

atha_r10_ent <- Atha_Ent %>% group_by(Name) %>% slice_max(Ent,n=10) %>% slice_min(Ent,n=1,with_ties = FALSE) %>% 
  ungroup %>% mutate(Name = Name%>%str_remove("MRNA.*"))
atha_r10_ent <- left_join(atha_r10_ent,Atha_Common %>% transmute(Name = Gene, HV = HV))

atha_r10_ent %>% ggplot(aes(x=Ent,
                            fill = as.factor(-HV)
                            
))+
  geom_histogram(color = "black") + 
  xlab("Entropy (bits)")+
  ylab("Number of NLRs")+
  theme_classic() +
  geom_vline(xintercept = 1.5,linetype="dotted")+
  scale_fill_discrete(name ="", labels = c("HV", "non-HV"))
