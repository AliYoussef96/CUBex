## ----eval=F---------------------------------------------------------------
## if (!requireNamespace("devtools", quietly=TRUE)){
##         install.packages("devtools")}
## devtools::install_github('AliYoussef96/vhcub')
## 


## ----eval=F---------------------------------------------------------------
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("coRdon")
## BiocManager::install("Biostrings")
## 
## 
## install.packages("vhcub")
## 


## ----eval=F---------------------------------------------------------------
## install.packages(ggplot2)
## 


## -------------------------------------------------------------------------
library("vhcub")
library("ggplot2")

fasta <- fasta.read("Ecoli.fasta","Ecoli.fasta")
fasta.org <- fasta[[1]]


## -------------------------------------------------------------------------
enc.df <- ENc.values(fasta.org)
head(enc.df, 10)



## -------------------------------------------------------------------------
enc.df$index <- row.names(enc.df)
ggplot(data = enc.df , aes(x = index, y = ENc)) + 
  geom_point(size = 2) + geom_hline(yintercept=35, color = "red")


## -------------------------------------------------------------------------
pa <- read.csv("E.coli_k12.csv")
pa_sorted <- pa[with(pa,order(-abundance)),]
nrow(pa_sorted) * 0.05
ref_pa <- pa_sorted[1:205,]$string_external_id
ref_pa


## -------------------------------------------------------------------------
ref_genes <- data.frame()
for (i_ref_pa in ref_pa){
  df <- fasta.org[grepl(i_ref_pa, fasta.org$seq_name, fixed = T),]
  ref_genes <- rbind(df,ref_genes)
}

enc.ref.gene <- ENc.values(ref_genes)


## -------------------------------------------------------------------------

cai.df <- CAI.values(fasta.org, enc.ref.gene, fasta.org, set.len = 100 ,genetic.code="1")
cai.df <- cai.df[cai.df$CAI <= 1,]

head(cai.df, 10)


## -------------------------------------------------------------------------
ggplot(data = cai.df , aes (y = CAI)) + geom_boxplot()


## -------------------------------------------------------------------------
id.ext <- function(id){
  id.new <- regmatches(id,regexec("b....",id))[[1]]
  return(id.new)
}
id.new  <- lapply(cai.df$gene.name, id.ext)
cai.df$string_external_id <- id.new

pa.cai <- merge(pa, cai.df[2:3] , by= "string_external_id")

ggplot(data = pa.cai, aes(x = CAI , y = abundance))+
 geom_point(size = 2)+ geom_smooth()


## -------------------------------------------------------------------------
pa.cai$abundance <- log2(pa.cai$abundance + 1)
ggplot(data = pa.cai, aes(x = CAI , y = abundance))+
 geom_point(size = 2)+ geom_smooth()


## -------------------------------------------------------------------------
rscu.df <- RSCU.values(fasta.org) 


## -------------------------------------------------------------------------

rscu.mean <- data.frame(apply(rscu.df, 2, mean))
colnames(rscu.mean) <- "rscu"


## -------------------------------------------------------------------------
rscu.mean$codons <- row.names(rscu.mean)
ggplot(data= rscu.mean , aes(x = reorder(codons, rscu), y = rscu))+
  geom_bar(stat='identity') + geom_hline(yintercept=1, color = "red")


## -------------------------------------------------------------------------
rscu.ref.genes <- RSCU.values(ref_genes) 

rscu.ref.genes.mean <- data.frame(apply(rscu.ref.genes, 2, mean))
colnames(rscu.ref.genes.mean) <- "rscu"
rscu.ref.genes.mean$codons <- row.names(rscu.ref.genes.mean)


ggplot(data= rscu.ref.genes.mean , aes(x = reorder(codons, rscu), y = rscu))+
  geom_bar(stat='identity') + geom_hline(yintercept=1, color = "red")


## -------------------------------------------------------------------------
library(Biostrings)
rscu.ref.genes.mean$codons <- toupper(rscu.ref.genes.mean$codons)

codon.to.aa <- function(codon){
  return(GENETIC_CODE[[codon]])
}

aa <- lapply(rscu.ref.genes.mean$codons ,codon.to.aa)
aa <- unlist(aa)
rscu.ref.genes.mean$aa <- aa


## -------------------------------------------------------------------------
library(dplyr)

p.codons <- rscu.ref.genes.mean %>% group_by(aa) %>% filter(rscu == max(rscu))
p.codons<- p.codons %>% filter(!aa %in% c("*","M","W"))
p.codons


## -------------------------------------------------------------------------

ggplot(data= p.codons , aes(x = reorder(codons, rscu), y = rscu))+
  geom_bar(stat='identity') + geom_hline(yintercept=1, color = "red")


## -------------------------------------------------------------------------

p.codons <- p.codons %>% filter(rscu >= 1)
ggplot(data= p.codons , aes(x = reorder(codons, rscu), y = rscu))+
  geom_bar(stat='identity') + geom_hline(yintercept=1, color = "red")

