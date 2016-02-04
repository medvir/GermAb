#analyze readcounts per VH family and per VH gene for both MiSeq runs

# libaries
library(dplyr)
library(tidyr)
library(ggplot2)

#data 1st run
folder <- "/Volumes/data/AbX/germline/151126_combine_mut/"
filenames <- list.files(path = folder, pattern = "_R_alleles_comb.txt")
out.gene1<-""
out.family1<-""

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=T)
  
  readcount_gene = data %>%
    group_by(gene) %>%
    summarise(reads_gene=sum(readcount)) %>%
    mutate(patient=substr(filenames[i],1,5)) 
  
  readcount_family = data %>%
    separate(gene, into=c("family", "b"), sep="\\-") %>%
    group_by(family) %>%
    summarise(reads_family=sum(readcount)) %>%
    mutate(patient=substr(filenames[i],1,5))
  
  out.gene1<-rbind(out.gene1, readcount_gene)
  out.family1<-rbind(out.family1, readcount_family)
}

out.gene1=out.gene1[-1,] %>%
  mutate(reads_gene=as.numeric(reads_gene)) %>%
  mutate(run="1")

out.family1=out.family1[-1,] %>%
  mutate(reads_family=as.numeric(reads_family)) %>%
  mutate(run="1")

#data 2nd run
folder <- "/Volumes/data/AbX/germline/160122/"
filenames <- list.files(path = folder, pattern = "_R_alleles_comb.txt")
out.gene2<-""
out.family2<-""

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=T)

readcount_gene = data %>%
  group_by(gene) %>%
  summarise(reads_gene=sum(readcount)) %>%
  mutate(patient=substr(filenames[i],1,5)) 

readcount_family = data %>%
  separate(gene, into=c("family", "b"), sep="\\-") %>%
  group_by(family) %>%
  summarise(reads_family=sum(readcount)) %>%
  mutate(patient=substr(filenames[i],1,5))
  
out.gene2<-rbind(out.gene2, readcount_gene)
out.family2<-rbind(out.family2, readcount_family)
}

out.gene2=out.gene2[-1,] %>%
  mutate(reads_gene=as.numeric(reads_gene)) %>%
  mutate(run="2")

out.family2=out.family2[-1,] %>%
  mutate(reads_family=as.numeric(reads_family)) %>%
  mutate(run="2")

#combine data
all.gene=bind_rows(out.gene1, out.gene2)
all.family=bind_rows(out.family1, out.family2)

#introduce 0 for families that are not represented
family_wide=all.family %>%
  spread(family, reads_family) %>%
  replace(is.na(.), 0)

family_long=family_wide %>%
  gather("family", "reads_family", 3:9) %>%
  arrange(patient)

#introduce 0 for genes that are not represented
gene_wide=all.gene %>%
  spread(gene, reads_gene) %>%
  replace(is.na(.), 0)

gene_long=gene_wide %>%
  gather("gene", "reads_gene", 3:55) %>%
  arrange(patient)

#write output
write.table(gene_long, "/Volumes/data/AbX/germline/results/reads_per_gene.txt", row.names=F, quote=F)
write.table(family_long, "/Volumes/data/AbX/germline/results/reads_per_family.txt", row.names=F, quote=F)

#figures (without 0`s, log scale)
fig_genes <- ggplot(all.gene, aes(x=gene, y=reads_gene, color=run)) +
  geom_boxplot() +
  xlab(" \nGene") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
        axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
  ylab("# Reads") +
  ggtitle("Reads (>0) per Gene") +
  scale_y_log10()

fig_genes
ggsave(filename="/Volumes/data/AbX/germline/results/Readspergene_log.pdf", plot=fig_genes, width = 30/2.54, height = 21/2.54)

fig_family <- ggplot(all.family, aes(x=family, y=reads_family, color=run)) +
  geom_boxplot() +
  xlab(" \nGene") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
        axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
  ylab("# Reads") +
  ggtitle("Reads (>0) per VH-Family") +
    scale_y_log10()
  
  fig_family
  ggsave(filename="/Volumes/data/AbX/germline/results/Readsperfamily_log.pdf", plot=fig_family, width = 30/2.54, height = 21/2.54)
  
#determine number of patients with 0 reads for VH families
  missing_family=family_long %>%
    filter(reads_family==0) %>%
    group_by(run, family) %>%
    summarise(n=n()) 
  
  fig_missing <- ggplot(missing_family, aes(x=family, y=n, fill=run)) +
    geom_bar(stat="identity", position="dodge") +
    xlab(" \nVH Family") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
          axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
    ylab("# Patients") +
    ggtitle("Number of patients with missing VH family reads")
  
  fig_missing
  ggsave(filename="/Volumes/data/AbX/germline/results/Patients_incomplete_families.pdf", plot=fig_missing, width = 30/2.54, height = 21/2.54)
  
  #determine number of patients with 0 reads for genes
  missing_gene=gene_long %>%
    filter(reads_gene==0) %>%
    group_by(run, gene) %>%
    summarise(n=n()) 
  
  fig_missing_gene <- ggplot(missing_gene, aes(x=gene, y=n, fill=run)) +
    geom_bar(stat="identity", position="dodge") +
    xlab(" \nGene") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
          axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
    ylab("# Patients") +
    ggtitle("Number of patients with missing Genes")
  
  fig_missing_gene
  ggsave(filename="/Volumes/data/AbX/germline/results/Patients_incomplete_genes.pdf", plot=fig_missing_gene, width = 30/2.54, height = 21/2.54)
  