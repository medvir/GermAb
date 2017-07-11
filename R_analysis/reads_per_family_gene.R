#analyze readcounts per VH family and per VH gene for both MiSeq runs

# libaries
library(dplyr)
library(tidyr)
library(ggplot2)

#data 1st run
folder <- "/Volumes/data/AbX/germline/first_run_151008/"
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
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
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
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

#data 3rd run
folder <- "/Volumes/data/AbX/germline/160624/"
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
out.gene3<-""
out.family3<-""

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
  
  out.gene3<-rbind(out.gene3, readcount_gene)
  out.family3<-rbind(out.family3, readcount_family)
}

out.gene3=out.gene3[-1,] %>%
  mutate(reads_gene=as.numeric(reads_gene)) %>%
  filter(!grepl("4-59-", patient)) %>%
  filter(!grepl("4-28-", patient)) %>%
  filter(!grepl("mix", patient)) %>%
  mutate(run="3")

out.family3=out.family3[-1,] %>%
  mutate(reads_family=as.numeric(reads_family)) %>%
  filter(!grepl("4-59-", patient)) %>%
  filter(!grepl("4-28-", patient)) %>%
  filter(!grepl("mix", patient)) %>%
  mutate(run="3")

#combine data
all.gene=bind_rows(out.gene1, out.gene2, out.gene3)
all.family=bind_rows(out.family1, out.family2, out.family3)

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
  gather("gene", "reads_gene", 3:52) %>%
  mutate(patient=as.integer(patient)) %>%
  arrange(patient)

#write output
write.table(gene_long, "/Volumes/data/AbX/germline/results/reads_per_gene.txt", row.names=F, quote=F)
write.table(family_long, "/Volumes/data/AbX/germline/results/reads_per_family.txt", row.names=F, quote=F)

#figures (without 0`s, log scale)
fig_genes <- ggplot(all.gene, aes(x=gene, y=reads_gene, color=run)) +
  geom_boxplot() +
  xlab(" \nGene") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1, size=10), axis.text.y=element_text(size=16, color="black"), axis.title = element_text(size = 18), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(),
        legend.title=element_text(size=14), legend.text=element_text(size=14), axis.line=element_line(size=0.5)) +
  ylab("# Reads") +
  ggtitle("Reads per Gene") +
  scale_x_discrete(limits=c("", "IGHV6-1", "IGHV1-2", "IGHV1-3", "IGHV4-4", "IGHV7-4-1", "IGHV2-5", "IGHV3-7", "IGHV3-64D", "IGHV1-8", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV1-18", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV1-24", "IGHV2-26", "IGHV4-28", "IGHV3-30", "IGHV4-30-2", "IGHV3-30-3", "IGHV4-30-4", "IGHV3-30-5", "IGHV4-31", "IGHV3-33", "IGHV4-34", "IGHV4-38-2", "IGHV3-43D", "IGHV4-39", "IGHV3-43", "IGHV1-45", "IGHV1-46", "IGHV3-48", "IGHV3-49", "IGHV5-10-1", "IGHV5-51", "IGHV3-53", "IGHV1-58", "IGHV4-59", "IGHV4-61", "IGHV3-64", "IGHV3-66", "IGHV1-69", "IGHV1-69-2", "IGHV2-70", "IGHV3-72", "IGHV3-73", "IGHV3-74", "IGHV1-2")) +
  scale_y_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE), breaks=c(10,100,1000,10000)) +
  scale_color_discrete(name="Run")

fig_genes
ggsave(filename="/Volumes/data/AbX/germline/results/figures/Readspergene_log_order.pdf", plot=fig_genes, width = 30/2.54, height = 14/2.54)

fig_family <- ggplot(all.family, aes(x=family, y=reads_family, color=run)) +
  geom_boxplot() +
  xlab(" \nVH Family") +
  theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1, size=16, color="black"), axis.text.y=element_text(size=16, color="black"), axis.title = element_text(size = 18), 
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(),
        legend.title=element_text(size=14), legend.text=element_text(size=14), axis.line=element_line(size=0.5)) +
  ylab("# Reads") +
  ggtitle("Reads per VH-Family") +
  scale_y_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE), breaks=c(10,100,1000,10000)) +
  scale_color_discrete(name="Run")
  
  fig_family
  ggsave(filename="/Volumes/data/AbX/germline/results/figures/Readsperfamily_log.pdf", plot=fig_family, width = 30/2.54, height = 16/2.54)
  
#determine number of patients with 0 reads for VH families
  missing_family=family_long %>%
    filter(reads_family==0) %>%
    group_by(run, family) %>%
    summarise(n=n()) 
  
  fig_missing <- ggplot(missing_family, aes(x=family, y=n, fill=run)) +
    geom_bar(stat="identity", position="dodge") +
    xlab(" \nVH Family") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
          axis.text.x=element_text(size=12), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
    ylab("# Patients") +
    ggtitle("Number of patients with missing VH family reads")
  
  fig_missing
  ggsave(filename="/Volumes/data/AbX/germline/results/figures/Patients_missing_families.pdf", plot=fig_missing, width = 30/2.54, height = 21/2.54)
  
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
  ggsave(filename="/Volumes/data/AbX/germline/results/figures/Patients_missing_genes.pdf", plot=fig_missing_gene, width = 30/2.54, height = 21/2.54)
  
#exclude patients with <10000 reads
  gene_excl = gene_long %>%
    filter(!grepl("17811|18322|15504|18357|15224|18669|18311|18418|19138|13853|25478|17241", patient))
  
  write.table(gene_excl, "/Volumes/data/AbX/germline/results/reads_per_gene_10000reads.txt", row.names=F, quote=F)
  
  family_excl = family_long %>%
    filter(!grepl("17811|18322|15504|18357|15224|18669|18311|18418|19138|13853|25478|17241", patient))
  
  write.table(family_excl, "/Volumes/data/AbX/germline/results/reads_per_family_10000reads.txt", row.names=F, quote=F)
  
  #determine number of patients with 0 reads for genes using only patients with >10000 reads
  missing_gene_excl=gene_excl %>%
    filter(reads_gene==0) %>%
    group_by(run, gene) %>%
    summarise(n=n()) 
  
  fig_missing_gene_excl <- ggplot(missing_gene_excl, aes(x=gene, y=n, fill=run)) +
    geom_bar(stat="identity", position="dodge") +
    xlab(" \nGene") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
          axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
    ylab("# Patients") +
    ggtitle("Number of patients (>10000 reads only) with missing Genes")
  
  fig_missing_gene_excl
  ggsave(filename="/Volumes/data/AbX/germline/results/figures/Patients_missing_genes_10000reads.pdf", plot=fig_missing_gene_excl, width = 30/2.54, height = 21/2.54)
  
  #determine number of patients with 0 reads for VH families using only patients >10000 reads
  missing_family_excl=family_excl %>%
    filter(reads_family==0) %>%
    group_by(run, family) %>%
    summarise(n=n()) 
  
  fig_missing_family_excl <- ggplot(missing_family_excl, aes(x=family, y=n, fill=run)) +
    geom_bar(stat="identity", position="dodge") +
    xlab(" \nVH Family") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=20), 
          axis.text.x=element_text(size=8), axis.text.x=element_text(size=20), panel.grid.major = element_line(colour = "grey")) +
    ylab("# Patients\n") +
    ggtitle("Number of patients (>10000 reads) with missing VH family reads")
  
  fig_missing_family_excl
  ggsave(filename="/Volumes/data/AbX/germline/results/figures/Patients_missing_families_10000reads.pdf", plot=fig_missing_family_excl, width = 30/2.54, height = 21/2.54)
  \n
  
  #analyze readount by ethnicity
  patients = read.table("/Volumes/data/AbX/germline/patient_characteristics.txt", header=T, sep="\t")
  
  readcount_ethn = full_join(gene_long,patients, by="patient") %>%
    select(patient, run, gene, reads_gene, ethnicity) %>%
    filter(run!="2") %>%
    filter(ethnicity!="other") %>%
    filter(complete.cases(.))
  
  fig_ethn_count <- ggplot(readcount_ethn, aes(x=gene, y=reads_gene, color=ethnicity)) +
    geom_boxplot() +
    xlab(" \nGene") +
    theme(plot.title=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1, size=10), axis.text.y=element_text(size=16, color="black"), axis.title = element_text(size = 18), 
          panel.grid.major = element_line(colour = "grey", linetype = "dotted"), panel.background = element_blank(),
          legend.title=element_text(size=14), legend.text=element_text(size=14), axis.line=element_line(size=0.5)) +
    ylab("# Reads") +
    ggtitle("Reads per Gene") +
    #scale_x_discrete(limits=c("", "IGHV6-1", "IGHV1-2", "IGHV1-3", "IGHV4-4", "IGHV7-4-1", "IGHV2-5", "IGHV3-7", "IGHV3-64D", "IGHV1-8", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV1-18", "IGHV3-20", "IGHV3-21", "IGHV3-23", "IGHV1-24", "IGHV2-26", "IGHV4-28", "IGHV3-30", "IGHV4-30-2", "IGHV3-30-3", "IGHV4-30-4", "IGHV3-30-5", "IGHV4-31", "IGHV3-33", "IGHV4-34", "IGHV4-38-2", "IGHV3-43D", "IGHV4-39", "IGHV3-43", "IGHV1-45", "IGHV1-46", "IGHV3-48", "IGHV3-49", "IGHV5-10-1", "IGHV5-51", "IGHV3-53", "IGHV1-58", "IGHV4-59", "IGHV4-61", "IGHV3-64", "IGHV3-66", "IGHV1-69", "IGHV1-69-2", "IGHV2-70", "IGHV3-72", "IGHV3-73", "IGHV3-74", "IGHV1-2")) +
    scale_y_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE), breaks=c(10,100,1000,10000)) +
    scale_color_discrete(name="Run")
  
  fig_ethn_count
  