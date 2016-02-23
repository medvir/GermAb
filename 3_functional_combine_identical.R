#input data: _aligned.txt (fastq data were stitched, trimmed, collapsed into unique sequences, seqs <10 reads discarded, aligned to IMGT reference, relevant columns from aligned sam file extracted into _aligned_comb.txt)
#this script filters functional seqs and combines identical seqs with 0, 1 and 2 mutations from reference using cigar strings; output is _R_alleles_comb.txt, contains readcount, allele, number mut from reference, gene

# libaries
library(dplyr)
library(tidyr)

folder <- "/Volumes/data/AbX/germline/first_run_151008/"
filenames <- list.files(path = folder, pattern = "_aligned.txt")

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=F, col.names = c("reads", "allele", "start", "sequence", "mutations", "MD"))
  
  #clean data
  clean = data %>% 
    separate(mutations, into= c("a", "b", "mutations"), sep="\\:") %>%
    separate(MD, into= c("w", "x", "MD"), sep="\\:") %>%
    select(reads, allele, start, mutations, MD) %>%
    filter(!grepl("ORF", allele)) %>%
    filter(!grepl("P", allele)) %>%
    filter(grepl("IGHV", allele)) %>%
    group_by(allele) %>%
    arrange(allele)
  
  #combine 0  
  zero = clean %>%
    filter(mutations==0) %>%
    summarise(readcount=sum(reads)) %>%
    mutate(mutations=0) %>%
    mutate(mutations=as.numeric(mutations))
  
  #combine 1
  one = clean %>%
    filter(mutations==1) %>%
    mutate(MD_orig=MD) %>%
    separate(MD, into=c("first", "second"), sep="[A-Z]") %>%
    mutate(start=as.numeric(start)) %>%
    mutate(first=as.numeric(first)) %>%
    mutate(position=start+first) %>%
    mutate(nt=gsub("[0-9]", "", MD_orig)) %>%
    group_by(allele, position, nt) %>%
    mutate(mutations=as.numeric(mutations)) %>%
    group_by(mutations, allele, position, nt) %>%
    summarise(readcount=sum(reads))
  
  two = clean %>%
    filter(mutations==2) %>%
    mutate(MD_orig=MD) %>%
    separate(MD, into=c("first", "second", "third"), sep="[A-Z]") %>%
    mutate(start=as.numeric(start)) %>%
    mutate(first=as.numeric(first)) %>%
    mutate(second=as.numeric(second)) %>%
    mutate(position=start+first) %>%
    mutate(nt=gsub("[0-9]", "", MD_orig)) %>%
    group_by(allele, position, nt) %>%
    mutate(mutations=as.numeric(mutations)) %>%
    group_by(mutations, allele, position, nt) %>%
    summarise(readcount=sum(reads))
  
  #select >2 mut, then bind
  all = clean %>%
    filter(mutations>2) %>%
    rename(readcount=reads) %>%
    mutate(mutations=as.numeric(mutations)) %>%
    bind_rows(zero, one, two) %>%
    mutate(gene = sub("\\*.*", "", allele)) %>%
    select(readcount, allele, mutations, gene) %>%
    group_by(gene) %>%
    arrange(gene, desc(readcount)) #%>%
    #filter(readcount>25)
  
  #write output
    new.filename<-paste0("/Volumes/data/AbX/germline/first_run_151008/",substr(filenames[i],1,9),"_alleles_comb.txt")
    write.table(all,new.filename, sep = "\t",row.names=F)
}
  
  