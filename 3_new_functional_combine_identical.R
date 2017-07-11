library(tidyr)
library(dplyr)

folder <- "/data/AbX/germline/160624/"
filenames <- list.files(path = folder, pattern = "_aligned.txt")

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=F, col.names = c("readcount", "allele", "start", "sequence", "mutations", "MD", "XA"))

in_data <- data %>%
  separate(mutations, into= c("a", "b", "mutations"), sep="\\:") %>%
  separate(MD, into= c("w", "x", "MD"), sep="\\:") %>%
  select(readcount, allele, sequence, start, mutations, MD) %>%
  filter(!grepl("ORF", allele)) %>%
  filter(!grepl("P", allele)) %>%
  filter(grepl("IGHV", allele)) %>%
  group_by(allele) %>%
  arrange(allele) %>%
  mutate(mutations=as.numeric(mutations))

collapse_short_sequences = function(df) {
    # function to collapse shorter sequences contained within identical longer ones
    # input is a data frame with the first three colums containing read count, allele and nucleotide sequence
    # sums read count, keeps allele and all other columns of the longer sequence, deletes entire row of shorter sequence

    colnames(df) = c("readcount", "allele", "sequence", "start", "mutations", "MD")
    df$sequence = as.character(df$sequence)
    df = df %>% mutate(length = nchar(sequence)) %>% arrange(desc(length))

    for (i in nrow(df):2) {
        for (j in (i-1):1) {
            if (grepl(df[i,3], df[j,3])) {
                df[j,1] = df[j,1] + df[i,1]
                df = df[-i,]
                break
            } 
        }
    }
    df = df %>% select(-matches("length")) %>% arrange(desc(readcount))
    return(df)
}

out_data_bwaL7 = collapse_short_sequences(in_data) %>%
  mutate(mutations=as.numeric(mutations))

#combine 0  
zero = out_data_bwaL7 %>%
  filter(mutations==0) %>%
  group_by(allele) %>%
  summarise(readcount=sum(readcount)) %>%
  mutate(mutations=0) %>%
  mutate(mutations=as.numeric(mutations))

#combine 1
one = out_data_bwaL7 %>%
  filter(mutations==1) %>%
  mutate(MD_orig=MD) %>%
  separate(MD, into=c("first", "second"), sep="[A-Z]") %>%
  mutate(start=as.numeric(start)) %>%
  mutate(first=as.numeric(first)) %>%
  mutate(position=start+first) %>%
  mutate(nt=gsub("[0-9]", "", MD_orig)) %>%
  group_by(mutations, allele, position, nt) %>%
  summarise(readcount=sum(readcount)) %>%      #the following deletes reads with mutation at position 229 (or 226, depends on primers) (wt: CCAAGAACCAGTT, mut: CCAAGACCCAGTT)
  filter(position != 230 | !grepl("IGHV4", allele) | !grepl("A", nt)) %>%
  filter(position != 227 | !grepl("IGHV4", allele) | !grepl("A", nt)) %>%
  filter(position != 233 | !grepl("IGHV4", allele) | !grepl("A", nt)) %>%
  filter(position != 40 | !grepl("IGHV4", allele) | !grepl("C", nt)) 

two = out_data_bwaL7 %>%
  filter(mutations==2) %>%
  mutate(MD_orig=MD) %>%
  separate(MD, into=c("first", "second", "third"), sep="[A-Z]") %>%
  mutate(start=as.numeric(start)) %>%
  mutate(first=as.numeric(first)) %>%
  mutate(second=as.numeric(second)) %>%
  mutate(position=start+first) %>%
  mutate(nt=gsub("[0-9]", "", MD_orig)) %>%
  group_by(mutations, allele, position, nt) %>%
  summarise(readcount=sum(readcount))

#select >2 mut, then bind
all = out_data_bwaL7 %>%
  filter(mutations>2) %>%
  bind_rows(zero, one, two) %>%
  mutate(gene = sub("\\*.*", "", allele)) %>%
  select(readcount, allele, mutations, gene) %>%
  group_by(gene) %>%
  arrange(gene, desc(readcount))

#write output
new.filename<-paste0("/data/AbX/germline/160624/",substr(filenames[i],1,9),"_alleles_comb.txt")
write.table(all,new.filename, sep = "\t",row.names=F, quote=F)
}



