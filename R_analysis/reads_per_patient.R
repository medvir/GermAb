#determines total readcount per patient

# libaries
library(dplyr)
library(tidyr)
library(ggplot2)

#get readcounts from 2nd run 
folder <- "/Volumes/data/AbX/germline/160122/"
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
out.file2<-""

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=T)

data= data %>%
  mutate(patient=substr(filenames[i],1,5)) %>%
  select(readcount, patient)

out.file2<-rbind(out.file2, data)
}

#get readcounts from 1st run 
folder <- "/Volumes/data/AbX/germline/first_run_151008/"
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
out.file1<-""

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=T)
  
  data= data %>%
    mutate(patient=substr(filenames[i],1,5)) %>%
    select(readcount, patient)
  
  out.file1<-rbind(out.file1, data)
}

#get readcounts from 3rd run 
folder <- "/Volumes/data/AbX/germline/160624/"
filenames <- list.files(path = folder, pattern = "_alleles_comb.txt")
out.file3<-""

for (i in 1:length(filenames))
{
  original.file<-paste(folder,filenames[i],sep="")
  data <- read.delim(original.file, header=T)
  
  data= data %>%
    mutate(patient=substr(filenames[i],1,5)) %>%
    select(readcount, patient)
  
  out.file3<-rbind(out.file3, data)
}

#clean 1st run
reads1_patient <- out.file1[-1,] %>%
  group_by(patient) %>%
  filter(!grepl("HD", patient)) %>%
  filter(!grepl("Hy", patient)) %>%
  mutate(readcount=as.numeric(readcount)) %>%
  summarise(reads=sum(readcount)) %>%
  mutate(run="1")
 
 #clean 2nd run
 reads2_patient <- out.file2[-1,] %>%
   group_by(patient) %>%
   mutate(readcount=as.numeric(readcount)) %>%
   summarise(reads=sum(readcount)) %>%
   mutate(run="2")
 
 #clean 3rd run
 reads3_patient <- out.file3[-1,] %>%
   group_by(patient) %>%
   mutate(readcount=as.numeric(readcount)) %>%
   summarise(reads=sum(readcount)) %>%
   mutate(run="3")
 
 #combine data 
 all=bind_rows(reads1_patient, reads2_patient, reads3_patient)
 write.table(all, "/Volumes/data/AbX/germline/results/reads_per_patient.txt", row.names=F, quote=F)
 
 fig = ggplot(all, aes(x=patient, y=reads, color=run)) +
   geom_point() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
   scale_y_continuous(limits=c(0,150000)) +
   xlab("Sample") +
   ylab("Reads") +
   ggtitle("Reads per Patient")
 
 fig
 ggsave(filename="/Volumes/data/AbX/germline/results/figures/reads_per_patient.pdf", plot=fig, width = 30/2.54, height = 21/2.54)
 