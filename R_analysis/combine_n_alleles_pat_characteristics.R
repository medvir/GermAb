# libaries
library(dplyr)
library(tidyr)

data1 = read.table("/Volumes/data/AbX/germline/first_run_151008/all_results.txt", header=F, sep = " ", col.names=c("patient_ID", "gene", "n_alleles"))
data2 = read.table("/Volumes/data/AbX/germline/160122/all_results.txt", header=F, sep = " ", col.names=c("patient_ID", "gene", "n_alleles"))
data3 = read.table("/Volumes/data/AbX/germline/160624/all_results.txt", header=F, sep = " ", col.names=c("patient_ID", "gene", "n_alleles"))
patients = read.table("/Volumes/data/AbX/germline/patient_characteristics.txt", header=T, sep="\t")

data1 = data1 %>%
  mutate(run="1") %>%
  filter(!grepl("46179_S1_|Hy|HD|AK170|41895|17420|26500|26586|34545|42335", patient_ID))

data2 = data2 %>%
  mutate(run="2") %>%
  filter(!grepl("17811|18322|15504|18357|15224|18669|18311|18418|19138|13853|25478|17241|41895|31822", patient_ID))

data3 = data3 %>%
  mutate(run="3") %>%
  filter(!grepl("4-59-|4-28-|mix|31396", patient_ID))

#reformat data long to wide and back to introduce 0`s for missing alleles
data12 = bind_rows(data1, data2, data3) %>%
  separate(patient_ID, into=c("patient", "ID"), sep="\\_") %>%
  select(patient, gene, n_alleles, run) %>%
  mutate(patient=as.integer(patient)) %>%
  spread(gene, n_alleles) %>%
  replace(is.na(.), 0) %>%
  gather("gene", "n_alleles", 3:52)

#join data with patient info
number_alleles = full_join(data12,patients, by="patient") %>%
  filter(!is.na(run))

write.table(number_alleles, "/Volumes/data/AbX/germline/results/patients_n_alleles_ethn_neut_subtype.txt", row.names=F, quote=F)
#after this step check alleles>5 and correct if necessary, save as patients_n_alleles_corr_ethn_neut_subtype.txt

#reformat to create "wide" dataset by patient (eg needed for heatmap) and combine all relevant patient information in one column
data4=read.table("/Volumes/data/AbX/germline/results/patients_n_alleles_ethn_neut_subtype.txt", header=T, sep=" ")
data_wide = data4 %>%
  mutate(n_alleles = replace(n_alleles, n_alleles>4, NA)) %>%
  spread(gene, n_alleles) %>%
  unite(patient_ID, patient, ethnicity, bnAb, type, run) 

write.table(data_wide, "/Volumes/data/AbX/germline/results/patients_number_alleles_4_wide.txt", row.names=F, quote=F)
