# libaries
library(dplyr)
library(tidyr)

data1 = read.table("/Volumes/data/AbX/germline/151126_combine_mut/all_results_comb.txt", header=F, sep = " ", col.names=c("patient_ID", "gene", "n_alleles"))
data2 = read.table("/Volumes/data/AbX/germline/160122/all_results.txt", header=F, sep = " ", col.names=c("patient_ID", "gene", "n_alleles"))
patients = read.table("/Volumes/data/AbX/germline/patient_characteristics.txt", header=T, sep="\t")

data1 = data1 %>%
  mutate(run="1") %>%
  filter(!grepl("46179_S1_", patient_ID))

data2 = data2 %>%
  mutate(run="2")

data12 = bind_rows(data1, data2) %>%
  separate(patient_ID, into=c("patient", "ID"), sep="\\_") %>%
  select(patient, gene, n_alleles, run) %>%
  mutate(patient=as.integer(patient))

number_alleles = full_join(data12,patients, by="patient")

write.table(number_alleles, "/Volumes/data/AbX/germline/results/patients_n_alleles_ethn_neut_subtype.txt", row.names=F, quote=F)
