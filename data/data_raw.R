library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

rm(list = ls())
load("result_data.rda")
Yulab_data <- subset(dataset, !(target %in% c("Unspecified", "n.a.")))
save(Yulab_data, file = "./Yulab_data.rda")


load("herb_data.rda")
herb_data <- dataset %>%
  select("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name", "molecule", "pubchem_id", 
         "PubChem_CID", "target")
save(herb_data, file = "./Yulab_data2.rda")
