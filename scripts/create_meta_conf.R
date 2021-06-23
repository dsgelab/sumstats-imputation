rm(list = ls())

library(rjson)
library(data.table)
library(dplyr)
library(plyr)

# conf files located in "gs://covid19-hg-analysis/20201215/conf/2/json/B2_ALL.json"
setwd('~')

f <- 'Projects/covid19-hgi/META_ANALYSIS/20201215_conf_5_filtered_json_B2_ALL.json'
d <- fromJSON(file = f)

d <- d$meta

df <- as.data.frame(do.call(rbind, lapply(d, unlist)), stringsAsFactors = F)

df <- df %>%
  mutate(ana = case_when(grepl("ANA(A1|_A1)", toupper(file)) ~ "A1",
                         (grepl("ANA(A2|_A2)", toupper(file)) | grepl(".A2", toupper(file))) ~ "A2",
                         (grepl("ANA(B1|_B1)", toupper(file)) | grepl(".B1", toupper(file))) ~ "B1",
                         (grepl("ANA(B2|_B2)", toupper(file)) | grepl(".B2", toupper(file))) ~ "B2",
                         grepl("ANA(C1|_C1)", toupper(file)) ~ "C1",
                         (grepl("ANA(C2|_C2)", toupper(file)) | grepl(".C2", toupper(file)) | grepl("_C2", toupper(file))) ~ "C2"),
         file = paste0('/cromwell_root/covid19-hg-imputation-sumstats/20201215/sumstats_imputed/20210114/B2_ALL/',name,"_",ana,'_imputed.txt.gz/'),
         effect = "Z",
         effect_type = "z",
         extra_cols = I(list(c("raiss.imputed", "imputationInfo", "N", "N_cases", "N_controls")))) %>% 
         select(-ana, -extra_cols1, -extra_cols2, -extra_cols3, -extra_cols4)

dd <- dlply(df,1,c)
