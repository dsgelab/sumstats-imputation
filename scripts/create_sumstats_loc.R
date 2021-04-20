rm(list = ls())

library(rjson)
library(data.table)
library(dplyr)

# conf files located in "gs://covid19-hg-analysis/20201215/conf/2/json/B2_ALL.json"
setwd('~')

f <- 'Projects/covid19-hgi/META_ANALYSIS/20201215_conf_5_filtered_json_B2_ALL.json'
d <- fromJSON(file = f)

d <- d$meta

df <- as.data.frame(do.call(rbind, lapply(d, unlist)), stringsAsFactors = F)

df$ancestry <- tolower(sub(".*_", "", df$name))
df$ancestry <- gsub('eur', 'nfe', df$ancestry)
df$ancestry <- gsub('sas', 'eas', df$ancestry)
df$ancestry <- gsub('arab', 'afr', df$ancestry)
df$ancestry <- gsub('his', 'amr', df$ancestry)

df$file <- sub('/cromwell_root/', 'gs://', df$file)

print(f)
studies <- gsub('.{4}$', '', df$name)

unique(studies[order(studies)])
length(unique(studies))

df <- df %>%
  mutate(n = as.numeric(n_cases) + as.numeric(n_controls),
         ana = case_when(grepl("ANA(A1|_A1)", toupper(file)) ~ "A1",
                         (grepl("ANA(A2|_A2)", toupper(file)) | grepl(".A2", toupper(file))) ~ "A2",
                         (grepl("ANA(B1|_B1)", toupper(file)) | grepl(".B1", toupper(file))) ~ "B1",
                         (grepl("ANA(B2|_B2)", toupper(file)) | grepl(".B2", toupper(file))) ~ "B2",
                         grepl("ANA(C1|_C1)", toupper(file)) ~ "C1",
                         (grepl("ANA(C2|_C2)", toupper(file)) | grepl(".C2", toupper(file)) | grepl("_C2", toupper(file))) ~ "C2"),
         name = paste(name, ana, sep = "_")) %>%
  select(file, ancestry, n, name)

fwrite(df[2:3], 'Projects/covid19-hgi/sumstats-imputation/conf/sumstats_loc_B2_ALL_test.txt', sep = '\t', col.names = F)
