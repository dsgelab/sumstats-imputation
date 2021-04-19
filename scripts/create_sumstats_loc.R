rm(list = ls())

library(rjson)
library(data.table)
library(dplyr)

# conf files located in "gs://covid19-hg-analysis/20201215/conf/2/json/B2_ALL.json"
setwd('~')

f <- 'Projects/covid19-hgi/META_ANALYSIS/20201215_conf_5_filtered_json_C2_ALL.json'
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
  mutate(n = as.numeric(n_cases) + as.numeric(n_controls)) %>%
  select(file, ancestry, n, name)

# fwrite(df, 'Projects/covid19-hgi/META_ANALYSIS/sumstats_loc_C2_ALL.txt', sep = '\t', col.names = F)



 
# df <- df %>%
#   mutate(n = as.numeric(n_cases) + as.numeric(n_controls),
#          n_eff = 4 * as.numeric(n_cases) * as.numeric(n_controls)/n) %>%
#   select(name, n_cases, n_controls, n, n_eff)
# fwrite(df, 'Projects/covid19-hgi/META_ANALYSIS/neff_A2_ALL.tsv', sep = '\t')


#fwrite(df, 'Projects/covid19-hgi/META_ANALYSIS/info_A2_ALL.txt', sep = '\t', col.names = T)


# # # # # # # # ## 
s <- fread('Projects/covid19-hgi/META_ANALYSIS/20201512/imputation_sumstats/plots/B2_ALL/success.txt', header = F)

s$V1 <- sub('.combined_imputed.txt.gz_p.value_manhattan.png','',s$V1)

done <- df %>%
  filter(name %in% s$V1)

to_do <- df %>%
  filter(!name %in% s$V1)

fwrite(done, 'Projects/covid19-hgi/META_ANALYSIS/sumstats_loc_B2_ALL_done.txt', sep = '\t', col.names = F)
fwrite(to_do, 'Projects/covid19-hgi/META_ANALYSIS/sumstats_loc_B2_ALL_to_do.txt', sep = '\t', col.names = F)


# shards that actually output 22 chr
shards_ok <- c(0, 2, 3, 4, 6, 7, 10, 12, 14, 16, 18, 19, 20, 21, 25, 28, 29, 31, 32, 33, 37) + 1

done_ok <- df$name[shards_ok]

df_done_ok <- df %>%
  filter(name %in% done_ok)

df_done_not <- done %>%
  filter(!name %in% done_ok)

fwrite(df_done_ok, 'Projects/covid19-hgi/META_ANALYSIS/sumstats_loc_B2_ALL_done_ok.txt', sep = '\t', col.names = F)


# create conf json for meta

dff <- df %>%
  select(-extra_cols1, -extra_cols2, -extra_cols3, -extra_cols4, -ancestry, -se)

dff$file <- paste0('gs://covid19-hg-imputation-sumstats/20201215/sumstats_imputed/20210114/C2_ALL/', dff$name, '.combined_imputed.txt.gz') 
dff$effect <- 'Z'
dff$effect_type <- 'z'
dff$extra_cols <- ''

dff$file <- sub('gs://', '/cromwell_root/', dff$file)

jsonlite::write_json(dff, 'Projects/covid19-hgi/META_ANALYSIS/meta_conf/C2_ALL.json', pretty = T, auto_unbox = T)


# create sumstats loc for add_N

dff <- df %>%
  filter(ancestry == "nfe" | ancestry == "fin")
dff$file <- paste0('gs://covid19-hg-imputation-sumstats/20201215/sumstats_imputed/B2_ALL/', dff$name, '.combined_imputed.txt.gz')
dff$file

fwrite(dff, 'Projects/covid19-hgi/META_ANALYSIS/sumstats_loc_B2_ALL_eur_imputed.txt', sep = '\t', col.names = F)
