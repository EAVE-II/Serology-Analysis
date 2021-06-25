
library(tidyverse)


df_serology <- readRDS("/conf/EAVE/GPanalysis/data/Serology_all.rds")
df_demographics <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
df_comorbid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/Qcovid_update_Jan21.rds")

nrow(df_serology)
nrow(df_demographics)
nrow(df_comorbid)

colnames(df_comorbid)

#head(df_comorbid[order(df_comorbid$EAVE_LINKNO), ],10)

head(df_demographics,10)
head(df_serology,10)
head(df_comorbid,10)

df_merged <- left_join(df_serology,df_demographics)
df_merged

nrow(df_merged)

colnames(df_serology)
colnames(df_merged)

ids <- deframe(df_serology[1:10,'EAVE_LINKNO'])
ids

df_comorbid_small <- df_comorbid %>%
  select(EAVE_LINKNO,cluster) %>%
  filter(EAVE_LINKNO %in% ids)
df_comorbid_small

t <- df_demographics %>% filter(EAVE_LINKNO %in% ids)
t

