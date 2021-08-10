
library(dplyr)
library("tidyr")
library(tibble)


df_serology <- readRDS("/conf/EAVE/GPanalysis/data/Serology_all.rds")
df_demographics <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
df_comorbid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/Qcovid_update_Jan21.rds")
df_vaccine <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine.rds")


nrow(df_serology)
nrow(df_demographics)
nrow(df_comorbid)
nrow(df_vaccine)

colnames(df_serology)

cols_to_keep = colnames(df_serology)

df_serology_small <- df_serology %>%
                     select(cols_to_keep) %>%
                     mutate(ID = readr::parse_number(as.character(EAVE_LINKNO))) %>%
                     select(-c("EAVE_LINKNO")) %>%
                     select(ID,everything()) %>%
                     arrange(ID)

head(df_serology_small)

ids <- deframe(df_serology_small[,'ID'])

colnames(df_demographics)
df_demographics_small <- df_demographics %>%
  select(c("EAVE_LINKNO","DataZone","simd2020v2_sc_quintile")) %>%
  mutate(ID = readr::parse_number(as.character(EAVE_LINKNO))) %>%
  select(-c("EAVE_LINKNO")) %>%
  select(ID,everything()) %>%
  arrange(ID) %>%
  filter(ID %in% ids)

head(df_demographics_small)
nrow(df_demographics_small)

colnames(df_vaccine)
df_vaccine_small <- df_vaccine %>%
  mutate(ID = readr::parse_number(as.character(EAVE_LINKNO))) %>%
  select(-c("EAVE_LINKNO")) %>%
  select(ID,everything()) %>%
  arrange(ID) %>%
  filter(ID %in% ids)

head(df_vaccine_small)
nrow(df_vaccine_small)

colnames(df_comorbid)
df_comorbid_small <- df_comorbid %>%
  select(c("EAVE_LINKNO","cluster","Value")) %>%
  mutate(ID = readr::parse_number(as.character(EAVE_LINKNO))) %>%
  select(-c("EAVE_LINKNO")) %>%
  select(ID,everything()) %>%
  arrange(ID) %>%
  filter(ID %in% ids)
  
head(df_comorbid_small)
nrow(df_comorbid_small)


write.csv(df_serology_small,"serology.csv")
write.csv(df_demographics_small,"demo.csv")
write.csv(df_comorbid_small,"comorbid.csv")

df_temp <- df_serology %>%
     left_join(y=df_demographics)

head(df_temp)

df_merged <- df_comorbid_small %>%
  left_join(y=df_temp) %>%
  group_by(EAVE_LINKNO)

nrow(df_merged)
tail(df_merged)

df <- df_merged %>% drop_na()
nrow(df)
head(df)

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

df_merged <- df_comorbid_small %>%
  left_join(y=df_demographics) %>%
  left_join(y=df_serology) %>%
  group_by(cluster)
df_merged
