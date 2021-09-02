
library(dplyr)
library(tibble)
library("tidyr")

# * open the comorbidity file
# * remove duplicates of comorbidity measurements
#   by ordering on EventDate and taking the most recent
# * keep only the EAVE_LINKNO and comorbidity type (cluster) and value
df_comorbid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/Qcovid_update_Jan21.rds") %>%
              arrange(desc(EventDate)) %>% 
              distinct(EAVE_LINKNO,cluster,.keep_all=TRUE) %>%
              select(EAVE_LINKNO,cluster,Value) %>%
              arrange(EAVE_LINKNO)

nrow(df_comorbid)

# Note: 
# - can we handle comorbidity measurement duplications better?
# - some comorbidities are longer term than others,
#   can we assume someone with a condition in 2015 still has it when Serology test 
#   has been taken?

# Make a better dataframe so there is one row per EAVE_LINKNO
# seperate columns now for different comorbidities
# NA means the person is not classified with said comorbidity 
# None NA means the person is classified (value is comorbidity value, e.g. BMI)
df_comorbid_flat <- df_comorbid %>% 
  group_by(EAVE_LINKNO) %>%
  spread(cluster, Value)

nrow(df_comorbid_flat)

# load the EAVE-II demographics and pull out just the Sex/Age of all studies
df_demographics <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds") %>%
  as_tibble() %>%
  select(EAVE_LINKNO,Sex,ageYear) 

# join the flat comorbidities with these demographics
# if a study is not in the comorbidities, then all clusters with be NA 
df_comorbid_flat <- df_demographics %>%
  left_join(df_comorbid_flat) %>%
  arrange(EAVE_LINKNO)

nrow(df_comorbid_flat)

# save this new dataframe
saveRDS(df_comorbid_flat,"/home/calumm09/data/comorbidities.rds")
