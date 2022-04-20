
# Initial step to load data

df_serology_pc <- readRDS("/conf/EAVE/GPanalysis/data/serology_primcare_march22.rds")
colnames(df_serology_pc)

df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") %>%
             select(-Age,-Sex) 

df_vac <- readRDS("/conf/EAVE/GPanalysis/data/temp/sero_vacc_time.rds")


df_ana <- df_serology_pc %>% 
  left_join(df_vac %>% filter(serology_source=='primary_care_serology')) %>%
  mutate(
         IgG = readr::parse_number(test_result_quant),
         SpecimenDate=as.Date(SpecimenDate),
         QualResult=test_result_qual
        ) %>% 
  filter(as.Date(Sampledate_iso) >= d1_date)  %>% 
  mutate(
         days_since_vac1= as.numeric(as.Date(Sampledate_iso) - d1_date,units='days'),
         days_since_vac2= as.numeric(as.Date(Sampledate_iso) - d2_date,units='days'),
         days_since_vac3= as.numeric(as.Date(Sampledate_iso) - d3_date,units='days'),
         days_since_vac4= as.numeric(as.Date(Sampledate_iso) - d4_date,units='days')
         ) %>%
  mutate(
         stage = ifelse(days_since_vac4>0 & !is.na(days_since_vac4),4,
                      ifelse(days_since_vac3>0 & !is.na(days_since_vac3),3,
                             ifelse(days_since_vac2>0 & !is.na(days_since_vac2),2,1)))
  ) %>%
  mutate(
    days_since_last_vac = get(paste0("days_since_vac",stage))
  ) %>% 
  group_by(EAVE_LINKNO,stage) %>% 
  filter(days_since_last_vac>15) %>% 
  arrange(EAVE_LINKNO,stage,desc(QualResult)) %>%
  filter(row_number() ==1) %>% 
  ungroup() %>% 
  left_join(df_qcovid)

rm(df_qcovid)
rm(df_serology_pc)
rm(df_vac)

df_ana