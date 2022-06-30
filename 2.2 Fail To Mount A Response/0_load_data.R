
# Initial step to load data

library(lubridate)
library(tibble)
library(dplyr)
library(tidyr)



# serloogy datasets, primary care and blood donors
df_serology_pc <- readRDS("/conf/EAVE/GPanalysis/data/serology_primcare_march22.rds")
df_serology_bd <- readRDS("/conf/EAVE/GPanalysis/data/serology_snbts_march22.rds")

nrow(df_serology_bd)
nrow(df_serology_bd %>% distinct(EAVE_LINKNO))
min(df_serology_bd$Sampledate_iso)
max(df_serology_bd$Sampledate_iso)

nrow(df_serology_bd %>% filter(Sampledate_iso >= '2021-02-01' & Sampledate_iso < '2021-02-14'))
nrow(df_serology_bd %>% filter(Sampledate_iso >= '2021-01-01' & Sampledate_iso < '2021-01-07'))
nrow(df_serology_bd %>% filter(Sampledate_iso >= '2021-07-01' & Sampledate_iso < '2021-07-07'))
nrow(df_serology_bd %>% filter(Sampledate_iso >= '2021-07-07' & Sampledate_iso < '2021-07-14'))




# get the latest QCOVID, remove the Age and Sex as we'll use these from the more accurate demographics file
df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") %>%
             select(-Age,-Sex) 


# get the dataframe created by Elliot
# - this is needed because it contains the vaccine time stamp (a variable we were investigating)
#ndays <- 18
#df_vac <- readRDS("/conf/EAVE/GPanalysis/data/temp/sero_vacc_time.rds") %>%
#          mutate(d2d1 = as.numeric(d2_datetime - d1_datetime,units='days'), # calculate the time between the first two vaccines
#                d3d2 = as.numeric(d3_datetime - d2_datetime,units='days'),
#                d4d3 = as.numeric(d4_datetime - d3_datetime,units='days')) 

#df_vac <- df_vac %>%
#          mutate(
#                d2_datetime = if_else(d2d1>ndays,d2_datetime,d3_datetime), # use d3_datetime for d2_datetime if we think d2 is bad
#                d3_datetime = if_else(d2d1<ndays | d3d2<ndays, d4_datetime,d3_datetime),
#                d4_datetime = case_when((d2d1<ndays | d3d2<ndays | d4d3<ndays) ~ NA)
#                ) %>%
#          select(-d2d1,-d3d2,-d4d3) 


#df_vac_clean <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/c19vaccine.rds") %>% as_tibble()
#df_vac_clean

df_vac <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned.rds") %>% as_tibble()
#df_vac_dvprod
#colnames(df_vac_dvprod)



#df_temp <- df_vac_clean %>% head(100) %>% 
#           group_by(EAVE_LINKNO) %>% 
#           arrange(EAVE_LINKNO,occurrence_time) %>% 
#           select(EAVE_LINKNO,occurrence_time,type,stage,batch) %>% 
#           mutate(n=row_number())

#,shielding,ch_resident,underlying_cond,sererely_immuno_supp,vacc_booster)
df_temp <- df_vac %>% distinct(EAVE_LINKNO,vacc_dose_number,.keep_all = TRUE) %>% 
           select(EAVE_LINKNO,vacc_occurence_date,vacc_product_name,vacc_dose_number,shielding,ch_resident,severely_immuno_supp,immuno_supp) %>%
           rename(datetime=vacc_occurence_date,product=vacc_product_name,n=vacc_dose_number) %>% 
           pivot_wider(names_from=n,values_from=c(datetime,product),names_glue="d{n}_{.value}")

df_vac <- df_temp
rm(df_temp)

# find PCR positive tests
df_cdw <- readRDS("/conf/EAVE/GPanalysis/data/CDW_deduped.rds")

df_cdw <- df_cdw %>% 
  filter(result==1) %>%
  select(EAVE_LINKNO,SpecimenDate,result,death28) %>% 
  mutate(CdwDate=SpecimenDate) %>%
  select(-SpecimenDate,-result) %>% arrange(EAVE_LINKNO) %>% as_tibble()

# load the demographics
df_demo <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")


build_analysis_frame <- function(df_serology, recalculate_qual_result = FALSE){
  
  #df_serology <- df_serology_pc
  #recalculate_qual_result <- TRUE
  
  # create breaks to be used for the vaccine time factor
  breaks <- hour(hm("00:00", "6:00", "12:00", "18:00", "23:59"))
  breaks <- c(-1,6,12,18,24)
  # labels for the breaks
  labels <- c("Night", "Morning", "Afternoon", "Evening")
  
  
  #get some cleaner variables for the IgG levels
  df_temp <- df_serology %>%
    mutate(
      IgG = readr::parse_number(test_result_quant),
      QualResult=test_result_qual
    ) 
  
  ### recacluate the IgG/Qual result
  # - this is needed for the PC which used a different assay 
  # --should not greatly affect these results are not many people had a vaccine by this point in time
  if (recalculate_qual_result){
    print ('recalculating qual result')
    df_temp <- df_temp %>%
      mutate(IgG = ifelse(Sampledate_iso < "2021-04-12",2.6*IgG,IgG),
             QualResult = ifelse(IgG<33.8,"Negative", "Positive")) 
  }
  
  #define the outcome
  df_temp <- df_temp %>% mutate(outcome=ifelse(QualResult=='Negative',1,0))
  
  print ('initial')
  print(df_temp %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
  #start off by joining the serology dataframe with the CDW (PCR tests)
  #calculate a days since injection
  df_temp <- df_temp %>% 
    left_join(df_cdw) %>%
    mutate(days_since_infection = as.numeric(as.Date(Sampledate_iso) - as.Date(CdwDate),units='days')) %>%
    mutate(days_since_infection = ifelse(days_since_infection>=0,days_since_infection,NA)) %>%
    select(-CdwDate) 
    
  print (paste0('+cdw= ',nrow(df_temp)))
  
  #merge now with the vaccine status
  df_temp <- df_temp %>%  left_join(df_vac)
  
  print (paste0('+vac= ',nrow(df_temp)))
 
  #start putting things together
  ndays_to_mount <- 15
  

  df_temp <- df_temp %>%
      mutate(
           days_since_vac1= as.numeric(as.Date(Sampledate_iso) - d1_datetime,units='days'), #calculate days since 1st/2nd/etc. vaccines
           days_since_vac2= as.numeric(as.Date(Sampledate_iso) - d2_datetime,units='days'),
           days_since_vac3= as.numeric(as.Date(Sampledate_iso) - d3_datetime,units='days'),
           days_since_vac4= as.numeric(as.Date(Sampledate_iso) - d4_datetime,units='days'),
           days_between_d1_d2 = ifelse(days_since_vac2>ndays_to_mount,days_since_vac1 - days_since_vac2,'N/A') # temp variable, not used
           )  %>% 
      filter(days_since_vac1 > ndays_to_mount) #initial filter down on samples taken after the first vaccine

  df_temp <- df_temp %>% mutate( #calculate stages: stage{1,2,3,4} if the sample has been taken more than 15 days since dose{1,2,3,4}
           stage = ifelse(days_since_vac4>0 & !is.na(days_since_vac4)& days_since_vac4 > ndays_to_mount,4, 
                        ifelse(days_since_vac3>0 & !is.na(days_since_vac3 & days_since_vac3 > ndays_to_mount),3,
                               ifelse(days_since_vac2>0 & !is.na(days_since_vac2) & days_since_vac2 > ndays_to_mount,2,1)))
    ) %>%
    mutate( #using the new stage variable, get last vaccine and last product 
      days_since_last_vac = get(paste0("days_since_vac",stage)), 
      last_vac_product = get(paste0("d",stage,"_product"))
    ) 
  
  
  print (paste0('+sample after 1st vacc= ',nrow(df_temp)))
  
  #df_temp %>% select(EAVE_LINKNO,QualResult,stage,days_since_vac1,days_since_vac2,days_since_last_vac) %>% group_by(EAVE_LINKNO,stage)
  df_temp <- df_temp %>% 
                group_by(EAVE_LINKNO,stage) %>% 
                arrange(EAVE_LINKNO,stage,desc(QualResult)) %>%
                filter(row_number() ==1) %>% 
                ungroup() %>% 
                mutate(
                  time_of_last_vac= cut(hour(get(paste0("d",stage,"_datetime"))),breaks=breaks,labels=labels)
                ) 
  print (paste0('making single person events 15 days after v1= ',nrow(df_temp)))
  print(df_temp %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
  
  #finally join now with QCOVID 
  df_temp <- df_temp %>%
               left_join(df_qcovid) %>%
               select(-age,-sex) 
  print (paste0('+qcovid= ',nrow(df_temp)))
  print(df_temp %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
    
  #join with demographics
  df_temp <- df_temp %>% 
               left_join(df_demo %>% 
               select(EAVE_LINKNO,Sex,ageYear,DataZone,simd2020v2_sc_quintile))
  print (paste0('+demo= ',nrow(df_temp)))
  return (df_temp)
}


df_temp <- build_analysis_frame(df_serology_pc,recalculate_qual_result = TRUE)
colnames(df_temp)

#df_temp %>% select(d1_datetime,d2_datetime,days_between_d1_d2)
#hist(as.numeric(format(strptime(df_temp$d1_time,"%H:%M:%S"),format="%H")),breaks=24)


df_ana_pc <- df_temp #%>% select(EAVE_LINKNO,Sampledate_iso,Sex,ageYear,last_vac_product,IgG,QualResult,stage,
                      #       days_since_last_vac,days_between_d1_d2,time_of_last_vac,days_since_infection,n_risk_gps,
                      #       contains("Q_"),contains("Cdw"),simd2020v2_sc_quintile,d1_product,d2_product,d3_product)

saveRDS(df_ana_pc,"/conf/EAVE/GPanalysis/data/temp/fail_to_mount_pc.rds")
rm(df_ana_pc)

colnames(df_serology_bd)
df_temp <- build_analysis_frame(df_serology_bd)
colnames(df_temp)

df_ana_bd <- df_temp %>% select(EAVE_LINKNO,outcome,Sampledate_iso,Sex,ageYear,last_vac_product,IgG,QualResult,stage,
                                days_since_last_vac,time_of_last_vac,days_since_infection,n_risk_gps,
                                contains("Q_"),contains("Cdw"),death28,DataZone,simd2020v2_sc_quintile,d1_product,d2_product,d3_product)

saveRDS(df_ana_bd,"/conf/EAVE/GPanalysis/data/temp/fail_to_mount_bd.rds")

rm(df_ana_bd)



colnames(df_ana)


rm(df_temp)
rm(df_qcovid)
rm(df_serology_pc)
rm(df_serology_bd)
rm(df_vac)
rm(df_cdw)
rm(df_demo)
rm(temp)

get_analysis_df <- function(df){
  print (nrow(df))
  print(df %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
  df <- df %>% filter(stage>0 & !is.na(Sex))
  print (nrow(df))
  print(df %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
  df <- df %>% filter(!is.na(simd2020v2_sc_quintile))
  print (nrow(df))
  print(df %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
    
  df <- df %>%  mutate(
    #outcome=ifelse(QualResult=='Negative',1,0),
                  days_since_first_measurement = as.numeric(Sampledate_iso - min(df_ana$Sampledate_iso) ,units='days'),
                  prior_infection = ifelse(is.na(days_since_infection),0,1),
                  #product_time = as.integer(factor(time_of_last_vac,labels=c(1,2,3,4))),
                  product = last_vac_product) %>%
           mutate(product_binary =ifelse(product=="Covid-19 Vaccine AstraZeneca",1,0)) %>% 
           mutate(across(starts_with("Q_"),.fns=as.numeric)) %>%
           mutate(Q_DIAG_CKD_LEVEL=ifelse(Q_DIAG_CKD_LEVEL==0,0,1)) %>% filter(!is.na(Q_BMI))
  print (nrow(df))
  print(df %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  df <- df %>% filter(Q_BMI<80 & Q_BMI>0) 
  print(df %>% group_by(outcome) %>% summarise(n=n()) %>% mutate(per=100*n/sum(n)))
  
  print (nrow(df))
  return (df);
}

df_ana <- readRDS("/conf/EAVE/GPanalysis/data/temp/fail_to_mount_pc.rds")
nrow(df_ana)
df_ana <- df_ana %>% get_analysis_df
df_ana

df_ana$Sampledate_iso


df_ana <- readRDS("/conf/EAVE/GPanalysis/data/temp/fail_to_mount_bd.rds")
nrow(df_serology_bd)
nrow(df_serology_bd %>% distinct(EAVE_LINKNO))
nrow(df_ana)
df_ana <- df_ana %>% get_analysis_df
df_ana
#df_ana <- df_ana %>% get_analysis_df


View(df_ana %>% filter(stage>2) %>% filter(d1_product == 'Covid-19 Vaccine AstraZeneca' & d1_product == d2_product & d1_product == d3_product) %>%
       filter(QualResult=='Negative'))

df_ana %>% group_by(outcome) %>% summarise(n=n())
nrow(df_ana)

df_temp <- df_ana %>%
           mutate(
             y = outcome,
             x1 = ageYear,
             x2 = Q_BMI,
             x3 = as.integer(n_risk_gps) -1,
             x4 = stage,
             x5 = Q_DIAG_BLOOD_CANCER,
             x6 = days_since_last_vac
           ) %>%
           #mutate(
          #   x1 = (x1 - mean(x1)) / mean(x1),
          #   x2 = (x2 - mean(x2)) / mean(x2)
           #) %>%
           select(y,x1,x2,x3,x4,x5,x6)
df_temp 

write.csv(df_temp,'temp.csv')

df_ana %>% group_by(product,product_binary) %>% summarise(n=n())

colnames(df_ana_2)

hist(df_ana_2$Q_BMI)


nrow(df_ana_2)

colnames(df_ana_2)

nrow(df_ana %>% filter(stage>0))
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)))
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>% filter(!is.na(Q_BMI)) )
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>% filter(!is.na(Q_BMI) & Q_BMI>0 & Q_BMI<80) )
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>%filter(!is.na(Q_BMI) & Q_BMI>0 & Q_BMI<80 & !is.na(simd2020v2_sc_quintile)))

df_ana_2 %>% group_by(stage) %>% summarise(n=n())

lookup <- c(
  "days_since_vac2"="Days Since 2nd Vaccination", 
  "n_risk_gps1"="1 Risk",
  "n_risk_gps2"="2 Risk",    
  "n_risk_gps3"="3-4 Risk",     
  "n_risk_gps4"="5+ Risk", 
  "n_risk_gps"="Number of Risk Groups",
  "days_since_last_vac"="Days Since Last Vaccination",
  "days_since_last_vac200+"="200+ Days Since Last Dose",
  "days_since_last_vac150-200"="150-200 Days Since Last Dose",
  "days_since_last_vac100-150"="100-150 Days Since Last Dose",
  "days_since_last_vac50-100"="50-100 Days Since Last Dose",
  "days_since_last_vac30-50"="30-50 Days Since Last Dose",
  "days_since_last_vac15-30"="15-30 Days Since Last Dose",
  "days_since_infection"="Days Since Last Recorded Infection",
  "days_between_d1_d20-25"="18-25 days between Doses",
  "days_between_d1_d225-50"="25-50 days between Doses",
  "days_between_d1_d275+"="More than 75 days between Doses",
  #"days_between_d1_d220-25"="20-25 days between Doses",
  "time_of_last_vac"="Time of Day of Last Vaccination",
  "age"="Age",
  "SexM"="Sex (Male)",
  "ageYear"="Age",
  "ageYear0-20"="Age (0-20)",
  "ageYear20-30"="Age (20-30)",
  "ageYear30-50"="Age (30-50)",
  "ageYear70+"="Age (70+)",
  "ageYear0-19"="Age (0-19)",
  "ageYear20-39"="Age (20-39)",
  "ageYear40-59"="Age (40-59)",
  "ageYear60+"="Age (60+)",
  "age_group"="Age",
  "age_group20-40"="Age 20-40",
  "age_group40-60"="Age 40-60",
  "age_group60-80"="Age 60-80",
  "age_group80+"="Age 80+",
  "sex2"="Sex (Female)",
  "product"="Vaccine (AstraZeneca)",
  "product_binary"="Last Vaccine (AstraZeneca)",
  "product_binary1"="Last Vaccine (AstraZeneca)",
  "product_time"="Time of Day of Last Vaccine",
  "product_time1"="Last Vaccine (Given at Night)",
  "product_time3"="Last Vaccine (Given in the Afternoon)",
  "product_time4"="Last Vaccine (Given in the Evening)",
  "last_vac_product"="Vaccine",
  "Q_BMI"="BMI",
  "Q_BMI0-20"="BMI (0-20) Underweight",
  "Q_BMI30+"="BMI (30+) Obese",
  "Q_BMI25-30"="BMI (25-30) Overweight",
  "Q_DIAG_DIABETES_1"="Diabetes I&II ", 
  "Q_DIAG_DIABETES_2"="Diabetes (Other)",
  "Q_DIAG_SICKLE_CELL.tex"="Sickle Cell Disease",
  "Q_DIAG_RESP_CANCER"="Respitory Cancer",
  "Q_DIAG_ASTHMA"="Asthma",    
  "prior_infection"="Prior SARs-CoV-2 Infection",
  "Q_DIAG_BLOOD_CANCER"="Blood Cancer",
  "Q_DIAG_CHD"="Coronary Heart Disease ",
  "Q_DIAG_COPD"="Chronic Obstructive Pulmonary Disease",   
  "Q_DIAG_CKD_LEVEL"="Chronic Kidney Disease",
  "Q_DIAG_AF"="Atrial Fibrillation",
  "Q_DIAG_CCF"="Heart Failure",
  "Q_DIAG_EPILEPSY"="Epilepsy",     
  "Q_DIAG_FRACTURE"="A prior fracture of hip, wrist, spine or humerus",     
  "Q_DIAG_IMMU"="Immune Deficiency",         
  "Q_DIAG_NEURO"="Neurone Disease",       
  "Q_DIAG_PARKINSONS"="Parkinsons",  
  "Q_DIAG_PULM_RARE"="Cystic Fibrosis or Bronchiectasis or Alveolitis",   
  "Q_DIAG_PVD"="Peripheral Vascular Disease",         
  "Q_DIAG_RA_SLE"="Rheumatoid Arthritis",      
  "Q_DIAG_SEV_MENT_ILL"="Severe Mental Health Illness",
  "Q_DIAG_STROKE"="Stroke",      
  "Q_DIAG_SICKLE_CELL"="Sickle Cell Disease",
  "Q_DIAG_VTE"="Thrombosis or Pulmonary Embolus",
  "Q_DIAG_CEREBRALPALSY"="Cerebralpalsy",
  "Q_DIAG_CIRRHOSIS"="Cirrhosis",
  "Q_DIAG_CONGEN_HD"="Congenital Heart Disease",
  "Q_DIAG_HIV_AIDS"="HIV/Aids",
  "Q_DIAG_PULM_HYPER"="Pulmonary Hypertension",
  "Q_DIAG_DEMENTIA"="Dementia",
  "simd2020v2_sc_quintile"="Scottish Index of Multiple Deprivation",
  "simd2020v2_sc_quintile1"="SIMD 1",
  "simd2020v2_sc_quintile2"="SIMD 2",
  "simd2020v2_sc_quintile3"="SIMD 3",
  "simd2020v2_sc_quintile4"="SIMD 4",
  "simd2020v2_sc_quintile5"="SIMD 5",
  "stage2"="2 Doses",
  "stage1"="1 Dose",
  "stage2+"="More than 2 Doses",
  "immuno_supp"="Immuno Suppressed (Vaccine Record)",
  "shielding"="Shielding",
  "ch_resident"="Care Home Resident",
  "catPf-Pf-Pf"='Mixed 3 Doses (no AZ)',
  "catPf-Pf"='2 Doses Pfizer or Moderna',
  "catmixed 3 doses"='Mixed 3 Doses (including AZ)',
  "catAZ-AZ"='2 Doses AZ',
  "catAZ"='1 Dose AZ',
  "catPf"='1 Dose Pfizer or Moderna'
)

folder <- "/home/calumm09/SerologyAnalysis/2.2 Fail To Mount A Response/tables/all/"
results_folder <- "/home/calumm09/SerologyAnalysis/2.2 Fail To Mount A Response/results/all/"

