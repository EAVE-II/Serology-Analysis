library(tibble)
library(dplyr)
library(tidyr)
library(data.table) #!may not be required

# - variable manipulation 
library(forcats) #factors 
library(lubridate) #dates


clean_serology <- function (df){
  df <- df %>%
    mutate(
      IgG = readr::parse_number(test_result_quant),
      QualResult=test_result_qual
    )
  return (df);
}

# A function to load up and clean the two latest and recommended serology datasets
#' @return A list of dataframes corresponding to the datasets we have
#' @export
load_serology <- function() {
  
  #due to a slighly different Diasorin assay being used before the 12th of April 2021 for the primary care
  #we must modify the IgG level by a factor of 2.6 to return consistent units
  #also redefine the definition of a negative/positive by the assay cutoff
  #... more details can be found in the serology papers... 
  df_serology_pc <- readRDS("/conf/EAVE/GPanalysis/data/serology_primcare_july22_v3.rds") %>% 
    clean_serology() %>%
    mutate(IgG = ifelse(Sampledate_iso < "2021-04-12",2.6*IgG,IgG),
           QualResult = ifelse(IgG<33.8,"Negative", "Positive")) %>%
    select(EAVE_LINKNO,Sampledate_iso,QualResult,IgG)
  
  #no rescaling the the blood donors needed
  #Note: IgG not in the same units as for the serology
  df_serology_bd <- readRDS("/conf/EAVE/GPanalysis/data/serology_snbts_july22_v3.rds") %>%
    clean_serology() %>%
    select(EAVE_LINKNO,Sampledate_iso,QualResult,IgG)
  
  return (list(primary_care=df_serology_pc,blood_donors=df_serology_bd));
}

# A function to load up and clean the main eave-ii datasets
#' @return A list of dataframes corresponding to the datasets we have
#' @export
load_recommended_datasets <- function() {
  
  #find all deaths
  df_deaths <- readRDS("/conf/EAVE/GPanalysis/data/all_deaths.rds") %>% 
               select(EAVE_LINKNO,UNDERLYING_CAUSE_OF_DEATH,NRS.Date.Death)
  message(paste0('loaded ',nrow(df_deaths),' deaths'))
  
  #get hospitalisations
  df_smr01 <- readRDS("/conf/EAVE/GPanalysis/data/smr01_2022_06_06.rds")
  #• U07.1 COVID-19, virus identified
  #• U07.2 COVID-19, virus not identified
  #• U07.3 Personal history of COVID-19
  #• U07.4 Post COVID-19 condition
  #• U07.5 Multisystem inflammatory syndrome associated with COVID-19
  #• U07.6 Need for immunization against COVID-19
  #• U07.7 COVID-19 vaccines causing adverse effects in therapeutic use
  
  #select only those hospitalisations where COVID-19 was the main cause for hospitalisation 
  df_smr01_covid <- df_smr01 %>% select(EAVE_LINKNO,ADMISSION_DATE,MAIN_CONDITION) %>% 
    mutate(ADMISSION_COVID = ifelse(MAIN_CONDITION=='U071' | MAIN_CONDITION=='U072',1,0)
  )
  #only keep the smr01 data of covid-19 hospitalisations
  rm(df_smr01)
  message(paste0('loaded ',nrow(df_smr01_covid),' covid-19 hospitalisations'))
  
  
  #load the QCOVID2-3 risk groups and BMI
  df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") %>%
    select(-Age,-Sex) 
  
  message(paste0('loaded ',nrow(df_qcovid),' records of risks from QCOVID-2/3'))
  
  #get the vaccination records 
  #mutate the dataframe so that we have one row per person 
  df_vac <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned_incl_cohorts_20220926.rds") %>% 
    as_tibble() %>% 
    distinct(EAVE_LINKNO,vacc_dose_number,.keep_all = TRUE) %>% 
    select(EAVE_LINKNO,vacc_occurence_date,vacc_product_name,vacc_dose_number,shielding,ch_resident,severely_immuno_supp,immuno_supp) %>%
    rename(datetime=vacc_occurence_date,product=vacc_product_name,n=vacc_dose_number) %>% 
    pivot_wider(names_from=n,values_from=c(datetime,product),names_glue="d{n}_{.value}")
  
  #if there's been 2 doses, make sure the record is good by requiring the products to have been the same
  #if there products are different, something is not right, so throw away the record
  #if d2_product is not NA, require d1_product == d2_product
  df_vac <- df_vac %>%  filter(if(!is.na(d2_product)) d1_product==d2_product else TRUE )
  
  
  message(paste0('loaded ',nrow(df_vac),' vaccination records'))
  
  # find PCR positive tests
  df_pcr_ve <- readRDS("/conf/EAVE/GPanalysis/data/CDW_deduped.rds") %>% 
    filter(result==1) %>%
    select(EAVE_LINKNO,SpecimenDate,result,death28) %>% 
    mutate(CdwDate=SpecimenDate) %>%
    select(-SpecimenDate,-result) %>% arrange(EAVE_LINKNO) %>% as_tibble()
  
  message(paste0('loaded ',nrow(df_pcr_ve),' positive pcr test records'))
  
  # load the demographics
  df_demo <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
  #df_demo <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds")
  
  message(paste0('loaded ',nrow(df_demo),' individual demographics'))
  
  #load the serology data
  df_serology = load_serology()
  
  return(list(
    demographics=df_demo,
    serology_primary_care=df_serology$primary_care,
    serology_blood_donors=df_serology$blood_donors,
    deaths=df_deaths,
    smr01_covid=df_smr01_covid,
    qcovid=df_qcovid,
    c19vaccine=df_vac,
    pcr_ve=df_pcr_ve))
}



# A 
#' @return 
#' @export
serology_vaccine_analysis.create_dataframe <- function(eave.data,serology_dataset_name,ndays_to_mount=14){
  
  #if (serology_dataset_name)
  
  #if eave.data is null...
  
  #log some meta data for making dataflow charts
  meta <- list()
  
  # start the analysis dataframe define our standard outcome
  # Insufficient response is defined by negative test results 
  # NA values of QualResult are likely from 'Equivocal' definition
  #   - therefore we should define these as not an insufficient response
  df <- eave.data[[serology_dataset_name]] %>% 
        mutate(insufficient_response=ifelse((QualResult=='Negative') & (is.na(QualResult) == FALSE),1,0)) 
  
  meta[['Initial']] <- nrow(df)
  
  #select the demographics, retrieving sex,age and SIMD and join with this
  #simd2020_sc_quintile
  #
  df <- df %>% left_join(eave.data$demographics %>% 
                           select(EAVE_LINKNO,Sex,ageYear,simd2020v2_sc_quintile)) %>%
               filter(!is.na(Sex) & !is.na(ageYear))
  
  meta[['Valid Demographics']] <- nrow(df)#length(unique(df$EAVE_LINKNO))

  df <- df %>% left_join(eave.data$qcovid) %>% filter(!is.na(n_risk_gps))
  
  meta[['Valid QCOVID']] <-  nrow(df)
   
  df <- df %>% left_join(eave.data$pcr_ve) %>%
    mutate(days_since_infection = as.numeric(as.Date(Sampledate_iso) - as.Date(CdwDate),units='days')) %>%
    mutate(days_since_infection = ifelse(days_since_infection>=0,days_since_infection,NA)) %>%
    select(-CdwDate) 
    
  #merge now with the vaccine status
  # - filtering on na immuno_supp varriable, but any other variable could be used
  df <- df %>%  left_join(eave.data$c19vaccine) %>% 
    filter(!is.na(immuno_supp))
  
  meta[['Valid Vaccine Record']] <-  nrow(df)#length(unique(df$EAVE_LINKNO))
  
  #calculate days since 1st/2nd/etc. vaccines
  #only concern ourselves with up to 4 vaccines (suitable for this data )
  # - this caculation could be more dynamic, but it'll do...
  df <- df %>%
    mutate(
      days_since_vac1= as.numeric(as.Date(Sampledate_iso) - d1_datetime,units='days'),
      days_since_vac2= as.numeric(as.Date(Sampledate_iso) - d2_datetime,units='days'),
      days_since_vac3= as.numeric(as.Date(Sampledate_iso) - d3_datetime,units='days'),
      days_since_vac4= as.numeric(as.Date(Sampledate_iso) - d4_datetime,units='days')
    )  %>% 
    filter(days_since_vac1 > ndays_to_mount) #initially filter down on samples taken after the first vaccine
  
  
  meta[['Measured After Vaccination']] <- nrow(df)
  
  
  # Calculate the vaccine stage and days since last vaccination & last vaccine product
  # - bit hacky because we have to define a stage based upon the days_since_vacX variable
  # - this is because we need to allow for a number of days (nominally 14) to be a sufficient
  #   time for the host to mount a response to the vaccine
  #   e.g. if someone had a serology test a day after a vaccine,
  #        it's unfair to say they may have had an insufficient response to the vaccine.. 
  #
  # Note: since we used a GAM and a spline of the days_since_last_vac for vaccine waning in the final 
  #       result, this could have absorbed confounding from a lack of time to mount a response..
  #
  # Comment: we performed a sensitivity analysis changing the ndays_to_mount and the final results don't change
  df <- df %>% mutate( 
    stage = case_when(
      days_since_vac4>0 & !is.na(days_since_vac4) & days_since_vac4 > ndays_to_mount ~ 4, 
      days_since_vac3>0 & !is.na(days_since_vac3) & days_since_vac3 > ndays_to_mount ~ 3,
      days_since_vac2>0 & !is.na(days_since_vac2) & days_since_vac2 > ndays_to_mount ~ 2,
      TRUE ~ 1 #always at least stage 1, as we have already filtered on days_since_vac1 > ndays_to_mount
    )) %>%
    mutate(days_since_last_vac = case_when(
                                            stage==1 ~ days_since_vac1,
                                            stage==2 ~ days_since_vac2,
                                            stage==3 ~ days_since_vac3,
                                            TRUE ~ days_since_vac4),
           last_vac_product =  case_when(
                                      stage==1 ~ d1_product,
                                      stage==2 ~ d2_product,
                                      stage==3 ~ d3_product,
                                      TRUE ~ d4_product)
    )
  
  # Group by all measurements in a stage
  # arrange them by the insufficient response result, to see if any insufficient responses were recorded
  # - if an insufficient response was recorded, put this one first in the group
  # - therefore just one measurement per vaccine stage to help define a if an insufficient response registered
  # Note: this only applies to ~100 samples out of >20,000 where there are multiple measurements per stage 
  df <- df %>% 
    group_by(EAVE_LINKNO,stage) %>%
    arrange(EAVE_LINKNO,stage,desc(insufficient_response)) %>%
    filter(row_number()==1) %>% 
    ungroup()
  
  #get rid of some junk variable we end up with from the serology data
  #df <- df %>% select(-LabSpecimenNo,-year,-isoweek,-healthboard,-SpecimenOrigin,-DateCollected)
    
  meta[['After Rearranging Vaccines']] <-  nrow(df)
  
  #do some further quick calculations for convience variables to be used in the analysis 
  df <- df %>%  mutate(
                        days_since_first_measurement = as.numeric(Sampledate_iso - min(df$Sampledate_iso) ,units='days'),
                        prior_infection = ifelse(is.na(days_since_infection),0,1)
                        ) %>% 
                mutate(across(starts_with("Q_"),.fns=as.numeric)) %>%
                mutate(Q_DIAG_CKD_LEVEL=ifelse(Q_DIAG_CKD_LEVEL==0,0,1)) #convert CKD level 3-5 to a binary of having CKD
    
    
  #this shouldnt be needed, but make sure the demographics are valid
  df <- df %>% filter(!is.na(Sex))
  meta[['Valid Sex']] <- nrow(df)
    
      
  #finally fix up the SIMD variable
  df <- df %>% mutate(simd2020v2_sc_quintile=as.character(simd2020v2_sc_quintile)) %>%
               mutate(simd2020v2_sc_quintile=replace_na(simd2020v2_sc_quintile,'Unknown'))
     
  attr(df,'meta') <- meta
  return (df)
}
