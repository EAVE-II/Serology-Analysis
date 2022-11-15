
# Initial step to load data

library(lubridate)
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)
library(forcats)


clean_serology <- function (df,fix_igg=F){
  df <- df %>%
    mutate(
      IgG = readr::parse_number(test_result_quant),
      QualResult=test_result_qual
    ) 
  
  ### recacluate the IgG/Qual result
  # - this is needed for the PC which used a different assay 
  # --should not greatly affect these results are not many people had a vaccine by this point in time
  if (fix_igg){
      print ('fixing IgG')
      df <- df %>%
        mutate(IgG = ifelse(Sampledate_iso < "2021-04-12",2.6*IgG,IgG),
               QualResult = ifelse(IgG<33.8,"Negative", "Positive")) 
  }
  return (df);
}


# serloogy datasets, primary care and blood donors
df_serology_pc <- readRDS("/conf/EAVE/GPanalysis/data/serology_primcare_july22_v3.rds") %>% 
                  clean_serology(T) 
df_serology_bd <- readRDS("/conf/EAVE/GPanalysis/data/serology_snbts_july22_v3.rds") %>%
                  clean_serology(F)


#temp <- df_serology_pc %>% left_join(df_qcovid %>% select(EAVE_LINKNO,n_risk_gps))
#temp %>% select(LabSpecimenNo,n_risk_gps) %>% saveRDS('/conf/EAVE/GPanalysis/data/temp/serology_primcare_july22_nrisks.rds')

#temp <- df_serology_bd %>% left_join(df_qcovid %>% select(EAVE_LINKNO,n_risk_gps))
#temp %>% select(SampleID,n_risk_gps) %>% saveRDS('/conf/EAVE/GPanalysis/data/temp/serology_snbts_july22_nrisks.rds')


    

df_deaths <- readRDS("/conf/EAVE/GPanalysis/data/all_deaths.rds")
df_deaths <- df_deaths %>% 
             select(EAVE_LINKNO,UNDERLYING_CAUSE_OF_DEATH,NRS.Date.Death)
nrow(df_deaths)
#df_deaths

#df_deaths %>% group_by(UNDERLYING_CAUSE_OF_DEATH) %>% summarise(n=n()) %>% arrange(desc(n)) %>% filter(grepl('U07',UNDERLYING_CAUSE_OF_DEATH))

df_smr01 <- readRDS("/conf/EAVE/GPanalysis/data/smr01_2022_06_06.rds")

#• U07.1 COVID-19, virus identified
#• U07.2 COVID-19, virus not identified
#• U07.3 Personal history of COVID-19
#• U07.4 Post COVID-19 condition
#• U07.5 Multisystem inflammatory syndrome associated with COVID-19
#• U07.6 Need for immunization against COVID-19
#• U07.7 COVID-19 vaccines causing adverse effects in therapeutic use

df_smr01 <- df_smr01%>% select(EAVE_LINKNO,ADMISSION_DATE,MAIN_CONDITION) %>% 
            mutate(ADMISSION_COVID = ifelse(MAIN_CONDITION=='U071' | MAIN_CONDITION=='U072',1,0)#,
                   #ADMISSION_WITH_COVID=if_any(starts_with("OTHER_"),  ~ (. == 'U071' | .== 'U072'))
            )

df_severe <- readRDS("/conf/EAVE/GPanalysis/data/cases_severe_dates.rds")
df_severe <- df_severe %>% as_tibble %>% filter(icu | hdu | dead28) %>% select(EAVE_LINKNO,SPECIMENDATE,covid_ucod,covid_cod,icu,hdu,dead28)


#df_smr01 <- df_smr01 %>% mutate(ADMISSION_COVID = ifelse(ADMISSION_COVID==1 | ADMISSION_WITH_COVID==T,1,0))



#mutate(ADMISSION_WITH_COVID=if_any(starts_with("OTHER_"),  ~ (. == 'U071' | .== 'U072'))) %>% 
#  mutate(ADMISSION_WITH_COVID=as.integer(ADMISSION_WITH_COVID)) %>% 
#  mutate(ADMISSION_WITH_COVID=replace(is.na(), 0))



# get the latest QCOVID, remove the Age and Sex as we'll use these from the more accurate demographics file
df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") %>%
             select(-Age,-Sex) 

#df_qcovid <- df_diabetes %>% select(EAVE_LINKNO,contains('Q_'),n_risk_gps) %>% as_tibble
#df_qcovid
#df_diabetes = readRDS("/conf/EAVE/GPanalysis/progs/CR/Vaccine/output/temp/Qcovid.rds") #%>%
#colnames(df_diabetes)
#              select(EAVE_LINKNO, Q_DIAG_DIABETES_1, Q_DIAG_DIABETES_2) 

#df_qcovid = full_join(df_qcovid, df_diabetes) 

#df_qcovid %>% group_by(Q_DIAG_DIABETES_2) %>% summarise(n=n())

#df_qcovid <- df_qcovid %>% filter(!is.na(n_risk_gps)) %>%
#  mutate(n_risk_gps=rowSums(across(starts_with("Q_DIAG"), ~ .x>0 ))) %>%
#  mutate(n_risk_gps = factor(n_risk_gps,ordered = T,levels=c(0,1,2,3,5),labels=c('0','1','2','3-4','5+')))
#df_qcovid %>% group_by(Q_DIAG_DIABETES_2) %>% summarise(n=n())
#df_qcovid %>% group_by(n_risk_gps) %>% summarise(n=n())

df_vac <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/restored/C19vaccine_dvprod_cleaned_restored.rds") %>% as_tibble()
#df_vac <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/C19vaccine_dvprod_cleaned.rds") %>% as_tibble()
colnames(df_vac)
nrow(df_vac)
df_vac %>% select(shielding,immuno_supp)

#,shielding,ch_resident,underlying_cond,sererely_immuno_supp,vacc_booster)
df_temp <- df_vac %>% distinct(EAVE_LINKNO,vacc_dose_number,.keep_all = TRUE) %>% 
           select(EAVE_LINKNO,vacc_occurence_date,vacc_product_name,vacc_dose_number,shielding,ch_resident,severely_immuno_supp,immuno_supp) %>%
           rename(datetime=vacc_occurence_date,product=vacc_product_name,n=vacc_dose_number) %>% 
           pivot_wider(names_from=n,values_from=c(datetime,product),names_glue="d{n}_{.value}")

df_temp
nrow(df_temp)

df_vac <- df_temp
rm(df_temp)

# find PCR positive tests
df_cdw <- readRDS("/conf/EAVE/GPanalysis/data/CDW_deduped.rds")
colnames(df_cdw)

df_cdw <- df_cdw %>% 
  filter(result==1) %>%
  select(EAVE_LINKNO,SpecimenDate,result,death28) %>% 
  mutate(CdwDate=SpecimenDate) %>%
  select(-SpecimenDate,-result) %>% arrange(EAVE_LINKNO) %>% as_tibble()

# load the demographics
df_demo_orig <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
#df_demo <- readRDS("/conf/EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds")
#df_demo <- df_demo %>% rename(simd2020v2_sc_quintile = simd2020_sc_quintile) %>%  select(EAVE_LINKNO,Sex,ageYear,simd2020v2_sc_quintile) 

df_demo <- df_demo_orig
rm(df_demo_orig)
#df_demo <- df_demo_v3
#rm(df_demo_v2)
#rm(df_demo_v3)




#df_serology <- df_serology_pc 
#recalculate_qual_result <- TRUE

# create breaks to be used for the vaccine time factor
#breaks <- lubridate::hour(lubridate::hm("00:00", "6:00", "12:00", "18:00", "23:59"))
#breaks <- c(-1,6,12,18,24)
# labels for the breaks
#labels <- c("Night", "Morning", "Afternoon", "Evening")


build_analysis_frame <- function(df_serology){
  
  meta <- list()
  
  #define our standard outcome
  df <- df_serology %>% mutate(outcome=ifelse(QualResult=='Negative',1,0))
  
  meta[['Initial']] <- nrow(df)#ength(unique(df$EAVE_LINKNO))

  
  #select the demographics, retrieving sex,age and SIMD and join with this
  #simd2020_sc_quintile
  #
  df <- df %>% left_join(df_demo %>% 
                select(EAVE_LINKNO,Sex,ageYear,simd2020v2_sc_quintile)) %>%
                filter(!is.na(Sex) & !is.na(ageYear))
  
  meta[['Valid Demographics']] <- nrow(df)#length(unique(df$EAVE_LINKNO))
  
  #join with Qcovid and select people who have a valid record
  print (df_qcovid %>% filter(!is.na(n_risk_gps)) %>% group_by(Q_DIAG_DIABETES_2) %>% summarise(n=n()))
  print (df_qcovid %>% group_by(Q_DIAG_DIABETES_2) %>% summarise(n=n()))
  
  df <- df %>%
    left_join(df_qcovid) %>% filter(!is.na(n_risk_gps))
 
  meta[['Valid QCOVID']] <-  nrow(df)#length((df$EAVE_LINKNO))
  

  #join the dataframe with the CDW (PCR tests)
  #calculate a days since injection
  df <- df %>% 
    left_join(df_cdw) %>%
    mutate(days_since_infection = as.numeric(as.Date(Sampledate_iso) - as.Date(CdwDate),units='days')) %>%
    mutate(days_since_infection = ifelse(days_since_infection>=0,days_since_infection,NA)) %>%
    select(-CdwDate) 
  
  #meta[['Valid CdW']] <-  length(unique(df$EAVE_LINKNO))

  
  #merge now with the vaccine status
  df <- df %>%  left_join(df_vac) %>% filter(!is.na(immuno_supp))
  
  meta[['Valid Vaccine Record']] <-  nrow(df)#length(unique(df$EAVE_LINKNO))
  

  ndays_to_mount <- 14
  

  df <- df %>%
      mutate(
           days_since_vac1= as.numeric(as.Date(Sampledate_iso) - d1_datetime,units='days'), #calculate days since 1st/2nd/etc. vaccines
           days_since_vac2= as.numeric(as.Date(Sampledate_iso) - d2_datetime,units='days'),
           days_since_vac3= as.numeric(as.Date(Sampledate_iso) - d3_datetime,units='days'),
           days_since_vac4= as.numeric(as.Date(Sampledate_iso) - d4_datetime,units='days'),
           days_between_d1_d2 = ifelse(days_since_vac2>ndays_to_mount,days_since_vac1 - days_since_vac2,'N/A') # temp variable, not used
           )  %>% 
      filter(days_since_vac1 > ndays_to_mount) #initial filter down on samples taken after the first vaccine

 
  meta[['Measured After Vaccination']] <- nrow(df)# length(unique(df$EAVE_LINKNO))
  
  
  df <- df %>% mutate( #calculate stages: stage{1,2,3,4} if the sample has been taken more than 15 days since dose{1,2,3,4}
           stage = ifelse(days_since_vac4>0 & !is.na(days_since_vac4)& days_since_vac4 > ndays_to_mount,4, 
                        ifelse(days_since_vac3>0 & !is.na(days_since_vac3) & days_since_vac3 > ndays_to_mount,3,
                               ifelse(days_since_vac2>0 & !is.na(days_since_vac2) & days_since_vac2 > ndays_to_mount,2,1)))
    ) 
  
  
  df <- df %>%
    mutate(days_since_last_vac = if_else(stage==1,days_since_vac1,
                                         if_else(stage==2,days_since_vac2,
                                                 if_else(stage==3,days_since_vac3,days_since_vac4))),
           last_vac_product = if_else(stage==1,d1_product,
                                      if_else(stage==2,d2_product,
                                              if_else(stage==3,d3_product,d4_product)))
           )
  
  
  df <- df %>% 
                group_by(EAVE_LINKNO,stage) %>% 
                arrange(EAVE_LINKNO,stage,desc(QualResult)) %>%
                filter(row_number() ==1) %>% 
                ungroup() #%>% 
                #mutate(
                #  time_of_last_vac= cut(hour(get(paste0("d",stage,"_datetime"))),breaks=breaks,labels=labels)
                #) 

  
  meta[['After Rearranging Vaccines']] <-  length(unique(df$EAVE_LINKNO))
  
  
  attr(df,'meta') <- meta
  return (df)

  
  ##join with SMR01
  #min_sample_date <- min(df_temp$Sampledate_iso)
  ##number of days to use to define a hospitalisation outcome
  #ndays <- as.difftime(1500,units='days')
  ##join with SMR01 (filter SMR01 on admissions after the minimum sample date to make it quicker)
  
  #df_temp <- df_temp %>% left_join(df_smr01 %>% filter(ADMISSION_DATE >= min_sample_date) )
  
  ##set the admission date as NA if the admission happened :
  ## * after the sample date + ndays
  ## * OR before the sample date
  
  #df_temp <- df_temp %>%
  #  mutate(ADMISSION_DATE=fifelse(ADMISSION_DATE<(Sampledate_iso+ndays) & ADMISSION_DATE>=Sampledate_iso,
  #                                ADMISSION_DATE,
  #                                as.POSIXct(NA)))
  
  ##modify the ADMISSION variables (set to 0 if now admission date is NA)
  #df_temp <- df_temp %>% 
  #  mutate(ADMISSION_COVID=ifelse(!is.na(ADMISSION_DATE),ADMISSION_COVID,0),
  #         ADMISSION_WITH_COVID=ifelse(!is.na(ADMISSION_DATE),ADMISSION_WITH_COVID,0),
  #         days_hospitalisation = as.numeric(ADMISSION_DATE - Sampledate_iso,units='days')) 
  
  ## * groupby the serology samples
  ## * arrange them so admissions from covid come first
  ##   followed by the admission closest in date to the serology sample
  ## * filter to select the first hospitalisation found
  #df_temp <- df_temp %>% group_by(EAVE_LINKNO,Sampledate_iso) %>% 
  #  arrange(desc(ADMISSION_COVID),abs(days_hospitalisation)) %>% 
  #  filter(row_number()==1) %>% 
  #  ungroup %>% arrange(EAVE_LINKNO) 
  
  # finally calculate some outcome variables 
  # better in get analysis dataframe?
  #df_temp <- df_temp %>% 
  #  mutate(
  #    outcome_h_for_covid = ifelse(!is.na(ADMISSION_COVID),ADMISSION_COVID,0),
  #    outcome_h_with_covid = ifelse(!is.na(ADMISSION_WITH_COVID),ADMISSION_WITH_COVID,0)) %>% 
  #  mutate( 
  #    outcome_h_covid = ifelse(outcome_h_for_covid==0,outcome_h_with_covid,1)
  #    ) %>%
  #  mutate(
  #    outcome_hospitalisation = ifelse(outcome_h_covid==1,0,ifelse(is.na(days_hospitalisation),0,1))
  #  )
  
  #df_temp <- df_temp %>% left_join(df_deaths %>% mutate(death_date = NRS.Date.Death) %>% select(EAVE_LINKNO,death_date))
  #return (df_temp)
}


df_ana_pc <- build_analysis_frame(df_serology_pc)

meta.pc <- attributes(df_ana_pc)$meta
meta.pc


df_ana_bd <- build_analysis_frame(df_serology_bd)

meta.bd <- attributes(df_ana_bd)$meta
meta.bd

ntotal <- meta.pc[[1]] + meta.bd[[1]]

gr<- DiagrammeR::grViz("
digraph graph2  {

graph [layout = dot, overlap = false]

# node definitions with substituted label text
node [shape = rectangle, style=filled, color='@@7'] 
start [color='white', label='Primary Care' ]
a [label = '@@1']
b [label = '@@2']
c [label = '@@3']
d [label = '@@4']
e [label = '@@5']
#f [label = '@@6']


start -> a -> b -> c -> d -> e;


start2 [color='white', label='Blood Donors' ]
a2 [label = '@@8', color='@@13']
b2 [label = '@@9', color='@@13']
c2 [label = '@@10', color='@@13']
d2 [label = '@@11', color='@@13']
e2 [label = '@@12', color='@@13']

start2 -> a2 -> b2 -> c2 -> d2 -> e2;

init[label = '@@14', color='@@15']
init -> start,start2

}

[1]: paste0(labels(meta.pc)[[1]], ': ',meta.pc[[1]])
[2]: paste0(labels(meta.pc)[[2]], ': ',meta.pc[[2]])
[3]: paste0(labels(meta.pc)[[3]], ': ',meta.pc[[3]])
[4]: paste0(labels(meta.pc)[[4]], ': ',meta.pc[[4]])
[5]: paste0(labels(meta.pc)[[5]], ': ',meta.pc[[5]])
[6]: paste0(labels(meta.pc)[[6]], ': ',meta.pc[[6]])
[7]: phsstyles::phs_colours('phs-liberty-10')
[8]: paste0(labels(meta.bd)[[1]], ': ',meta.bd[[1]])
[9]: paste0(labels(meta.bd)[[2]], ': ',meta.bd[[2]])
[10]: paste0(labels(meta.bd)[[3]], ': ',meta.bd[[3]])
[11]: paste0(labels(meta.bd)[[4]], ': ',meta.bd[[4]])
[12]: paste0(labels(meta.bd)[[5]], ': ',meta.bd[[5]])
[13]: phsstyles::phs_colours('phs-liberty-30')
[14]: paste0('NSamples: ',ntotal)
[15]: phsstyles::phs_colours('phs-rust-50')

")

gr

# 2. Convert to SVG, then save as png
tmp = DiagrammeRsvg::export_svg(gr)
tmp = charToRaw(tmp) # flatten
tmp


saveRDS(df_ana_pc,"/conf/EAVE/GPanalysis/data/temp/fail_to_mount_pc_2.rds")
rm(df_ana_pc)


df_ana_bd <- build_analysis_frame(df_serology_bd)
saveRDS(df_ana_bd,"/conf/EAVE/GPanalysis/data/temp/fail_to_mount_bd_2.rds")
rm(df_ana_bd)


rm(df_temp)
rm(df_qcovid)
rm(df_serology_pc)
rm(df_serology_bd)
rm(df_vac)
rm(df_cdw)
rm(df_demo)
rm(temp)

get_summary <- function (df){
  summary <- df %>% group_by(outcome) %>% summarise(nfail=n()) %>% mutate(per=100*nfail/sum(nfail),total=sum(nfail)) %>% filter(outcome==1)
  return (summary %>% select(-outcome));
}

get_analysis_df <- function(df){
  
  df <- df %>%  mutate(
    days_since_first_measurement = as.numeric(Sampledate_iso - min(df$Sampledate_iso) ,units='days'),
    prior_infection = ifelse(is.na(days_since_infection),0,1),
    #product_time = as.integer(factor(time_of_last_vac,labels=c(1,2,3,4))),
    product = last_vac_product) %>%
    mutate(product_binary =ifelse(product=="Covid-19 Vaccine AstraZeneca",1,0)) %>% 
    mutate(across(starts_with("Q_"),.fns=as.numeric)) %>%
    mutate(Q_DIAG_CKD_LEVEL=ifelse(Q_DIAG_CKD_LEVEL==0,0,1))
  
  meta <- list()
  meta[['Total']] <- df %>% get_summary
  
  df <- df %>% filter(!is.na(Sex))
  meta[['Valid Sex']] <- df %>% get_summary
  
  #require a valid QCOVID record
  df <- df %>% filter(!is.na(n_risk_gps))
  meta[['Valid QCOVID']] <- df %>% get_summary

  #allow simd, qcovid and bmi to be NA
  ddf <- df %>% filter(!is.na(simd2020v2_sc_quintile))
  df <- df %>% mutate(simd2020v2_sc_quintile=as.character(simd2020v2_sc_quintile)) %>%
               mutate(simd2020v2_sc_quintile=replace_na(simd2020v2_sc_quintile,'Unknown'))
                
  meta[['Valid SIMD']] <- ddf %>% get_summary
  
  
  ddf <- ddf %>% filter(!is.na(Q_BMI) & Q_BMI<80 & Q_BMI>0) 
  meta[['Valid BMI']] <- ddf %>% get_summary
  
  attr(df,'meta') <- meta
  return (df);
}


load_pc <- function(){
  df_ana <- readRDS("/conf/EAVE/GPanalysis/data/temp/fail_to_mount_pc.rds")
  df_ana <- df_ana %>% get_analysis_df
  return (df_ana)
}

load_bd <- function(){
  df_ana <- readRDS("/conf/EAVE/GPanalysis/data/temp/fail_to_mount_bd.rds")
  df_ana <- df_ana %>% get_analysis_df
  return (df_ana)
}


df_ana_pc <- load_pc()
nrow(df_ana_pc)
data <- attributes(df_ana_pc)$meta
data


temp <- df_ana_bd %>% left_join(df_demo %>% as_tibble %>% select(EAVE_LINKNO,ur6_2016_name) %>% 
                             mutate(ur6_2016_name=as.factor(ifelse(is.na(ur6_2016_name),'Unknown',ur6_2016_name)))) %>%
                       create_modelA %>% 
                       mutate(outcome = ifelse(is.na(outcome),1,outcome)) 
nrow(temp)

temp %>% group_by(outcome) %>% summarise(n=n())
temp %>% ggplot(aes(x=IgG,fill=ur6_2016_name)) + geom_histogram() + scale_x_log10() + scale_y_log10()


tab <- temp %>% group_by(ur6_2016_name) %>% summarise(n=n(),n2=sum(outcome==1)) %>% mutate(ntotal=round(100*n/sum(n),2)) %>%
                           mutate(N=paste0(n,' (',ntotal,')'),NF=paste0(n2," (",round(100*n2/n,2),")")) %>% 
                           select(-n,-ntotal,-n2)
tab

View(tab)

df_ana_bd <- load_bd()
data <- attributes(df_ana_bd)$meta
data


gr <- DiagrammeR::grViz("
digraph graph2  {

graph [layout = dot, overlap = false]

# node definitions with substituted label text
node [shape = rectangle, width = 4, style=filled, color='@@6'] 
start [style = invis ]
a [label = '@@1']
b [label = '@@2']
c [label = '@@3']
d [label = '@@4' color='@@7']
e [label = '@@5' color='@@7']

start -> a [label='15 days after at least one Vaccination '];
a -> b [label='Valid Sex & Age'];
b -> c [label='Valid QCOVID'];
c -> d [label='Valid SIMD'];
d -> e [label='Valid BMI'];


}

[1]: paste0('N(Total) = ', data$Total$total,' ; N(Impaired) = ',data$Total$nfail,' [',round(data$Total$per,2),'%]')
[2]: paste0('N(Total) = ', data[['Valid Sex']]$total,' ; N(Impaired) = ',data[['Valid Sex']]$nfail,' [',round(data[['Valid Sex']]$per,2),'%]')
[3]: paste0('N(Total) = ', data[['Valid QCOVID']]$total,' ; N(Impaired) = ',data[['Valid QCOVID']]$nfail,' [',round(data[['Valid QCOVID']]$per,2),'%]')
[4]: paste0('N(Total) = ', data[['Valid SIMD']]$total,' ; N(Impaired) = ',data[['Valid SIMD']]$nfail,' [',round(data[['Valid SIMD']]$per,2),'%]')
[5]: paste0('N(Total) = ', data[['Valid BMI']]$total,' ; N(Impaired) = ',data[['Valid BMI']]$nfail,' [',round(data[['Valid BMI']]$per,2),'%]')
[6]: phsstyles::phs_colours('phs-liberty-10')
[7]: phsstyles::phs_colours('phs-liberty-30')
")

gr



df_ana <- readRDS("/conf/EAVE/GPanalysis/data/temp/fail_to_mount_bd.rds")
df_ana <- df_ana %>% get_analysis_df
df_ana





#df_ana <- df_ana %>% get_analysis_df



nrow(df_ana_2)

colnames(df_ana_2)

nrow(df_ana %>% filter(stage>0))
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)))
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>% filter(!is.na(Q_BMI)) )
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>% filter(!is.na(Q_BMI) & Q_BMI>0 & Q_BMI<80) )
nrow(df_ana %>%  filter(stage>0 & !is.na(Sex)) %>%filter(!is.na(Q_BMI) & Q_BMI>0 & Q_BMI<80 & !is.na(simd2020v2_sc_quintile)))

df_ana_2 %>% group_by(stage) %>% summarise(n=n())

lookup <- c(
  "ur6_2016_name"="Location",
  "outcome_igg"='IgG Level',
  "additional"="Subsequent Vaccinations",
  "days_since_start"="Pandemic Period",
  "days_since_vac2"="Days Since 2nd Vaccination", 
  "n_risk_gpsOther"="Other Risks",
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
  "Q_BMIUnknown"="Unknown BMI",
  "Q_DIAG_DIABETES_1"="Diabetes I&II ", 
  "Q_DIAG_DIABETES_2"="Diabetes (Other)",
  #"Q_DIAG_DIABETES_1"="Diabetes I", 
  #"Q_DIAG_DIABETES_2"="Diabetes II",
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
  "simd2020v2_sc_quintile1 - High"="SIMD 1 (High)",
  "simd2020v2_sc_quintile1"="SIMD 1",
  "simd2020v2_sc_quintile2"="SIMD 2",
  "simd2020v2_sc_quintile3"="SIMD 3",
  "simd2020v2_sc_quintile4"="SIMD 4",
  "simd2020v2_sc_quintile5"="SIMD 5",
  "simd2020v2_sc_quintile5-Low"="SIMD 5 (Low)",
  "simd2020v2_sc_quintileNA"="Unknown SIMD",
  "simd2020v2_sc_quintileUnknown"="SIMD Unknown",
  "stage2"="2 Doses",
  "stage1"="1 Dose",
  "stage2+"="More than 2 Doses",
  "immuno_supp"="Immuno Suppressed (Vaccine Record)",
  "immuno100"="Immunosuppressed",
  "immuno110"="Severely Immunosuppressed",
  "immuno111"="Immunodeficient",
  "shielding"="Shielding",
  "ch_resident"="Care Home Resident",
  "catPf-Pf-Pf"='Mixed 3-4 Doses (no AZ)',
  "catPf-Pf"='2 Doses Pfizer or Moderna',
  "catmixed 3 doses"='Mixed 3-4 Doses (including AZ)',
  "catmixed 4 doses"='Mixed 4 Doses',
  "catAZ-AZ"='2 Doses AZ',
  "catAZ"='1 Dose AZ',
  "cat2"='Dose(s)',
  "catPf"='1 Dose Pfizer or Moderna',
  "outcome"="Insufficient IgG",
  "outcome_igg75%"="IgG level 75% quantile",
  "outcome_igg50%"="IgG level 50% quantile",
  "outcome_igg25%"="IgG level 25% quantile",
  "outcome_igg10%"="IgG level 10% quantile",
  "outcome_igg5%"="IgG level 5% quantile",
  "days_since_first_measurement"="Days Since 20th December 2020",
  "days_since_last_vac"="Days Since Last Vaccination",
  "cat3 doses"='Mixed 3-4 Doses (no AZ)'
)

folder <- "/home/calumm09/SerologyAnalysis/2.2 Fail To Mount A Response/FinalPlots/"
results_folder <- "/home/calumm09/SerologyAnalysis/2.2 Fail To Mount A Response/results/all/"

