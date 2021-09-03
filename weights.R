
library(dplyr)
library("tidyr")
library(tibble)
library(survey)

# load healthboard population estimates for Scotland for 2019
df_pop_estimate <- readRDS("/conf/EAVE/GPanalysis/data/lookups/HB2019_pop_lookup.rds")
# extract the population distribution for age and sex at the same time
df_pop_age_sex <- df_pop_estimate %>% group_by(age,sex,sex_name) %>% summarise(pop = sum(pop))
df_pop_age_sex

# summarise/get the population age distribution only
pop_age_dist <- df_pop_age_sex %>% group_by(age) %>% summarise(Freq=sum(pop))
pop_age_dist


# load the serology dataset (first 6 month of the pandemic)
df_serology_first_6months <- readRDS("/conf/EAVE/GPanalysis/data/Serology_all.rds")

# load the rest of the serology data
df_serology_Oct20May21 <- readRDS("/conf/EAVE/GPanalysis/data/Serology_Oct20May21.rds")

# print column names as they're different 
colnames(df_serology_first_6months)
colnames(df_serology_Oct20May21)



# - comorbidities from QCOVID + demographics 
#   created by prep_datasets.R
# - need some manipulation to filter out bad ages:
#   some ages are -1, some ages are >130, can't be possible 
#   therefore filter so ages are what are possible in the HB 
#   population estimates
df_comorbid <- readRDS("/home/calumm09/data/comorbidities.rds") %>%
               rename(age = ageYear) %>%
               filter(age %in% pop_age_dist$age )


df_serology_first_6months.unique_people <- df_serology_first_6months %>% 
                                           distinct(EAVE_LINKNO) %>%
                                           add_column(had_serology_test=1L)
df_serology_first_6months.unique_people

df_serology_Oct20May21.unique_people <- df_serology_Oct20May21 %>% 
                                        distinct(EAVE_LINKNO) %>%
                                        add_column(had_serology_test=1L)
df_serology_Oct20May21.unique_people



df_comorbid.serology_first_6months <- df_comorbid %>% 
                                      left_join(df_serology_first_6months.unique_people) %>%
                                      replace(is.na(.),0L) 

df_comorbid.serology_first_6months


calculate_weights <- function (dataset, subset_filter, comorbidity) {
  freq <- count(dataset,!!as.name(comorbidity) ) %>%
          rename(Freq=n)

  subset <- dataset %>% 
             filter(!!as.name(subset_filter) > 0 ) %>%
             select(!!as.name(comorbidity))
  
  freq.subset <- count(subset, !!as.name(comorbidity))  %>% rename(Freq=n)

  print(freq$Freq / sum(freq$Freq))
  print(freq.subset$Freq / sum(freq.subset$Freq))
  
  
  # Use the survey package to create the initial unweighted data
  svy.unweighted <- svydesign(ids=~1, 
                              data=subset,
                              weights=NULL)
  
  svy.rake <- rake(design = svy.unweighted,
                            sample.margins = list(as.formula(paste("~",comorbidity)) ),
                            population.margins = list(freq))
  
  
  # Using the survey package svymean, get the mean/SE weights for Age & Sex when there is no weighting
  print (svymean( subset ,
                  design = svy.unweighted))
  
  # Again, now get the actual weights when given (via the svy.rake) the weighted EAVE-II population via demographics 
  print (svymean( subset,
                  design = svy.rake))

}

for(x in colnames(df_comorbid.serology_first_6months))  {
  if(!grepl('Q_DIAG', x )){ next}

  calculate_weights(df_comorbid.serology_first_6months,"had_serology_test", x)
  
}


