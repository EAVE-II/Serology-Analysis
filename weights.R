
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

# get unique people (studies) that are in the serology dataset (first 6 months)
df_serology_first_6months.unique_people <- df_serology_first_6months %>% 
                                           distinct(EAVE_LINKNO) %>%
                                           add_column(had_serology_test=1L)

# get unique people (studies) that are in the serology dataset (next 6 months+)
# !note: currently this is not used
df_serology_Oct20May21.unique_people <- df_serology_Oct20May21 %>% 
                                        distinct(EAVE_LINKNO) %>%
                                        add_column(had_serology_test=1L)

# join the comorbidity table with the unique people serology table
# replace nan values with 0L (0 integer) for people who are in the serology dataset but not in the comorbid+demographics 
df_comorbid.serology_first_6months <- df_comorbid %>% 
                                      left_join(df_serology_first_6months.unique_people) %>%
                                      replace(is.na(.),0L) 


# Build a function to calculate some weights
# - dataset: is the RDataFrame for the comorbid+demographics+serology dataset (will be the first/next 6months, or a combination)
# - subset_filter: name of flag to filter the dataset into a subset (will be comorbid+demographics for people who have serology data)
# - comorbidity: name of the comorbidity to calculate weights for
calculate_weights <- function (dataset, subset_filter, comorbidity) {
  
  #first calculate the frequency of occurrence of a comorbidity in the full dataset
  freq <- count(dataset,!!as.name(comorbidity) ) %>%
          rename(Freq=n)

  #get the subset from the dataset and just select the comorbidity we are interested in 
  subset <- dataset %>% 
             filter(!!as.name(subset_filter) > 0 ) %>%
             select(!!as.name(comorbidity))
  
  # calculate the frequency of comorbidities for this too
  freq.subset <- count(subset, !!as.name(comorbidity))  %>% rename(Freq=n)

  # do a poor man's calculation of the weights:
  # print the fraction of those who do/don't have the comorbidity in the full dataset
  print(freq$Freq / sum(freq$Freq))
  # print the fraction of those who do/don't have the comorbidity in the subset dataset
  print(freq.subset$Freq / sum(freq.subset$Freq))
  
  # Use the survey package to create the initial unweighted data
  # !NOTE - to-do, this may need to include weighting of the full dataset to match the Scottish population
  #         the full dataset has ~5.8 million people and needs to be pre-weights to match the ~5.5 million people 
  #         therefore we may need to take this into account using df_pop_estimate
  svy.unweighted <- svydesign(ids=~1, 
                              data=subset,
                              weights=NULL)
  
  # Use the rake method https://r-survey.r-forge.r-project.org/survey/html/rake.html
  # to weight for the actual frequencies of comorbidities
  # !NOTE - just doing this for one comorbidity at a time, but could be done for multiple
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

# loop over all column names in the main dataset
for(x in colnames(df_comorbid.serology_first_6months))  {
  # pick out the comorobidities that we want to test for (i.e. skip other comorbidities
  if(!grepl('Q_DIAG', x )){ next}
  #calculate the weights for these
  calculate_weights(df_comorbid.serology_first_6months,"had_serology_test", x)
  
}


