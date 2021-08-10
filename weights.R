
library(dplyr)
library("tidyr")
library(tibble)
library(survey)

# load the serology dataset (first 6montsh of the pandemic)
df_serology <- readRDS("/conf/EAVE/GPanalysis/data/Serology_all.rds")
# load the EAVE-II demographics 
df_demographics <- readRDS("/conf/EAVE/GPanalysis/data/EAVE_demographics_SK.rds")
# load the comorbidities from QCOviD
df_comorbid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/Qcovid_update_Jan21.rds")

# print to see how many rows/entries of data we have for each dataset
# demographics has ~5.8 million > scotish population 
# comorbidies can have multiple rows for a single study 
nrow(df_serology)
nrow(df_demographics)
nrow(df_comorbid)


# print to see the column names 
colnames(df_serology)
colnames(df_demographics)
colnames(df_comorbid)



# Perform some transformations on the Serology dataset: 
# * only select the study ID, the Age, the Sex and the quantative serology result
# * map Sex values M--> 1 (integer) and F--> 2 (integer). This is needed by the survey package
df_serology_analysis <- df_serology %>%
  #  select(EAVE_LINKNO,Age,Sex,Result_quant) %>%
  select(Age,Sex) %>%
  mutate(Sex = recode(Sex, `M` = 1L,
                           `F` = 2L )) 


df_serology_analysis$Sex <- as.integer(df_serology_analysis$Sex)

# convert age into an integer as it's a char
df_serology_analysis$Age <- as.integer(df_serology_analysis$Age)
# print to see this new dataframe
head(df_serology_analysis)

# Use the survey package to create the initial unweighted data
df_serology_analysis.svy.unweighted <- svydesign(ids=~1, 
                                                 data=df_serology_analysis,
                                                 weights=NULL)
# Print out of the number of rows to make sure things still make sense 
nrow(df_serology_analysis.svy.unweighted)

# Now open the EAVE-II demographics
# * extract the sex column
# * count how many unique values there are of each Sex
# * rename the count  column (n) to be called Freq (as needed by survey package)
# * convert the Sex into an integer again

# !note: this is not 100% correct because nrows(df_demographics) > population of scotland
#        we need to account for this bias somehow 
sex_dist <- count(df_demographics,Sex) %>% 
            rename(Freq=n) %>% 
            mutate(Sex = recode(Sex, `M` = 1L,
                                     `F` = 2L )) 
# Print to see the distribution
sex_dist



# Use the survey package to perform a rake
# "Raking uses iterative post-stratification to match marginal distributions of a survey sample to known population margins."
# https://www.rdocumentation.org/packages/survey/versions/4.1-1/topics/rake
df_serology_analysis.svy.rake <- rake(design = df_serology_analysis.svy.unweighted,
                                      sample.margins = list(~Sex),
                                      population.margins = list(sex_dist))


# Using the survey package svymean, get the mean/SE weights for Age & Sex when there is no weighting
svymean(df_serology_analysis,
        design = df_serology_analysis.svy.unweighted)
# Again, now get the actual weights when given (via the svy.rake) the weighted EAVE-II population via demographics 
svymean(df_serology_analysis,
        design = df_serology_analysis.svy.rake)


# As with sex, do the same with age
age_dist <- count(df_demographics,ageYear) %>% 
  rename(Freq=n, Age=ageYear)
age_dist$Age = as.integer(age_dist$Age)

# Print to see the distribution
# !note: showing entries with age=-1 and people as old as 137!!
age_dist


# In consistent values was causing all sorts of problems
# Temp fix - list the age_dist on a list of ages that are actually in the dataset
age_list = sort(unique(df_serology_analysis$Age))
age_dist <- age_dist %>%
            filter(Age %in% age_list)


# Use the survey package to perform a rake
df_serology_analysis.svy.age_rake <- rake(design = df_serology_analysis.svy.unweighted,
                                      sample.margins = list(~Age),
                                      population.margins = list(age_dist))

# get the weights again
svymean(df_serology_analysis,
        design = df_serology_analysis.svy.age_rake)


# try both age and sex at the same time
df_serology_analysis.svy.combo_rake <- rake(design = df_serology_analysis.svy.unweighted,
                 sample.margins = list(~Sex,~Age),
                 population.margins = list(sex_dist, age_dist))

svymean(df_serology_analysis,
        design = df_serology_analysis.svy.combo_rake)

svymean(df_serology_analysis,
        design = df_serology_analysis.svy.unweighted)


# Print the unique comorbidity values 
unique(df_comorbid$cluster)

# Test playing with a new dataframe for people classified as having diabetes type-1
# !note: this assumes df_comorbid has everyone in EAVE-II, which is likely a false assumption
df_comorbid_diabetes_1 <- df_comorbid %>%
                          select(EAVE_LINKNO,cluster)  %>%
                          mutate(has_condition = if_else( cluster == 'Q_DIAG_DIABETES_1', 1L, 2L)) %>%
                          select(EAVE_LINKNO,has_condition)
head(df_comorbid_diabetes_1)


df_serology_analysis_diabetes <- df_serology %>%
                                 select(EAVE_LINKNO,Age,Sex) %>%
                                 left_join(df_comorbid_diabetes_1) %>%
                                 select(Age,has_condition) %>%
                                 mutate(Age = as.integer(Age)) %>%
                                 na.omit()
head(df_serology_analysis_diabetes)


df_serology_analysis_diabetes.svy.unweighted <- svydesign(ids=~1, 
                                                 data=df_serology_analysis_diabetes,
                                                 weights=NULL)

svymean(df_serology_analysis_diabetes,
        design = df_serology_analysis_diabetes.svy.unweighted)

condition_dist <- count(df_comorbid_diabetes_1,has_condition) %>% 
                  rename(Freq=n)
         

df_serology_analysis_diabetes.svy.rake <- rake(design = df_serology_analysis_diabetes.svy.unweighted,
                                       sample.margins = list(~has_condition),
                                       population.margins = list(condition_dist))


svymean(df_serology_analysis_diabetes,
        design = df_serology_analysis_diabetes.svy.rake)
