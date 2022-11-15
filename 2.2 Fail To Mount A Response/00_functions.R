######################################################################
## Title: Serology Analysis - Insufficient IgG Response
##
## Code author: Calum Macdonald (calmacx@gmail.com)
##
## Description: collection of functions used while performing the analysis 
##
##
## > sessionInfo()
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux
## 
## Matrix products: default
## BLAS:   /opt/R/3.6.1/lib64/R/lib/libRblas.so
## LAPACK: /opt/R/3.6.1/lib64/R/lib/libRlapack.so
##
## Last Updated: 02/11/2022
##
######################################################################

#Load libraries that we need
# - dataframe creation and manipulation
library(tibble)
library(dplyr)
library(tidyr)
library(data.table) #!may not be required

# - variable manipulation 
library(forcats) #factors 
library(lubridate) #dates

library(namespace)

eave.clean_serology <- function (df,fix_igg=F){
  df <- df %>%
    mutate(
      IgG = readr::parse_number(test_result_quant),
      QualResult=test_result_qual
    ) 
  
  ### recacluate the IgG/Qual result
  # - this is needed for the PC which used a different assay 
  # --should not greatly affect these results are not many people had a vaccine by this point in time
  if (fix_igg){
    df <- df %>%
      mutate(IgG = ifelse(Sampledate_iso < "2021-04-12",2.6*IgG,IgG),
             QualResult = ifelse(IgG<33.8,"Negative", "Positive")) 
  }
  return (df);
}







