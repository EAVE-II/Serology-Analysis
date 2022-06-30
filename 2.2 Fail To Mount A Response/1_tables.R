library(xtable)
library(tidyverse)



create_table <- function(df,variable) {

  tab1 <- df %>% filter(outcome==0) %>% 
          group_by(!!as.name(variable)) %>% 
          summarise(n = n()) %>%
          mutate(per = round(100.*n / sum(n),2))
  
  tab2 <- df %>% filter(outcome==1) %>% 
          group_by(!!as.name(variable)) %>% 
          summarise(n = n()) %>%
          mutate(per = round(100.*n / sum(n),2) ) 
  
                   
  tab <- merge(tab1, tab2, by=variable, all = T) %>%
         rename( !!as.character(lookup[variable]) := !!as.name(variable))
  
  file <- paste0(folder,variable,".csv")
  write.csv(tab,file,row.names=FALSE)
                                           
  return(tab)
}

create_table_mean <- function(df,variable) {
  
  tab1 <- df %>% filter(outcome==0) 
  tab1 <- tab1[variable] %>% filter(!is.na(!!as.name(variable))) %>% 
          summarise_if(is.numeric, list(mean,sd))
  
  tab2 <- df %>% filter(outcome==1) 
  tab2 <- tab2[variable] %>% filter(!is.na(!!as.name(variable))) %>%
          summarise_if(is.numeric, list(mean,sd))
  
  new <- c("Mean")
  tab <- cbind(new,cbind(tab1, tab2))
  names(tab)[1] <- variable
  
  
  file <- paste0(folder,variable,"_mean.csv")
  write.csv(tab,file,row.names=FALSE)
  

  return(tab)
 
}


tab <- create_table(df_ana,"simd2020v2_sc_quintile")
tab

tab <- create_table(df_ana_2,"n_risk_gps")
tab

tab <- create_table(df_ana_2,"Sex")
tab

tab <- create_table(df_ana_2,"last_vac_product")
tab



colnames(df_ana_2)
tab <- create_table_mean(df_ana_2,"ageYear")
tab

create_table_mean(df_ana_2,"days_since_infection")
create_table_mean(df_ana_2,"days_since_last_vac")




tab <- create_table(df_ana_2,"time_of_last_vac")
tab

temp <- df_ana_2 %>% mutate(prior_infection = ifelse(is.na(days_since_infection),0,1))
tab <- create_table(temp,"prior_infection")
rm(temp)
tab


q_names <- colnames(df_ana_2 %>% select(contains('Q_DIAG')))

for (name in q_names)
{
  print (create_table(df_ana_2,name))
}




