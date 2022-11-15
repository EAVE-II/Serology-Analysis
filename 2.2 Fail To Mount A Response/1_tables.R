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



get_table_2 <- function(model){
  model$days_since_last_vac = as.factor(cut(model$days_since_last_vac, breaks = c(0,30,50,100,150,200,10000), right = T, 
                                            labels = c("15-30","30-50","50-100","100-150","150-200","200+")))#,
  #levels=c("50-100","15-30","30-50","100-150","150-200","200+"))
  
  model$days_since_first_measurement = as.factor(cut(model$days_since_first_measurement, breaks = c(0,90,180,270,360,1000), right = T, 
                                                     labels = c("0-90","90-180","180-270","270-360","360-1000")))#,
  #levels=c("100-200","0-100","100-300","300-400"))
  
  
  qnames <- model %>% get_qnames(n=0) %>% paste0(collapse = " + ")
  qnames
  
  formula <- paste0('outcome ~ ageYear + Sex + Q_BMI + simd2020v2_sc_quintile + n_risk_gps + ',qnames,' + days_since_last_vac 
                  + shielding + immuno + prior_infection + days_since_first_measurement + cat')
  
  formula <- paste0('outcome ~ ageYear + Sex + Q_BMI + simd2020v2_sc_quintile + n_risk_gps + immuno ')
  
  #+ ch_resident
  #formula <- 'outcome ~ outcome_hospitalisation'
  #formula <- 'outcome ~ days_since_first_measurement'
  
  #tab <- perform_glm(model,formula,do_adjusted = FALSE)
  
  
  tab <- tab %>% replace(is.na(.), 1) %>% 
    mutate(or=paste0(round(uOR,2)," (", round(uLCL,2), " - ", round(uUCL,2), " ); ",ifelse(p=='<0.05',paste0('p',p),paste0('p=',p)) )) %>% 
    mutate(impaired=paste0(fail," (",round(percentage,2)," %)")) %>% 
    select(total,impaired,or)
  return (tab);
}

days_to_factors <- function(model){
  
  model$days_since_last_vac = as.factor(cut(model$days_since_last_vac, breaks = c(0,30,50,100,150,200,10000), right = T, 
                                            labels = c("15-30","30-50","50-100","100-150","150-200","200+")))
  
  model$days_since_first_measurement = as.factor(cut(model$days_since_first_measurement, breaks = c(-1,90,180,270,360,10000), right = T, 
                                                     labels = c("0-90","90-180","180-270","270-360","360-1000")))
  return(model);
}

get_table_1 <- function(model,vars){
  tab <- list()
  for (name in vars){
    orig <- name
    print (orig)
    temp <- model %>% group_by(!!sym(name)) %>% summarise(orig=name,name=lookup[name],N=n(),NF=sum(outcome==1)) %>% ungroup %>% 
      rename(var=!!sym(name)) %>% 
      mutate(pN = round(100*N/sum(N),2),pNF=round(100.*NF/N,2)) %>% 
      mutate(N=paste0(N," (",pN,")"), NF=paste0(NF," (",pNF,")")) %>%
      mutate(var=ifelse(is.na(lookup[paste0(orig,var)]),as.character(var),lookup[paste0(orig,var)])) %>%
      select(name,var,N,NF) 
    #mutate(var=!!sym(name)) %>%  select(-!!sym(name)) %>% select(outcome,name,var,n) %>%
    #mutate(name=ifelse(is.na(name),orig,name,)) 
    tab[[name]] <- temp
  }
  tab1 <- rbindlist(tab) %>% as_tibble
  return (tab1);
}


pc.modelA.tab1_vars <- pc.modelA %>% select(ageYear,Sex,simd2020v2_sc_quintile,Q_BMI,
                                            immuno,shielding,ch_resident,
                                            contains("Q_DIAG"),n_risk_gps,
                                            days_since_last_vac,days_since_first_measurement,prior_infection,cat) %>% colnames()
pc.modelA.tab1_vars

pc.tab1 <- pc.modelA %>% days_to_factors %>% select(outcome,pc.modelA.tab1_vars) 
pc.tab1
pc.tab1 <- get_table_1(pc.tab1,pc.modelA.tab1_vars)
pc.tab1


bd.modelA.tab1_vars <- bd.modelA %>% select(ageYear,Sex,simd2020v2_sc_quintile,Q_BMI,
                                            immuno,shielding,ch_resident,
                                            contains("Q_DIAG"),n_risk_gps,
                                            days_since_last_vac,days_since_first_measurement,prior_infection,cat) %>% colnames()

bd.modelA.tab1_vars 
bd.tab1 <- bd.modelA %>% days_to_factors %>% select(outcome,bd.modelA.tab1_vars) 
bd.tab1 <- get_table_1(bd.tab1,bd.modelA.tab1_vars)


pc.tab1 %>% merge(bd.tab1,by=c('name','var'))


View(tab1)


tab_pc <- get_table_2(model)
tab_pc

tab_bd <- load_bd() %>% constrast_code() %>% get_table_2()
tab_bd


tab_pc$id  <- 1:nrow(tab_pc)
tab_pc

tab <- tab_pc %>% merge(tab_bd,by=0,all=T) %>% arrange(id) %>% select(-id) %>% replace(is.na(.), '') %>%
       replace(. == '1 (1 - 1 ); p=1','')
tab
View(tab)


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




