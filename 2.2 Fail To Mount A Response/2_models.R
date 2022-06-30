library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(xtable)


get_qnames <- function(df,n=5){
  
  qnames <- df %>% filter(outcome==1) %>% select(contains("Q_DIA")) %>% 
       summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% gather() %>% 
       filter(value>n)
  qnames <- qnames$key
  return (qnames)
}

get_vars_to_use <- function(model,thres=0.1){

  variables <- model %>% select(-any_of(c("outcome","ageYear","Sex","n_risk_gps"))) %>% names()
  print (variables)
  
  dfp <- NULL
  for (name in variables){
  
    formula <- paste0("outcome ~ ",name)
    print (name)
    
    glmFit <- glm(formula,family=binomial,data=model)
    glmFit.summary <- summary(glmFit)
    print (glmFit.summary)
    
    p = coef(glmFit.summary)[,"Pr(>|z|)"][[2]]
    test = p < 0.05
    temp <- data.frame(name,as.character(lookup[name]),p)
    names(temp) <- c('name','Variable','p-Value')
    dfp = rbind(dfp, temp)
  }
  #write.csv(dfp,paste0(folder,"unadjusted_productother.csv"))
  
  vars_to_use <- as.vector((dfp %>% filter(`p-Value`< thres))$name)
  return (vars_to_use)
}

calc_prod <- function (stage,d1,d2,d3){
  
  d1 = ifelse(d1 == 'Covid-19 Vaccine AstraZeneca',1,2)
  d2 = ifelse(is.na(d2),0,ifelse(d2 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d3 = ifelse(is.na(d3),0,ifelse(d3 == 'Covid-19 Vaccine AstraZeneca',1,2))
  
  d2 <- ifelse(stage<2,0,d2)
  d3 <- ifelse(stage<3,0,d3)
  
  #cat <- d1 + d2*10 + d3*100
  cat <- ifelse(d1==2 & d2==2 & d3==0,0,
                ifelse(stage==2 & !(d1==d2),'mixed 2 doses',
                  ifelse(stage>2 & !(d1==d2 & d1==d3),'mixed 3 doses',
                          paste0(
                            ifelse(d1==2,'Pf','AZ'),
                            ifelse(d2>0,ifelse(d2==2,'-Pf','-AZ'),''),
                            ifelse(d3>0,ifelse(d3==2,'-Pf','-AZ'),'')
                          )
                  )
                )
  )
  #cat <- ifelse(stage>2 & d1==d2 & d1==d3,'same 3 doses',cat)
  cat <- ifelse(cat==0,'Pf-Pf',cat)
  return (cat)
}

constrast_code <- function(model){
  return (model %>% mutate(Sex = as.factor(Sex)#,
                           #simd2020v2_sc_quintile = ifelse(simd2020v2_sc_quintile==5,0,simd2020v2_sc_quintile),
                           #product_time = ifelse(product_time==2,0,product_time)
                           ) %>%
                    mutate(
                      age=ageYear,
                      #product_time=as.factor(product_time),
                      n_risk_gps=as.factor(as.integer(n_risk_gps)-1)) %>%
                    mutate(cat = relevel(as.factor(calc_prod(stage,d1_product,d2_product,d3_product)),'Pf-Pf')) %>%
                    mutate_at(vars(one_of('simd2020v2_sc_quintile')), ~relevel(as.factor(.x),3)) %>%
                    mutate_at(vars(one_of('days_between_d1_d2')), ~as.integer(.x)) %>%
                    #mutate_at(vars(one_of('days_between_d1_d2')), ~relevel(as.factor(cut(.x,breaks=c(-1,25,50,75,1000),right=T,
                    #                                                          labels=c("0-25","25-50","50-75","75+"))),"50-75")) %>%
                              
                    #mutate_at(vars(one_of('days_between_d1_d2')), ~as.factor(cut(.x, breaks = c(-1,2,10,20,50,100,10000), right = T))) %>%
                    #mutate_at(vars(one_of('d1_product')), ~calc_prod(.x)) %>%
                    #mutate_at(vars(one_of('d2_product')),   ~calc_prod(.x) ) %>%
                    #mutate_at(vars(one_of('d3_product')),   ~calc_prod(.x) ) %>%
                    mutate_at(vars(one_of('stage')), ~ factor(ifelse(.x>2,'2+',as.character(.x)),levels=c('2','1','2+'))) %>%
                    #mutate_at(vars(one_of('ageYear')), ~factor(cut(.x, breaks = c(0,20,30,50,70,150), right = T, 
                    #                                                  labels = c("0-20","20-30","30-50","50-70","70+")),
                    #                                           levels=c("50-70","0-20","20-30","30-50","70+")
                    ##)) %>%
                    #mutate_at(vars(one_of('days_since_last_vac')), ~factor(cut(.x, breaks = c(0,30,50,100,150,200,10000), right = T, 
                    #                                               labels = c("15-30","30-50","50-100","100-150","150-200","200+")),
                    #                                           levels=c("50-100","15-30","30-50","100-150","150-200","200+")
                    #)) %>%
                    mutate_at(vars(one_of('ageYear')), ~factor(cut(.x, breaks = c(0,20,40,60,150), right = T, 
                                                                    labels = c("0-19","20-39","40-59","60+")),
                                                               levels=c("40-59","0-19","20-39","60+")
                    )) %>%
                    mutate_at(vars(one_of('Q_BMI')), ~factor(cut(.x, breaks = c(0,20,25,30,150), right = T, 
                                                                             labels = c("0-20","20-25","25-30","30+")),
                                                                         levels=c("20-25","0-20","25-30","30+")
                              ))
  )
}


nrow(df_ana)
modelA <- df_ana %>% select("outcome","ageYear","Sex","Q_BMI","days_since_last_vac","n_risk_gps","stage", 
                            "days_between_d1_d2","days_since_first_measurement",
                              "prior_infection","product_binary","product_time","simd2020v2_sc_quintile",
                              "d1_product","d2_product","d3_product") %>%
                   mutate(age=ageYear,
                          bmi=Q_BMI) %>% 
                   constrast_code() #%>% filter(stage==2)  #%>% select(-d1_product,-d2_product,-d3_product) %>% drop_na() 

#modelA %>% mutate(temp = ifelse(is.na(days_between_d1_d2),0,1)) %>% group_by(temp) %>% summarise(n=n())
#modelA %>% group_by(days_between_d1_d2) %>% summarise(n=n())
#hist((modelA %>% filter(product_binary==0))$days_between_d1_d2, col = 'skyblue3', breaks = 150,xlim=c(0,150))
#hist(modelA$days_between_d1_d2, col = 'skyblue3', breaks = 150,xlim=c(0,150))
#nrow(modelA)
modelA.vars <- get_vars_to_use(modelA)
modelA.formula <- paste0("outcome ~ ageYear + Sex + n_risk_gps + ",paste0(modelA.vars,collapse = " + "))
modelA.formula <- paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + Q_BMI + days_since_last_vac  + product_binary + stage ")
modelA.formula <- paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + Q_BMI + days_since_last_vac + days_between_d1_d2 + product_binary")
#modelA.formula <- paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + Q_BMI + days_since_last_vac  + cat ")
#modelA.formula <- paste0("outcome ~ ageYear + Sex + n_risk_gps + Q_BMI + prior_infection + cat ")
modelA.formula

modelA.formula2 <- paste0("outcome ~ bs(ageYear,df=5,intercept=TRUE) + Sex + n_risk_gps + ",paste0(modelA.vars,collapse = " + "))
modelA.formula2 

modelA.formula3 <- paste0("outcome ~ ageYear")
modelA.formula4 <- paste0("outcome ~ pspline(ageYear)")




modelA_2 <- df_ana_2 %>% select("outcome","ageYear","Sex","days_since_last_vac","n_risk_gps",                            
                              "prior_infection","product_binary","product_time","simd2020v2_sc_quintile") %>%
  constrast_code() %>% drop_na()
modelA_2.formula <- paste0("outcome ~ ageYear + Sex + n_risk_gps + prior_infection + days_since_last_vac + product_time*product_binary")
modelA_2.formula 

modelA_3 <- df_ana_2 %>% select("outcome","ageYear","Sex","days_since_last_vac","n_risk_gps",                            
                                "prior_infection","product_binary","product_time","simd2020v2_sc_quintile") %>%
                          constrast_code() %>% drop_na()
modelA_3.formula <- paste0("outcome ~ ageYear + Sex + prior_infection + days_since_last_vac + n_risk_gps*product_time + product_binary")
modelA_3.formula 


colnames(df_ana)

df_ana %>% group_by(immuno_supp) %>% summarise(n=n())

modelB <- df_ana #%>% filter(n_risk_gps>0) 
nrow(modelB)
qnames <- get_qnames(modelB,n=3)
modelB <- modelB %>% select("outcome","ageYear","Sex","simd2020v2_sc_quintile","Q_BMI","n_risk_gps",qnames,  
                            "shielding","ch_resident","immuno_supp", # "days_between_d1_d2"
                            "days_since_first_measurement",#product_time
                            "prior_infection","days_since_last_vac","product_binary","stage","d1_product","d2_product","d3_product") %>%
                            constrast_code() %>% select(-d1_product,-d2_product,-d3_product) %>% drop_na()# %>%
                            #select(-product_binary,-stage)
                            #select(-cat)#,-age)
#modelB %>% group_by(Q_DIAG_DIABETES_2) %>% summarise(n=n())

colnames(modelB)
modelB %>% group_by(stage) %>% summarise(n=n()) 


modelB <- modelB %>% filter(cat!='AZ-AZ-AZ' & cat!='mixed 2 doses')# & cat!='Pf-Pf-Pf') 
modelB$cat <- droplevels(modelB$cat)
modelB %>%  group_by(cat) %>% summarise(n=n())


#modelB.vars <- get_vars_to_use(modelB %>% select(outcome,qnames),thres=0.5) 
#modelB.formula <- paste0("outcome ~ ageYear + Sex  +  Q_BMI + ",paste0(modelB.vars,collapse = " + "))
modelB.formula <- paste0("outcome ~ ageYear + Sex  +  Q_BMI + ",paste0(qnames,collapse = " + "))
modelB.formula 


nrow(modelB)

modelB.formula2 <- gsub("\\+ product_binary ","",modelB.formula)
#modelB.formula2 <- gsub("\\+ stage ","",modelB.formula2)
modelB.formula2

modelB %>% group_by(cat,product_binary) %>% summarise(n=n())


modelC <- df_ana_2 %>% filter(n_risk_gps>0) 
modelC
qnames <- get_qnames(modelC)

modelC <- modelC %>% select("outcome","ageYear","Sex","Q_BMI","days_since_last_vac","n_risk_gps",qnames,                     
                            "prior_infection","product_binary","product_time","stage",contains("_product")) %>%
                     constrast_code() %>% 
                     select(-n_risk_gps,-age,-product_time,-product_binary,-stage,-days_since_last_vac,-last_vac_product) %>% 
                     select(-prior_infection) %>%
                     drop_na() 

colnames(modelC)
modelC.vars <- colnames(modelC %>% select(-outcome))#get_vars_to_use(modelC)
modelC.vars
modelC.formula <- paste0("outcome ~  ",paste0(modelC.vars,collapse = " + "))
modelC.formula 

modelC %>% group_by(d1_product) %>% summarise(n=n())

library(survival)
modelB.formula2 <- paste0("outcome ~ pspline(ageYear,df=2) + Sex + ",paste0(modelB.vars,collapse = " + "))
modelB.formula2


modelC_1 <- df_ana_2 %>% filter(product_binary==0 & n_risk_gps>0)
modelC_1

qnames <- get_qnames(modelC_1)

modelC_1 <- modelC_1 %>% select("outcome","ageYear","Sex","days_since_last_vac","n_risk_gps",qnames,                         
                                                       "prior_infection","product_time","simd2020v2_sc_quintile") %>%
                   constrast_code() %>% select(-n_risk_gps) %>% drop_na() 

modelC_1.vars <- get_vars_to_use(modelC_1)
modelC_1.formula <- paste0("outcome ~ ageYear + Sex + product_time + ",paste0(modelC_1.vars,collapse = " + "))
modelC_1.formula 



modelC_2 <- df_ana_2 %>% filter(product_binary>0 & n_risk_gps>0)


qnames <- get_qnames(modelC_2)
modelC_2 <- modelC_2 %>% select("outcome","ageYear","Sex","days_since_last_vac","n_risk_gps",qnames,                         
                                "prior_infection","product_time","simd2020v2_sc_quintile") %>%
                     constrast_code() %>% select(-n_risk_gps) %>% drop_na() 

modelC_2





