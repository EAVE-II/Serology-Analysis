library(tibble)
library(dplyr)
library(ggplot2)
library(tidyr)
library(xtable)
library(phsstyles)


get_qnames <- function(df,var='outcome',n=5){
  
  qnames <- df %>% filter(!!sym(var)==1) %>% select(contains("Q_DIA")) %>% 
       summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) 
  qnames <- qnames %>% gather() %>% 
       filter(value>=n)
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

calc_prod <- function (stage,d1,d2,d3,d4){
  
  d1 = ifelse(d1 == 'Covid-19 Vaccine AstraZeneca',1,2)
  d2 = ifelse(is.na(d2),0,ifelse(d2 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d3 = ifelse(is.na(d3),0,ifelse(d3 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d4 = ifelse(is.na(d4),0,ifelse(d4 == 'Covid-19 Vaccine AstraZeneca',1,2))
  
  
  d2 <- ifelse(stage<2,0,d2)
  d3 <- ifelse(stage<3,0,d3)
  d4 <- ifelse(stage<4,0,d4)
  
  #cat <- d1 + d2*10 + d3*100
  cat <- ifelse(d1==2 & d2==2 & d3==0 & d4==0,0,
                ifelse(stage==2 & !(d1==d2),'mixed 2 doses',
                  #ifelse(stage>3,'mixed 4 doses',
                  ifelse(stage>2 & !(d1==d2 & d1==d3),'mixed 3 doses',
                          paste0(
                            ifelse(d1==2,'Pf','AZ'),
                            ifelse(d2>0,ifelse(d2==2,'-Pf','-AZ'),''),
                            ifelse(d3>0,ifelse(d3==2,'-Pf','-AZ'),'')#,
                            #ifelse(d4>0,ifelse(d4==2,'-Pf','-AZ'),'')
                          )
                  )
                #)
            )
  )

  cat <- ifelse(cat==0,'Pf-Pf',cat)
  return (cat)
}

contrast_code <- function(model){
  return (model %>% mutate(Sex = as.factor(Sex)#,
                           #simd2020v2_sc_quintile = ifelse(simd2020v2_sc_quintile==5,0,simd2020v2_sc_quintile),
                           #product_time = ifelse(product_time==2,0,product_time)
                           ) %>%
                    mutate(
                      age=ageYear,
                      #product_time=as.factor(product_time),
                      n_risk_gps=as.factor(as.integer(n_risk_gps)-1)) %>%
                    mutate(cat = relevel(as.factor(calc_prod(stage,d1_product,d2_product,d3_product,d4_product)),'Pf-Pf')) %>%
                    mutate_at(vars(one_of('simd2020v2_sc_quintile')), ~relevel(as.factor(.x),'3')) %>%
                    mutate_at(vars(one_of('days_between_d1_d2')), ~as.integer(.x)) %>%
                    #mutate_at(vars(one_of('days_between_d1_d2')), ~relevel(as.factor(cut(.x,breaks=c(-1,25,50,75,1000),right=T,
                    #                                                          labels=c("0-25","25-50","50-75","75+"))),"50-75")) %>%
                              
                    #mutate_at(vars(one_of('days_between_d1_d2')), ~as.factor(cut(.x, breaks = c(-1,2,10,20,50,100,10000), right = T))) %>%
                    #mutate_at(vars(one_of('d1_product')), ~calc_prod(.x)) %>%
                    #mutate_at(vars(one_of('d2_product')),   ~calc_prod(.x) ) %>%
                    #mutate_at(vars(one_of('d3_product')),   ~calc_prod(.x) ) %>%
                    mutate_at(vars(one_of('stage')), ~ factor(ifelse(.x>3,'4+',as.character(.x)),levels=c('2','1','3','4+'))) %>%
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
                                                               levels=c("40-59","0-19","20-39","60+"))) %>%
                    #mutate_at(vars(one_of('Q_BMI')), ~factor(cut(.x, breaks = c(0,20,25,30,150), right = T, 
                    #                                                         labels = c("0-20","20-25","25-30","30+")),
                    #                                                     levels=c("20-25","0-20","25-30","30+"))) %>% 
                    mutate_at(vars(one_of('Q_BMI')), ~relevel(cut(ifelse((is.na(.x) | .x<5 | .x>80 ),-1,.x), breaks = c(-2,0,20,25,30,80), right = T, 
                                                                        labels = c("Unknown","0-20","20-25","25-30","30+")),ref='20-25')
                    ) %>%
                   #mutate(immuno = strtoi(paste0(immuno_supp,severely_immuno_supp,Q_DIAG_IMMU))) %>%
                   mutate(immuno = strtoi(paste0(immuno_supp,severely_immuno_supp,0))) %>%
                   mutate(immuno = as.factor(ifelse(immuno==1,0,immuno))) %>%
                   select(-immuno_supp,-severely_immuno_supp,-Q_DIAG_IMMU)
  )
}


create_modelA <- function(df_ana){
  modelA <- df_ana %>% 
                   mutate(age=ageYear,
                          bmi=Q_BMI) %>% 
                   contrast_code()  %>% 
                   filter(cat!='AZ-AZ-AZ' & cat!='mixed 2 doses') %>% droplevels 

  
  #cats <- (modelA %>% group_by(cat) %>% summarise(n=n()) %>% mutate(cat=as.character(cat)) %>% arrange(desc(n)))$cat
  cats <- c("Pf-Pf","AZ","Pf","AZ-AZ","mixed 3 doses","Pf-Pf-Pf")#,"mixed 4 doses")
  modelA <- modelA %>% mutate(cat=factor(cat,levels=cats)) %>% droplevels
  
  return (modelA)

}

pc.modelA <- create_modelA(df_ana_pc)
pc.modelA.T <- pc.modelA %>% mutate(outcome=ifelse(IgG<4.9,1,0))
pc.modelA.L <- pc.modelA %>% mutate(outcome=ifelse(IgG<50,1,0))
colnames(pc.modelA)


#temp <- pc.modelA %>% select(EAVE_LINKNO,Sampledate_iso) %>% 
#  left_join(df_smr01 %>% filter(ADMISSION_COVID!=1) %>% select(EAVE_LINKNO,ADMISSION_DATE) )%>% 
#  group_by(EAVE_LINKNO) %>% summarise(nprior=sum(abs(as.integer(ADMISSION_DATE-Sampledate_iso,units='days')<50)))
#pc.modelA <- pc.modelA %>% left_join(temp) %>% mutate(recent_hosp = ifelse(is.na(nprior),0,1)) %>% select(-nprior)



#tables here

#ch_resident
#recent_hosp
pc.modelA.vars <- pc.modelA %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding,ch_resident) %>% get_vars_to_use
pc.modelA.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(pc.modelA.vars,collapse = " + "))
pc.modelA <- pc.modelA %>% select(stage,all.vars(as.formula(pc.modelA.gam.formula)))
pc.modelA.gam.formula


pc.modelA.T.vars <- pc.modelA.T %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding,ch_resident) %>% get_vars_to_use
pc.modelA.T.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(pc.modelA.T.vars,collapse = " + "))
pc.modelA.T <- pc.modelA.T %>% select(stage,all.vars(as.formula(pc.modelA.T.gam.formula)))
pc.modelA.T.gam.formula

pc.modelA.L.vars <- pc.modelA.L %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding,ch_resident) %>% get_vars_to_use
pc.modelA.L.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(pc.modelA.L.vars,collapse = " + "))
pc.modelA.L <- pc.modelA.L %>% select(stage,all.vars(as.formula(pc.modelA.L.gam.formula)))
pc.modelA.L.gam.formula




bd.modelA <- create_modelA(df_ana_bd)
bd.modelA <- bd.modelA %>% mutate(outcome = ifelse(is.na(outcome),1,outcome)) %>%
  mutate(cat = fct_collapse(cat, `3 doses`=c("mixed 3 doses","Pf-Pf-Pf")))

bd.modelA.T <- bd.modelA %>% mutate(outcome=ifelse(IgG<0.12,1,0))
bd.modelA.L <- bd.modelA %>% mutate(outcome=ifelse(IgG<5,1,0))

nrow(bd.modelA)
#temp <- bd.modelA %>% select(EAVE_LINKNO,Sampledate_iso) %>% 
#  left_join(df_smr01 %>% filter(ADMISSION_COVID!=1) %>% select(EAVE_LINKNO,ADMISSION_DATE) )%>% 
#  group_by(EAVE_LINKNO) %>% summarise(nprior=sum(abs(as.integer(ADMISSION_DATE-Sampledate_iso,units='days')<50)))#
#bd.modelA <- bd.modelA %>% left_join(temp) %>% mutate(recent_hosp = ifelse(is.na(nprior),0,1)) %>% select(-nprior)


#ch_resident
bd.modelA.vars <- bd.modelA %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding) %>% get_vars_to_use
#recent_hosp
bd.modelA.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(bd.modelA.vars,collapse = " + "))

bd.modelA <- bd.modelA %>% select(stage,all.vars(as.formula(bd.modelA.gam.formula)))

bd.modelA.T.vars <- bd.modelA.T %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding) %>% get_vars_to_use
bd.modelA.T.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(bd.modelA.T.vars,collapse = " + "))
bd.modelA.T <- bd.modelA.T %>% select(stage,all.vars(as.formula(bd.modelA.T.gam.formula)))
bd.modelA.T.gam.formula

bd.modelA.L.vars <- bd.modelA.L %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding) %>% get_vars_to_use
bd.modelA.L.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat  + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + simd2020v2_sc_quintile + ",paste0(bd.modelA.L.vars,collapse = " + "))
bd.modelA.L <- bd.modelA.L %>% select(stage,all.vars(as.formula(bd.modelA.L.gam.formula)))
bd.modelA.L.gam.formula





nthres <- 2

pc.modelB <- create_modelA(df_ana_pc %>% filter(n_risk_gps!=0)) %>%
             mutate(outcome=ifelse(IgG<4.9,1,0))

#levels <- c(0.0,0.70)
#igg_quantile <- quantile(pc.modelB$IgG,levels)
#pc.modelB <- pc.modelB %>% mutate(outcome=ifelse(IgG>igg_quantile[[2]],1,0))


pc.modelB.vars <- pc.modelB %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding,ch_resident) %>% get_vars_to_use
pc.modelB.qnames <- get_qnames(pc.modelB,n=nthres)
pc.modelB.qnames

pc.modelB <- pc.modelB %>% select(stage,all.vars(as.formula(pc.modelB.gam.formula)))
pc.modelB.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(pc.modelB.qnames,collapse = " + "),'+',paste0(pc.modelB.vars,collapse = " + "))
#pc.modelB.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(pc.modelB.qnames,collapse = " + "),'+',paste0(pc.modelB.vars,collapse = " + "))

#pc.modelB.gam.formula <- paste0("outcome ~ s(age) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(pc.modelB.qnames,collapse = " + "),'+',paste0(pc.modelB.vars,collapse = " + "))
#pc.modelB.gam.formula <- paste0("outcome ~ ageYear + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI + ",paste0(pc.modelB.qnames,collapse = " + "),'+',paste0(pc.modelB.vars,collapse = " + "))

pc.modelB <- pc.modelB %>% select(stage,all.vars(as.formula(pc.modelB.gam.formula)))
pc.modelB.gam.formula




bd.modelB <- create_modelA(df_ana_bd %>% filter(n_risk_gps!=0))%>%
  mutate(outcome=ifelse(IgG<0.11,1,0))

#levels <- c(0.0,0.70)
#igg_quantile <- quantile(bd.modelB$IgG,levels)
#bd.modelB <- bd.modelB %>% mutate(outcome=ifelse(IgG>igg_quantile[[2]],1,0))

bd.modelB <- bd.modelB %>% mutate(outcome = ifelse(is.na(outcome),1,outcome)) %>%
  mutate(cat = fct_collapse(cat, `3 doses`=c("mixed 3 doses","Pf-Pf-Pf")))


bd.modelB.vars <- bd.modelB %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno,shielding) %>% get_vars_to_use
bd.modelB.qnames <- get_qnames(bd.modelB,n=nthres)
bd.modelB.qnames
bd.modelB.vars


bd.modelB.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(bd.modelB.qnames,collapse = " + "),'+',paste0(bd.modelB.vars,collapse = " + "))
#bd.modelB.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(bd.modelB.qnames,collapse = " + "),'+',paste0(bd.modelB.vars,collapse = " + "))

#bd.modelB.gam.formula <- paste0("outcome ~ s(age) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(bd.modelB.qnames,collapse = " + "),'+',paste0(bd.modelB.vars,collapse = " + "))
#bd.modelB.gam.formula <- paste0("outcome ~ ageYear + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI  + ",paste0(bd.modelB.qnames,collapse = " + "),'+',paste0(bd.modelB.vars,collapse = " + "))
bd.modelB <- bd.modelB %>% select(stage,all.vars(as.formula(bd.modelB.gam.formula)))
bd.modelB.gam.formula





#modelB <- modelB %>% mutate(outcome=ifelse(IgG<70,1,0))

if(F){
  modelB <- modelB %>% mutate(outcome = ifelse(is.na(outcome),1,outcome)) %>%
    mutate(cat = fct_collapse(cat, `3 doses`=c("mixed 3 doses","Pf-Pf-Pf")))
}
modelB

qnames <- get_qnames(modelB,n=5)
qnames

#ch_resident
modelB.vars <- modelB %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno_supp,shielding,ch_resident) %>% get_vars_to_use
modelB.vars

#s(days_since_first_measurement)
modelB.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + simd2020v2_sc_quintile + Q_BMI + ch_resident + ",paste0(qnames,collapse = " + "),'+',paste0(modelB.vars,collapse = " + "))
modelB.gam.formula



modelB.tight <- modelB %>% mutate(outcome=ifelse(IgG<15,1,0))
qnames <- get_qnames(modelB.tight,n=1)


modelB.tight.vars <- modelB.tight %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno_supp,shielding,ch_resident) %>% get_vars_to_use

#s(days_since_first_measurement)
modelB.tight.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat +  + s(days_since_last_vac) + Sex + ",paste0(qnames,collapse = " + "),'+',paste0(modelB.tight.vars,collapse = " + "))
modelB.tight.gam.formula



modelC.A <- modelA %>% filter(stage==3 | stage=='4+') %>% droplevels()
modelC.A.vars <- modelC.B %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno_supp,shielding,ch_resident) %>% get_vars_to_use
modelC.A.vars


modelC.A.gam.formula <- paste0("outcome ~ s(age,by=cat) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + n_risk_gps + ",paste0(modelC.A.vars,collapse = " + "))
modelC.A.gam.formula



modelC.B <- modelB %>% filter(stage==3 | stage=='4+') %>% droplevels()
modelC.B.vars <- modelC.B %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno_supp,shielding,ch_resident) %>% get_vars_to_use
modelC.B.vars

qnames <- get_qnames(modelC.B,n=1)
qnames



modelC.B.gam.formula <- paste0("outcome ~ s(age,by=cat) + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + ",paste0(qnames,collapse = " + "),'+',paste0(modelC.B.vars,collapse = " + "))
modelC.B.gam.formula






modelD <- modelB %>% filter(stage==3 | stage=='4+') %>% filter(Q_DIAG_BLOOD_CANCER==1) %>% droplevels()
modelD
modelD.vars <- modelD %>% select(outcome,prior_infection,simd2020v2_sc_quintile,immuno_supp,shielding) %>% get_vars_to_use
modelD.vars

qnames <- get_qnames(modelD,n=1)
qnames



modelD.gam.formula <- paste0("outcome ~ s(age,by=cat) + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + Q_BMI + ",paste0(qnames,collapse = " + "))#,'+',paste0(modelD.vars,collapse = " + "))
modelD.gam.formula






modelC.gam.formula <- "outcome ~ s(age) + cat + s(days_since_first_measurement) + s(days_since_last_vac) + Sex + n_risk_gps + prior_infection + simd2020v2_sc_quintile + immuno_supp + shielding "







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

modelB <- df_ana %>% filter(n_risk_gps>0) 
nrow(modelB)
qnames <- get_qnames(modelB,n=0)
modelB <- modelB %>% select("outcome","ageYear","Sex","simd2020v2_sc_quintile","Q_BMI","n_risk_gps",qnames,  
                            , # "days_between_d1_d2"
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




colors <- c("Negative" = phs_colours(c("phs-blue")),"Equivocal" = phs_colours(c("phs-magenta")),"Positive" = phs_colours(c("phs-purple")))


modelC %>% group_by(test_result_qual) %>% summarise(n=n())

df_ana %>% group_by(stage,outcome) %>% summarise(n=n())

df <- df_ana %>% filter(stage==1)

nrow(df_ana)


p1 <- ggplot(df, aes(x=IgG)) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Negative'), aes(fill='Negative')) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Positive'), aes(fill='Positive')) +
  labs(title='Dose 1',x="Antibody Levels [IgG]", y="Number of Samples", fill="Result") +
  scale_fill_manual(values = colors) +
  theme(aspect.ratio = 0.5)  +
  scale_y_log10() +
  scale_x_log10()

p1

df <- df_ana %>% filter(stage==2)

p2 <- ggplot(df, aes(x=IgG)) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Negative'), aes(fill='Negative')) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Positive'), aes(fill='Positive')) +
  labs(title='Dose 2',x="Antibody Levels [IgG]", y="Number of Samples", fill="Result") +
  scale_fill_manual(values = colors) +
  theme(aspect.ratio = 0.5)  +
  scale_y_log10() +
  scale_x_log10()

df <- df_ana %>% filter(stage==3)
p3 <- ggplot(df, aes(x=IgG)) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Negative'), aes(fill='Negative')) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Positive'), aes(fill='Positive')) +
  labs(title='Dose 3',x="Antibody Levels [IgG]", y="Number of Samples", fill="Result") +
  scale_fill_manual(values = colors) +
  theme(aspect.ratio = 0.5)  +
  scale_y_log10() +
  scale_x_log10()


df <- df_ana %>% filter(stage==4)
df

p4 <- ggplot(df, aes(x=IgG)) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Negative'), aes(fill='Negative')) +
  geom_histogram(position = "stack", bins=30, data=subset(df,test_result_qual == 'Positive'), aes(fill='Positive')) +
  labs(title='Dose 4',x="Antibody Levels [IgG]", y="Number of Samples", fill="Result") +
  scale_fill_manual(values = colors) +
  theme(aspect.ratio = 0.5)  +
  scale_y_log10() +
  scale_x_log10()




grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2)
grid.arrange(p1,p2,p3,ncol=2,nrow=2)



df %>% group_by(stage,outcome) %>% summarise(n=n())



qnames <- df_ana %>% filter(stage==3) %>% get_qnames(n=1)
qnames

name <- qnames[[1]]
name

data <- NULL
for (name in qnames){
  df <- df_ana %>% filter(!!as.name(name)==1) %>% group_by(outcome,stage) %>% summarise(n=n())%>% ungroup %>%
          group_by(stage) %>% summarise(outcome=outcome,total=sum(n),n=n,e=1/sqrt(n),p=100*n/sum(n)) %>% 
           mutate(pu=p+p*e,pd=p-p*e) %>% filter(outcome==1) %>% select(-outcome) %>%
          mutate(name=name)
  data <- data %>% rbind(df)
  
}


df <- df_ana %>% group_by(outcome,stage) %>% summarise(n=n())%>% ungroup %>%
  group_by(stage) %>% summarise(outcome=outcome,total=sum(n),n=n,e=1/sqrt(n),p=100*n/sum(n)) %>% 
  mutate(pu=p+p*e,pd=p-p*e) %>% filter(outcome==1) %>% select(-outcome) %>%
  mutate(name='All')
df
data <-  df %>% rbind(data) 
data

tab <- data %>% mutate(name=lookup[name],n=paste0(n," (",round(p,2)," %)")) %>% select(name,stage,total,n)
tab
View(tab)



data <- NULL
for (lev in levels(df_ana$n_risk_gps) ){
  df <- df_ana %>% filter(n_risk_gps==lev) %>% group_by(outcome,stage) %>% summarise(n=n())%>% ungroup %>%
          group_by(stage) %>% summarise(outcome=outcome,total=sum(n),n=n,e=1/sqrt(n),p=100*n/sum(n)) %>% 
           mutate(pu=p+p*e,pd=p-p*e) %>% filter(outcome==1) %>% select(-outcome) %>%
          mutate(name=paste0(lev,' Risks'))
  data <- data %>% rbind(df)
  
}
data
tab <- data %>% mutate(name=name,n=paste0(n," (",round(p,2)," %)")) %>% select(name,stage,total,n)
tab
View(tab)







data %>% ggplot(aes(x=stage,y=p,color=name, fill=name, weight=e)) +
  geom_point() +
  geom_errorbar(aes(ymin=pd, ymax=pu), width=.0) +
  stat_smooth(method = "lm",formula = y ~ exp(-x), se = F) 









