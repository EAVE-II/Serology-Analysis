library(data.table)
library(forcats)
library(broom)
library(survival)
library("survminer")
library(forestplot)
library(ggtext)
library(kableExtra)
library(rlang)
library(visreg)

get_forest_df <- function(model,data){
  main <- "Hazard ratio"
  cpositions <- c(0.02, 0.22, 0.4)
  fontsize <- 0.7
  refLabel <- "reference"
  noDigits <- 2
  conf.high <- conf.low <- estimate <- NULL
  
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  
  
  
  allTerms <- lapply(seq_along(terms), function(i){
    var <- base::names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data),
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data),
                 pos = seq_along(vars))
    }
  })
  
  
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[,1:2], 1, paste0, collapse="")
  
  
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds,])[,c("var", "level","estimate", "conf.low", "conf.high", "pos")]
  toShow[is.na(toShow)] <- 0
  toShow <- toShow %>% as_tibble() %>% mutate_at(c('estimate','conf.low','conf.high'),exp) %>%
    mutate(name=var)
  #lookup[as.character(var)]) %>%
  #  mutate(name=ifelse(is.na(name),as.character(var),name))
  
  levels <- unique(toShow$name)
  toShow <- toShow %>% mutate(name=factor(name,levels=levels))
  
  toShow <- toShow %>% mutate(level = ifelse(pos==1,paste0('<b>',as.character(level),' (ref) </b>'),as.character(level)))
  return (toShow);
}

forestplot <- function (model,data,toShow=NULL,facet=T,ymin=0.,ymax=100.) {

  if(is.null(toShow)){
    toShow <- get_forest_df(model,data)
  }
  
  print (toShow)
  
  p<- toShow %>% 
    ggplot(aes(x=reorder(level,-pos)))+
    geom_pointrange(aes(y=estimate, ymin=conf.low, ymax=conf.high, fill='Adjusted'), size=1,width=.4,alpha=1.0,shape=21,
                    position=position_dodge(.9)) +
    scale_fill_manual(values=c(phs_colours("phs-blue"), phs_colours("phs-green"))) +
    ylim(ymin, ymax) +
    labs(title='',x='',y='Hazard Ratio') +
    geom_hline(yintercept=c(1), linetype="dashed") +
    geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
    #ylim(0., 10) +
    scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
    coord_flip()+#ylim=c(xmin,xmax)) +
    theme_classic() + guides(fill="none")
  
  if(!facet){
    return (p)
  }
  
  return (p + facet_grid(name ~ .,scale='free', switch='both') + 
            theme(
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.placement = 'outside',
              axis.text.y = element_markdown(angle=0),
              strip.text.y.left = element_text(angle = 0)));
}

haztab <- function(res.cox){
  tab <- exp(coef(res.cox)) %>% as.data.frame %>% merge(exp(confint(res.cox))  %>% as.data.frame, by=0)
  tab <- tab %>% as_tibble %>% filter(grepl('outcome',Row.names)) %>% 
    mutate(Row.names = gsub('outcome_igg','',Row.names)) %>% 
    rename(`IgG Quantile`='Row.names','Hazard Ratio'=".") %>% 
    arrange(`Hazard Ratio`) %>% kbl() %>% kable_styling()
  
  return (tab)
}


get_base_haz <- function(res.cox,ndf=4){
  H0 <- basehaz(res.cox,centered = F)
  #H0 %>% ggplot(aes(x=time,y=hazard)) + geom_point()  + geom_smooth(method = 'gam',formula = y~ns(x,df=100) ) +
  #  labs(x='Time (days)',y='Cumulative Baseline Hazard (H0(t))') +
  #  theme_minimal()
  H0.model <- gam(hazard ~ ns(time,df=4), data=H0)
  tseq <- seq(from=0, to=600, length=600)
  H0.pred <- predict(H0.model, newdata = data.frame(time=tseq), se=TRUE)
  H0.pred.df <- H0.pred %>% as_tibble %>% mutate(x=row_number(),y=fit,yerr=se.fit) %>% select(x,y,yerr)
  
  H0.pred.df2 <- H0.pred.df %>% mutate(diff=y-lag(y),differr=yerr-lag(yerr)) %>% filter(!is.na(y)) %>% 
    mutate(date=as.Date(x, origin = '2020-12-28'))
  
  return (H0.pred.df2);
  
}
#plot_base_haz <- function(df){
#  df %>% ggplot(aes(x=date,y=diff)) + 
#    #geom_ribbon(color='red',fill='red',alpha=0.4,aes(ymin=diff-differr,ymax=diff+differr)) +
#    geom_line(color='red',fill='red',alpha=0.9,size=2) +
#    #geom_vline(xintercept=c(as.Date('2021-07-15'),as.Date('2021-11-20')),linetype='dotted') +
#    labs(x='Date',y='Basesline Hazard') +
#    theme_minimal() -> p
#  return (p)
#}



#haz <- exp(diff(H0[, 'hazard'])*diff(H0[, 'time']))
#H <- H0$hazard * exp(sum(res.cox$coefficients))
#S <- exp(-H) 


df_ana_pc <- load_pc()
df_ana_bd <- load_bd() %>% mutate(outcome=ifelse(is.na(outcome),0,outcome))

min_date <- min(min(df_ana_pc$Sampledate_iso),min(df_ana_bd$Sampledate_iso))

df_smr01_covid <- df_smr01 %>% filter(ADMISSION_DATE > min_date) %>% filter(ADMISSION_COVID==1) %>% select(EAVE_LINKNO,ADMISSION_DATE) %>% group_by(EAVE_LINKNO,ADMISSION_DATE) %>% 
  filter(n()==1) %>% ungroup %>% arrange(EAVE_LINKNO)


df_deaths_covid <- df_deaths %>%  filter(UNDERLYING_CAUSE_OF_DEATH=='U071' | UNDERLYING_CAUSE_OF_DEATH =='U072') %>%
  mutate(death_date=NRS.Date.Death) %>% rename(ADMISSION_DATE=NRS.Date.Death) %>% 
  mutate(ADMISSION_DATE=as.POSIXct(ADMISSION_DATE))

df_deaths_covid

df_deaths_any <- df_deaths %>% select(EAVE_LINKNO,NRS.Date.Death)

df_severe_1 <- df_smr01_covid %>% rbind(df_deaths_covid %>% select(EAVE_LINKNO,ADMISSION_DATE)) %>% 
  arrange(EAVE_LINKNO) 

df_infect <- df_cdw %>% rename(ADMISSION_DATE=CdwDate) %>% select(-death28)



#levels <- c(0.01,0.05,0.1,0.25,0.5,0.75)
#labels <- c('75%','50%','25%','10%','5%','1%')
#levels <- c(0,0.02,0.05,0.1,0.25,0.5,0.75)
#levels <- c(0,0.1,0.2,0.3,0.4,0.6)
#labels <- c('10%','20%','30%','40%','60%')
#m_igg_stage <- quantile(modelC$IgG,levels)
#m_igg_stage

#max((modelC %>% filter(outcome==1))$IgG)
#mean((modelC %>% filter(outcome==0 & IgG<2000))$IgG)

#modelC %>% ggplot(aes(x=IgG,fill=outcome_igg)) + 
#  geom_histogram() +
#  scale_x_log10()

max((df_ana_bd %>% filter(outcome==1 & !is.na(IgG)))$IgG)

df_ana_bd %>% filter(outcome==1) %>% select(outcome,IgG)

df_ana_bd %>% ggplot(aes(x=IgG,fill=as.factor(outcome))) + 
  geom_histogram() +
  scale_x_log10()

create_modelC <- function(df_ana,f,df_outcome,split=T){
  

  modelC <- df_ana %>%  contrast_code() %>% select(EAVE_LINKNO,age,ageYear,Sex,IgG,insufficient_response,shielding,
                                                    ch_resident,immuno,simd2020v2_sc_quintile,
                                                    Sampledate_iso,Q_BMI,n_risk_gps,immuno,
                                                    d1_datetime,d2_datetime,d3_datetime,d4_datetime,
                                                    cat,stage,contains('Q_DIAG'))

  labels <- labels(f)

  calculate_igg_quartiles <- function(x){
    return (labels(f[sum(f<=x)])[1])
  }

  modelC$outcome_igg <- mapply(calculate_igg_quartiles, modelC$IgG)
  modelC <- modelC %>% mutate(outcome_igg = factor(outcome_igg,levels=labels)) 
 

  model <- modelC %>% left_join(df_outcome) %>% left_join(df_deaths_any)

  max_date <- max((model %>% filter(!is.na(ADMISSION_DATE)))$ADMISSION_DATE)
  min_date <- min(model$Sampledate_iso)
  
  print(max_date)
  print(min_date)

  #set the max date
  # - end of the experiment
  # - or the date of death (for any reason)
  model <- model %>% 
         mutate(max_date = case_when(is.na(NRS.Date.Death) ~ as.Date(max_date),
                                      TRUE ~ as.Date(NRS.Date.Death)))

  model <- model %>%
           mutate(days_since_sample=as.integer(as.Date(ADMISSION_DATE) - as.Date(Sampledate_iso),units='days')) %>%
           mutate(days_since_sample = ifelse(days_since_sample<0,NA,days_since_sample)) %>%
           group_by(EAVE_LINKNO,Sampledate_iso) %>%
           filter(row_number()==1) %>%
           ungroup %>%
           mutate(outcome_hosp = ifelse(is.na(days_since_sample),0,1))%>%
           mutate(days_since_measurement = as.numeric(as.Date(ADMISSION_DATE) - as.Date(Sampledate_iso),units='days' )) %>%
           mutate_at(c("ADMISSION_DATE","Sampledate_iso"),as.Date) %>%
           mutate_at(c("insufficient_response","ch_resident","shielding"),~as.factor(ifelse(.==1,'Yes','No'))) %>%
           mutate(
             t1=as.numeric(Sampledate_iso - as.Date(min_date),units='days'),
             t2=ifelse(outcome_hosp==1,
                   t1 + days_since_sample,
                   as.numeric(max_date - as.Date(min_date),units='days' )+1)) %>%
           mutate(tv2 = t1 + as.numeric(d2_datetime - as.Date(Sampledate_iso),units='days'),
                  tv3 = t1 + as.numeric(d3_datetime - as.Date(Sampledate_iso),units='days'),
                  tv4 = t1 + as.numeric(d4_datetime - as.Date(Sampledate_iso),units='days')
                 ) %>%
           mutate(
             days_since_last = case_when(
                                        !is.na(tv4) & tv4>t1 & tv4<t2 ~ t2-tv4,
                                        !is.na(tv3) & tv3>t1 & tv3<t2 ~ t2-tv3,
                                        !is.na(tv2) & tv2>t1 & tv2<t2 ~ t2-tv2,
                                        TRUE ~ 0),
             days_since_last = as.factor(case_when( days_since_last==0 ~ 'None',
                                          days_since_last<70 ~ '<150',
                                          #days_since_last<200 ~ '<200',
                                          TRUE ~ '>150')),
             additional = as.numeric(!is.na(tv2) & tv2>t1 & tv2<t2) + 
                               as.numeric(!is.na(tv3) & tv3>t1 & tv3<t2) + 
                               as.numeric(!is.na(tv4) & tv4>t1 & tv4<t2) )
          #%>%
          #mutate(cat2=fct_collapse(cat, 
          #                 "Pf-Pf" = c("Pf-Pf","mixed 2 doses"),
          #                 "3 doses" = c("mixed 3 doses","AZ-AZ-AZ","Pf-Pf-Pf")
          #)) %>%
          #mutate(cat2 = fct_relevel(cat2,'AZ-AZ')) %>%
          #mutate(stage2 = fct_collapse(stage,"3+"=c("3","4+"))) %>%
          #mutate(ageYear2 = fct_collapse(ageYear,"0-39"=c("0-19","20-39")))
          # %>
          #mutate(additional = fct_collapse(additional,"2+"=c("2","3")))
  
  if (split == F){
    return (model);
  }

  temp <- model %>% mutate(id=row_number()) %>% select(id,everything()) %>%
         pivot_longer(c("tv2","tv3","tv4"), names_to = "names", values_to = "days") %>%
         mutate(tt1=t1,tt2=t2,days=ifelse(days>t1 & days<t2,days,NA)) %>%
         group_by(id) %>%
         slice(c(1:n()), n()) %>%
         filter(row_number()==n() | !is.na(days)) %>%
         mutate(t1=ifelse(is.na(lag(days)),t1,lag(days)),t2=days) %>%
         #mutate(t1=tt1,t2=days) %>%
         mutate(
           additional = row_number()-1,
           outcome_hosp = ifelse(row_number()==n(),outcome_hosp,0),
           #ttt1 = ifelse(row_number()==n(),t2,t1),
           t2 = ifelse(row_number()==n(),tt2,t2)) %>% 
         select(-tt1,-tt2,-days,-names) %>%
         ungroup %>%
         mutate(additional=fct_collapse(as.factor(additional),"2+"=c("2","3")))%>% 
         mutate(outcome_igg = fct_relevel(outcome_igg,'Average'))

  return (temp);
}

f.pc <- base::c('Undetectable'=0.,'Insufficient'=4.82,'Below Average'=33.3,'Average'=230.,'Above Average'=2000)
f.bd <- base::c('Undetectable'=-1.,'Insufficient'=0.11,'Below Average'=0.7,'Average'=5,'Above Average'=8.5)


pc.modelC.infect <- create_modelC(df_ana_pc,f.pc,df_infect,split=F) %>%
                     left_join(df_demo %>% as_tibble %>% select(EAVE_LINKNO,ur6_2016_name) %>% 
                               filter(!is.na(ur6_2016_name)) %>%
                               mutate(ur6_2016_name=as.factor(ur6_2016_name)) %>% droplevels)


bd.modelC.infect <- create_modelC(df_ana_bd,f.bd,df_infect) %>%
                     left_join(df_demo %>% as_tibble %>% select(EAVE_LINKNO,ur6_2016_name) %>% 
                               filter(!is.na(ur6_2016_name)) %>%
                               mutate(ur6_2016_name=as.factor(ur6_2016_name)) %>% droplevels)

#pc.modelC.infect %>% filter(cat!='AZ-AZ-AZ' & cat!='mixed 2 doses') %>% group_by(outcome_hosp) %>% summarise(n=n())
#pc.modelC.infect %>% filter(cat!='AZ-AZ-AZ' & cat!='mixed 2 doses') %>% group_by(outcome_hosp,outcome) %>% summarise(n=n())
                     
                     
                     

#bd.modelC.infect %>% group_by(outcome_igg) %>% summarise(n=n())
#bd.modelC.infect %>% ggplot(aes(x=IgG,fill=as.factor(outcome_igg))) + geom_histogram()

cox.pc.modelC.infect <- with(pc.modelC.infect,
                             coxph(Surv(t1,t2, outcome_hosp) ~  outcome + n_risk_gps + ageYear + Sex + shielding + ch_resident + Q_BMI + simd2020v2_sc_quintile + ur6_2016_name + additional ))
    
forestplot(cox.pc.modelC.infect,data=pc.modelC.infect)

                         
cox.bd.modelC.infect <- with(bd.modelC.infect,
                             coxph(Surv(t1,t2, outcome_hosp) ~  outcome + n_risk_gps + ageYear + Sex + shielding  + Q_BMI + simd2020v2_sc_quintile + ur6_2016_name + additional ))



cox.df.pc.modelC.infect <- get_forest_df(cox.pc.modelC.infect,pc.modelC.infect)
cox.df.bd.modelC.infect <- get_forest_df(cox.bd.modelC.infect,bd.modelC.infect)



cox.df.modelC.infect <- cox.df.pc.modelC.infect %>% mutate(dataset='Primary Care') %>% 
                           rbind(cox.df.bd.modelC.infect %>% mutate(dataset='Blood Donors')) %>%
                                 mutate(dataset=fct_rev(as.factor(dataset)))


p.modelC.infect.HRs <- cox.df.modelC.infect %>% 
  ggplot(aes(x=reorder(level,-pos)))+
  geom_pointrange(aes(y=estimate, ymin=conf.low, ymax=conf.high, fill=dataset), size=1,width=.4,alpha=1.0,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-blue"), phs_colours("phs-rust"))) +
  #ylim(0., 16.0) +
  labs(title='',x='',y='Hazard Ratio') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.333,3), linetype="dotted") +
  #ylim(0., 10) +
  scale_y_log10(breaks=c(0.33,1,3)) +
  coord_flip()+#ylim=c(xmin,xmax)) +
  theme_classic() + guides(fill="none")+ 
  facet_grid(name ~ dataset,scale='free', switch='y') + theme(
  panel.background = element_rect(fill = NA, color = "black"),
  strip.background = element_blank(),
  strip.placement = 'outside',
  axis.text.y = element_markdown(angle=0),
  strip.text.y.left = element_text(angle = 0))


ggsave(paste0(folder,"modelC_infect_HRs.pdf"), p.modelC.infect.HRs, width=10, height=12, dpi=300, units="in")



df <-  get_base_haz(cox.pc.modelC.infect) %>% mutate(dataset='Primary Care')%>% 
  rbind(get_base_haz(cox.bd.modelC.infect) %>% mutate(dataset='Blood Donors')) %>%
  mutate(dataset=as.factor(dataset))

df %>% ggplot(aes(x=date,y=diff,color=dataset)) +
  geom_line(alpha=0.9,size=2) +
  labs(x='Date',y='Basesline Hazard',color='Dataset') +
  #scale_colour_discrete_phs(palette='main') +
  scale_colour_manual(values = phs_colours(c("phs-rust","phs-blue"))) +
  theme_classic() -> p.modelC.infect.basehaz

ggsave(paste0(folder,"modelC_infect_basehaz.pdf"), p.modelC.infect.basehaz, width=7, height=5, dpi=300, units="in")



cox.pc.modelC.infect.full <- with(pc.modelC.infect,
                             coxph(Surv(t1,t2, outcome_hosp) ~  outcome_igg + n_risk_gps + ageYear + Sex + shielding + ch_resident + Q_BMI + simd2020v2_sc_quintile + ur6_2016_name + additional ))

cox.bd.modelC.infect.full <- with(bd.modelC.infect,
                             coxph(Surv(t1,t2, outcome_hosp) ~  outcome_igg + n_risk_gps + ageYear + Sex + shielding  + Q_BMI + simd2020v2_sc_quintile + ur6_2016_name + additional ))



cox.df.pc.modelC.infect.full <- get_forest_df(cox.pc.modelC.infect.full,pc.modelC.infect)
cox.df.bd.modelC.infect.full <- get_forest_df(cox.bd.modelC.infect.full,bd.modelC.infect)


cox.df.modelC.infect.full <- cox.df.pc.modelC.infect.full %>% mutate(dataset='Primary Care') %>% 
  rbind(cox.df.bd.modelC.infect.full %>% mutate(dataset='Blood Donors')) %>%
  mutate(dataset=fct_rev(as.factor(dataset)))

p.modelC.infect.full.HRs <- cox.df.modelC.infect.full %>% 
  ggplot(aes(x=reorder(level,-pos)))+
  geom_pointrange(aes(y=estimate, ymin=conf.low, ymax=conf.high, fill=dataset), size=1,width=.4,alpha=1.0,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-blue"), phs_colours("phs-rust"))) +
  #ylim(0., 16.0) +
  labs(title='',x='',y='Hazard Ratio') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.333,3), linetype="dotted") +
  #ylim(0., 10) +
  scale_y_log10(breaks=c(0.33,1,3)) +
  coord_flip()+#ylim=c(xmin,xmax)) +
  theme_classic() + guides(fill="none")+ 
  facet_grid(name ~ dataset,scale='free', switch='y') + theme(
    panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.left = element_text(angle = 0))


ggsave(paste0(folder,"modelC_infect_full_HRs.pdf"), p.modelC.infect.full.HRs, width=10, height=12, dpi=300, units="in")






pc.modelC.hosp <- create_modelC(df_ana_pc,f.pc,df_severe_1) 
bd.modelC.hosp <- create_modelC(df_ana_bd,f.bd,df_severe_1) %>%
                  filter(ageYear!='0-19' & Q_BMI !='0-20' & Q_BMI !='25-30')

cox.pc.modelC.hosp <- with(pc.modelC.hosp,
                          coxph(Surv(t1,t2, outcome_hosp) ~ outcome + additional + n_risk_gps + ageYear + Sex + Q_BMI  ))



#View(df)



p.pc.modelC.hosp.HRs <- forestplot(cox.pc.modelC.hosp,data=pc.modelC.hosp)

cox.df.modelC.hosp <- get_forest_df(cox.pc.modelC.hosp,pc.modelC.hosp)

cox.df.modelC.hosp

ggsave(paste0(folder,"modelC_hosp_HRs.pdf"),p.pc.modelC.hosp.HRs, width=9, height=9, dpi=300, units="in")

cox.pc.modelC.hosp.full <- with(pc.modelC.hosp,
                           coxph(Surv(t1,t2, outcome_hosp) ~ outcome_igg + additional + n_risk_gps + ageYear + Sex + Q_BMI  ))

p.pc.modelC.hosp.full.HRs <- forestplot(cox.pc.modelC.hosp.full,data=pc.modelC.hosp)
ggsave(paste0(folder,"modelC_hosp_full_HRs.pdf"),p.pc.modelC.hosp.full.HRs, width=9, height=9, dpi=300, units="in")


qnames <- pc.modelC.hosp %>% get_qnames(var='outcome_hosp',n=1)
qnames

qnames2 <- c("Q_DIAG_BLOOD_CANCER","Q_DIAG_CHD","Q_DIAG_PULM_HYPER","Q_DIAG_COPD")
qnames2



formula1 <- paste0(paste0(qnames,collapse = ' + ')," - ",paste0(qnames2,collapse = ' - '))
formula1

#formula2 <- paste0("as.integer(n_risk_gps) - 1 - ",paste0(qnames,collapse = ' - ')) 
#formula2

temp <- pc.modelC.hosp %>%
  #mutate(n_risk_gpsOther=as.factor(as.integer(n_risk_gps)-1 - Q_DIAG_ASTHMA - Q_DIAG_DIABETES_1 ) ) %>%
  mutate(n_risk_gpsOther = !!parse_quo(formula1, env = caller_env())) %>%
  mutate(n_risk_gpsOther = as.factor(ifelse(n_risk_gpsOther>2,'3+',n_risk_gpsOther))) %>% 
  mutate_at(vars(contains("Q_DIAG")),~as.factor(ifelse(.==1,'Yes','No'))) 

temp %>% group_by(n_risk_gpsOther) %>% summarise(n=n())


formula <- paste0("Surv(t1,t2, outcome_hosp) ~ outcome_igg + additional  + ",#
                              paste0(qnames2,collapse=" + "),
                              "+ n_risk_gpsOther + ageYear + Sex + Q_BMI " )
formula

res.cox <- with(temp,coxph(as.formula(formula)))

p.pc.modelC.hosp.full.risks.HRs <- forestplot(res.cox,data=temp)
p.pc.modelC.hosp.full.risks.HRs


View(get_forest_df(res.cox,temp))

vars <- (attr(res.cox$terms,"variables") %>% as.list())
vars <- vars[3:length(vars)]
vars
df <- NULL
for (var in vars){
  cox.temp <- with(temp,coxph(as.formula(paste0('Surv(t1,t2, outcome_hosp) ~ ', var))))
  df.temp <- get_forest_df(cox.temp,temp)
  df <- df %>% rbind(df.temp)
}
View(df)


ggsave(paste0(folder,"modelC_hosp_full_risks_HRs.pdf"),p.pc.modelC.hosp.full.risks.HRs, width=9, height=12, dpi=300, units="in")


df <-  get_base_haz(cox.pc.modelC.hosp,ndf=4) %>% mutate(dataset='Primary Care')

df %>% ggplot(aes(x=date,y=diff,color=dataset)) +
  geom_line(alpha=0.9,size=2) +
  labs(x='Date',y='Basesline Hazard',color='Dataset') +
  #scale_colour_discrete_phs(palette='main') +
  scale_colour_manual(values = phs_colours(c("phs-rust","phs-blue"))) +
  theme_classic() -> p.modelC.hosp.basehaz

p.modelC.hosp.basehaz

ggsave(paste0(folder,"modelC_pc_hosp_basehaz.pdf"), p.modelC.hosp.basehaz, width=7, height=5, dpi=300, units="in")






temp <- bd.modelC.hosp %>%
        mutate(n_risk_gpsOther=as.factor(as.integer(n_risk_gps)-1 - Q_DIAG_ASTHMA - Q_DIAG_DIABETES_1 ) ) %>%
        mutate_at(vars(contains("Q_DIAG")),~as.factor(ifelse(.==1,'Yes','No'))) 

res.cox <- with(temp,
              coxph(Surv(t1,t2, outcome_hosp) ~ outcome + additional + n_risk_gpsOther + Q_DIAG_ASTHMA + Q_DIAG_DIABETES_1 + ageYear + Sex + Q_BMI  ))

forestplot(res.cox,data=temp)
 
View(temp %>% filter(outcome_hosp==1 & Q_DIAG_ASTHMA!='Yes'))

  
#cox.bd.modelC.hosp <- with(bd.modelC.hosp,
#                           coxph(Surv(t1,t2, outcome_hosp) ~ outcome + additional + n_risk_gps + ageYear + Sex + Q_BMI  ))
#View(bd.modelC.hosp %>% filter(outcome_hosp==1)) %>% group_by(n_risk_gps) %>% summarise(n=n())

forestplot(res.cox,data=temp)




#res.cox <- with(temp
#                ,
p <-forestplot(res.cox,data=temp,ymin=0.2,ymax=100.)
p 







haztab(res.cox)
summary(res.cox)
cox.zph(res.cox)
ggcoxzph(cox.zph(res.cox))
ggforest(res.cox,data=temp)




#model %>% select(t1,t2,tv2,tv3,tv4,Sampledate_iso,d2_datetime,additional)

  #mutate(
  #              event_date = Sampledate_iso + days(days_since_sample),
  #              additional = as.factor(as.numeric(days_v2>0) + as.numeric(days_v3>0))) %>%
  #         mutate(days_since_start = create_period_factor(Sampledate_iso,days_since_sample),
  #                days_since_start_3 = create_period_factor(event_date,days_since_sample),
  #                days_since_start_2 = as.factor(cut(as.Date(Sampledate_iso),
  #                              breaks=c(as.Date('2020-01-01'),as.Date('2021-05-01'),as.Date('2021-07-15'),
  #                                       as.Date('2021-11-01'),as.Date('2021-12-15'),as.Date('2022-02-01'),as.Date('2023-01-01')),
  #                              right=T,
  #                              labels=c('Winter 2020-2021','Spring 2021','Delta Period','Autumn 2021','Omicron Period','Spring 2022')))
  #         )


#model <- model %>% 

#%>%
         #mutate_at(c("days_since_start","days_since_start_2","days_since_start_3"), ~fct_relevel(.,'Delta Period')) 
         #mutate_at('days_since_start',~fct_relevel(.,'2021-08-02'))



#model %>% ggplot(aes(x=Sampledate_iso)) + geom_histogram()
#model %>% ggplot(aes(x=event_date)) + geom_histogram()


#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps + ageYear + Sex + Q_BMI + additional))
#prior_hospitalisation

temp <- model %>% mutate(ageYear=as.factor(cut(age,breaks=c(0,20,30,40,50,60,70,80,1000),right=T,
                                         labels=c('0-20','20-30','30-40','40-50','50-60','60-70','70-80','80')))) %>%
                  mutate(ageYear=fct_relevel(ageYear,'40-50'))##

#hist(model$age)

#temp <- model %>% mutate(ageYear=as.factor(cut(age,breaks=c(0,20,30,40,50,60,1000),right=T,
#                                               labels=c('0-20','20-30','30-40','40-50','50-60','60+')))) %>%
#  mutate(ageYear=fct_relevel(ageYear,'40-50'))##




#qnames2 <- c("Q_DIAG_BLOOD_CANCER","Q_DIAG_COPD","Q_DIAG_CHD")
#model <- model %>% 
#          mutate(n_risk_gpsOther=as.integer(n_risk_gps)-1 - Q_DIAG_BLOOD_CANCER  - Q_DIAG_COPD - Q_DIAG_CHD ) %>%
#          mutate(temp=paste(ifelse(Q_DIAG_BLOOD_CANCER==1,'BC',''),ifelse(Q_DIAG_CKD_LEVEL==1,'CKD',''),
#                             ifelse(Q_DIAG_COPD==1,'COPD',''),ifelse(Q_DIAG_CHD==1,'CHD',''),sep='')) %>%
#          mutate(temp = as.factor(temp), n_risk_gpsOther = as.factor(ifelse(n_risk_gpsOther>=3,'3+',n_risk_gpsOther)))
model <- model %>% mutate(cat2=fct_collapse(cat, 
                                  "Pf-Pf" = c("Pf-Pf","mixed 2 doses"),
                                  "3 doses" = c("mixed 3 doses","AZ-AZ-AZ","Pf-Pf-Pf")
                                  )) %>%
                   mutate(cat2 = fct_relevel(cat2,'AZ-AZ'))

#model %>% filter(outcome_hosp ==1 ) %>% group_by(cat2) %>% summarise(n=n())

#qnames <- get_qnames(model %>% filter(outcome_hosp==1))
#qnames


#model %>% head %>% select(qnames) %>%  pivot_longer(cols = qnames, names_to = 'cond', values_to = 'Values')



#res.cox <- with(temp,
#                coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps + ageYear + shielding + Sex + Q_BMI + simd2020v2_sc_quintile + ur6_2016_name + additional + days_since_start))

#bd.res.cox <- res.cox
#pc.res.cox <- res.cox
#forestplot(res.cox,temp)#+ scale_y_log10(lim=c(0.2,5), breaks=c(0.33,1,3)) 


#pc.df <- get_forest_df(res.cox,temp)
#pc.df

bd.df <- get_forest_df(res.cox,temp)


model %>% group_by(outcome_hosp) %>% summarise(n=n())
View(model %>% filter(outcome_hosp ==1))


df <- pc.df %>% mutate(dataset='Primary Care') %>% rbind(bd.df %>% mutate(dataset='Blood Donors')) %>%
      mutate(dataset=fct_rev(as.factor(dataset)))

df %>% group_by(dataset) %>% summarise(n=n())
             
p <- df %>% 
  ggplot(aes(x=reorder(level,-pos)))+
  geom_pointrange(aes(y=estimate, ymin=conf.low, ymax=conf.high, fill=dataset), size=1,width=.4,alpha=1.0,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-blue"), phs_colours("phs-rust"))) +
  ylim(0., 16.0) +
  labs(title='',x='',y='Hazard Ratio') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  #ylim(0., 10) +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  coord_flip()+#ylim=c(xmin,xmax)) +
  theme_classic() + guides(fill="none")
p


p + facet_grid(name ~ dataset,scale='free', switch='y') + theme(
    panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.left = element_text(angle = 0))

#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps_2 + Q_DIAG_BLOOD_CANCER + Q_DIAG_CKD_LEVEL + ageYear + Sex + Q_BMI + additional + days_since_start))

#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps + ageYear + Sex + Q_BMI + additional + days_since_start))



#formula <- paste0('Surv(time, outcome_hosp) ~  outcome_igg +', paste0(qnames,collapse = " + "),' + ageYear + Sex + Q_BMI + additional + days_since_start')
#formula
#res.cox <- with(model %>% mutate_at(qnames,as.factor),coxph(as.formula(formula)))

#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps_2 + temp + ageYear + Sex + Q_BMI + additional + days_since_start))

#formula <- paste0('Surv(time, outcome_hosp) ~  outcome_igg  + ', paste0(qnames2,collapse = " + "),' + n_risk_gpsOther + ageYear + Sex + Q_BMI + additional + days_since_start')
#formula
#res.cox <- with(model %>% mutate_at(qnames2,as.factor),coxph(as.formula(formula)))



#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + n_risk_gps + days_since_start + additional))
#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg + Q_DIAG_BLOOD_CANCER + Q_DIAG_CKD_LEVEL + Q_DIAG_CHD +  Q_DIAG_DIABETES_1 + days_since_start))

#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg  + additional))

#cat2
#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg  + n_risk_gps + ageYear + Sex + Q_BMI + additional + days_since_start))
#forestplot(res.cox,model)


#res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg  + n_risk_gps + ageYear + Sex + Q_BMI + additional))
#forestplot(res.cox,model)



#ggcoxdiagnostics(res.cox, type = "dfbeta",
#                 linear.predictions = FALSE, ggtheme = theme_bw())



ggcoxdiagnostics(res.cox1,
                 type = "deviance",
                 ox.scale = "linear.predictions")

ggcoxdiagnostics(res.cox1,
                 type = "schoenfeld",
                 ox.scale = "Time")


df1 <- get_forest_df(res.cox1,temp)
df2 <- get_forest_df(res.cox2,temp)

df <- df1 %>% mutate(dataset='Surv(t1,t2,event)') %>% rbind(df2 %>% mutate(dataset='Surv(t,event)'))


df %>% 
  ggplot(aes(x=reorder(level,-pos)))+
  geom_pointrange(aes(y=estimate, ymin=conf.low, ymax=conf.high, fill=dataset), size=1,width=.4,alpha=1.0,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-blue"), phs_colours("phs-rust"),phs_colours("phs-teal"))) +
  ylim(0., 16.0) +
  labs(title='',x='',y='Hazard Ratio',fill='Fit Type') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  #ylim(0., 10) +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  coord_flip()+#ylim=c(xmin,xmax)) +
  theme_classic() + guides() +
  facet_grid(name ~ .,scale='free', switch='y') + theme(
    panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.left = element_text(angle = 0))



#model <- model %>% mutate(days_since_start_3 = fct_explicit_na(days_since_start_3, "END")) 
fit <- survfit(Surv(time, outcome_hosp) ~  outcome + days_since_start ,data=temp)
ggsurvplot(fit,
           #newdata=new_df,
           facet.by='days_since_start',
           palette=phs_colours(c("phs-blue", "phs-purple","phs-magenta","phs-graphite","phs-teal","phs-rust")),
           conf.int=T,
           xlab = "Days since Serology Measurement", 
           legend.title = '',
           #legend.labs =
           #    c('Insufficient IgG','Sufficient IgG'),
              #c("IgG above assay max",'Average IgG','Bottom 25% IgG','Bottom 10% IgG','Insufficient IgG','Fail-to-Mount any IgG'),
           risk.table = F,
           #xlim=c(0,460),
           ylim=c(0.970,1)) 




#temp %>% select(t1,t2,time,outcome_hosp) %>% filter(outcome_hosp==1) 





plot(survfit(res.cox1))
plot(survfit(res.cox2))
plot(survfit(res.cox3))


res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg  + n_risk_gps + ageYear + Sex + Q_BMI + additional + days_since_start))
forestplot(res.cox,model)

res.cox <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome_igg  + n_risk_gps + ageYear + Sex + Q_BMI + additional + days_since_start_3))
forestplot(res.cox,model)




pc.res.cox <- res.cox



p <- ggforest(res.cox,data=model,fontsize = 0.7)
p




res.cox

p <- ggforest(res.cox,data=model,fontsize = 0.7)
p


hist(model$IgG)

ggplot(model,aes(x=IgG,fill=outcome_igg)) +  
       geom_histogram(bins=25) +
       scale_fill_discrete_phs(palette = 'all') +
       scale_x_log10() +
       labs(x='IgG [BAU]',y='Number of Measurements',fill='IgG Outcome')

ggplot(model, aes(x=IgG,fill=outcome_igg)) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '75%'), aes(fill='75%')) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '50%'), aes(fill='50%')) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '25%'), aes(fill='25%')) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '10%'), aes(fill='10%')) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '5%'), aes(fill='5%')) +
  geom_histogram(position = "stack", bins=30, data=subset(model,outcome_igg == '1%'), aes(fill='1%')) +
  labs(title='All Measurements',x="Antibody Levels [IgG]", y="Number of Samples", fill="Result") +
  scale_fill_discrete_phs(palette = 'all') +
  theme(aspect.ratio = 0.5)  +
  scale_y_log10() +
  scale_x_log10()


colnames(model)


#cfit <- with(model,coxph(Surv(time, outcome_hosp) ~  outcome))



new_data = data.frame(
  strata = c("stage1,outcome=0","stage1,outcome=1","stage2,outcome=0","stage2,outcome=1"),
  stage = c(1,1,'2+','2+'),
  outcome = c(0,1,0,1),
  outcome_hosp = c(1,1,1,1)
) %>% as_tibble

#fit <- survfit(Surv(time, outcome_hosp) ~  outcome + stage,data=model)
fit <- survfit(Surv(time, outcome_hosp) ~  outcome_igg  ,data=model)
#ggcompetingrisks(fit)




ggsurvplot(fit,
           newdata=new_data,
           palette=phs_colours(c("phs-blue", "phs-purple","phs-magenta","phs-graphite","phs-teal","phs-rust")),
           conf.int=T,
           xlab = "Days since Serology Measurement", 
           legend.title = '',
           #legend.labs =
          #   c("IgG above assay max",'Average IgG','Bottom 25% IgG','Bottom 10% IgG','Insufficient IgG','Fail-to-Mount any IgG'),
           risk.table = TRUE,
           xlim=c(0,460),
           ylim=c(0.970,1)) 

ggsurvtable(fit,data=model)

#ggsurvplot(survfit(cfit,data=new_data),data=new_data) %>% scale_y_log10()

#tidy(cfit)

plot(survfit(cfit), 
     xlab = "days since sample",
     ylab = "Proportion not hospitalised")


#newdat <- with(model, 
#               data.frame(
#                 outcome = c(0, 1),
#                 days_since_sample = c(50,50)
#               )
#)


data <-lung %>% as_tibble %>% select(time,status,sex)
data
fit <- survfit(Surv(time, status) ~ sex, data = data)
ggsurvplot(fit, data = data)


plot(survfit(cfit,newdata=newdat), 
     xlab = "days since first measurement",
     ylab = "Proportion hospitalised")





colnames(modelD)
modelD %>% group_by(outcome,stage,outcome_hospitalisation) %>% summarise(n=n())    

df_severe <- 

nrow(df_ana)



nrow(modelD)
modelD %>% filter(outcome_h_covid==1) %>%
  ggplot(aes(x=days_hospitalisation,fill='Stage3')) +
           geom_histogram(bins=20) +
           scale_fill_discrete_phs(palette='main')




igg_stages <- list()
for (dose in unique(levels(modelD$stage))){
   levels <- quantile((modelD %>% filter(stage==dose))$IgG,c(0.01,0.05,0.1,0.25,0.5)) 
   igg_stages[[dose]]= levels
}


levels <- c(0.01,0.05,0.1,0.25,0.5,0.75)
m_igg_stage <- quantile(modelD$IgG,levels) 

calculate_igg_quartiles <- function(x,s){
  #f <- igg_stages[[s]]
  #f <- igg_stages[['1']]
  f <- m_igg_stage 
  i <- sum(f<=x)
  retval <- names(f[i])[1]
  return (retval);#ifelse(is.na(retval),'50%',retval));
}


modelD$outcome_igg <- mapply(calculate_igg_quartiles, modelD$IgG,modelD$stage)

levels <- c('75%','50%','25%','10%','5%','1%')
modelD <- modelD %>% mutate(outcome_igg = factor(outcome_igg,levels=levels)) %>%
                     filter(!is.na(outcome_igg))

modelD <- modelD %>% mutate(outcome_igg = fct_rev(outcome_igg))
modelD <- modelD %>% mutate(outcome_death = ifelse(is.na(death_date),0,
                                                   ifelse(death_date>Sampledate_iso,1,0)))



table <- modelD %>% filter(outcome_h_covid==1) %>%  group_by(outcome_igg,outcome_death) %>% summarise(n=n())
table

modelD %>% group_by(Q_BMI) %>% summarise(n=n())

table <- modelD %>% filter(Q_BMI != 'Unknown') %>% group_by(outcome_h_covid,outcome_igg) %>% summarise(n=n()) %>%
  group_by(outcome_igg) %>% summarise(outcome_covid=outcome_h_covid,nall=sum(n),n=n,p=100*n/sum(n)) %>%
  mutate(n=paste0(n," (",round(p,2),"%)")) %>%
  #mutate(outcome_igg=poutcome_igg)]) %>%
  filter(outcome_covid==1) %>%
  select(outcome_igg,nall,n)
colnames(table) <- c("IgG Level Quantile","Total",paste0("Hospitalised Due to COVID-19 After Serology Measurement"))
table

View(table)




  
formula <- paste0('outcome_h_for_covid ~ ageYear + Sex  + Q_BMI + ',paste0(qnames,collapse = ' + '),' + outcome_igg ')
formula <- paste0('outcome_covid ~ ageYear + Sex + Q_BMI + n_risk_gps + outcome_igg ')
formula <- paste0('outcome_h_covid ~  ageYear + Sex + Q_BMI + n_risk_gps + outcome +  shielding  + immuno_supp + 
                      days_since_first_measurement + days_since_last_vac')
#formula <- paste0('outcome_h ~ ageYear + Sex + Q_BMI + n_risk_gps + outcome_igg ')
#formula <- paste0('outcome_d ~ ageYear + Sex   + Q_BMI + n_risk_gps + outcome + cat ')
#formula <- paste0('outcome_covid ~ ageYear + Sex   + Q_BMI + n_risk_gps + outcome2 + cat')
#formula <- paste0('outcome_covid ~ ageYear + Sex   + Q_BMI + ',paste0(qnames,collapse = ' + '),' ')
#formula <- paste0('outcome_covid ~ ageYear + Sex + Q_BMI + outcome')
formula

#formula <- "outcome_covid ~ ageYear + Sex  + Q_BMI  + Q_DIAG_BLOOD_CANCER*outcome + cat "

#modelD %>% group_by(outcome_covid,outcome) %>% summarise(n=n())

#modelD %>% group_by(cat) %>% summarise(n=n())

nrow(modelD)
tab <- perform_glm(modelD,formula ,do_unadjusted = F)

plot_simple(tab)


formula <- paste0('outcome_h_covid ~  s(age) + Sex + Q_BMI + n_risk_gps + outcome_igg +  shielding  + immuno_supp + 
                      s(days_since_first_measurement) + s(days_since_last_vac)')

gamFit <- perform_gam(modelD,as.formula(formula))

or <- get_or_from_gam(gamFit)
pc <- plot_simple(or)
pc

formula <- paste0('outcome_hospitalisation ~  s(age) + Sex + Q_BMI + n_risk_gps + outcome_igg +  shielding  + immuno_supp + 
                      s(days_since_first_measurement) + s(days_since_last_vac)')

gamFit <- perform_gam(modelD,as.formula(formula))

or <- get_or_from_gam(gamFit)
ph<- plot_simple(or)
ph

grid.arrange(pc,ph,nrow=1)





tab <- perform_glm(modelD,formula ,do_unadjusted = T)
plot_fit(tab,'hospt')




table <- modelD %>% group_by(stage,outcome_covid,outcome_igg) %>% summarise(n=n()) %>%
           group_by(stage,outcome_igg) %>% summarise(outcome_covid=outcome_covid,nall=sum(n),n=n,p=100*n/sum(n)) %>%
           mutate(n=paste0(n," (",round(p,2),"%)")) %>%
           #mutate(outcome_igg=poutcome_igg)]) %>%
           filter(outcome_covid==1) %>%
           mutate(stage = as.character(stage)) %>% arrange(stage) %>% 
           select(stage,outcome_igg,nall,n)
colnames(table) <- c("Vaccine Stage","Serology Outcome","Total","Hospitalised Due to COVID-19")
table
View(table)





modelD %>% filter(stage==3 & outcome==1 & outcome_covid==1) %>% select(days_hospitalisation)

tab <- perform_glm(modelD,formula ,do_unadjusted = F)
tab

plot_simple(tab)#,'hospt')



nrow(modelD)
modelF <- modelD %>% filter(!is.na(death_date))
nrow(modelF)

View(modelF)

modelF %>% group_by(outcome) %>% summarise(n=n())


modelD %>% filter(stage==3) %>% group_by(outcome_d,outcome) %>% summarise(n=n()) %>% group_by(outcome_d) %>% mutate(p=100*n/sum(n)) %>% ungroup

