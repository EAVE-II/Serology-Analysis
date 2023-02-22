
library(ggplot2)
library(nlstools)
library(mgcv)
library(phsstyles)
library(facetscales)
library(eavehelpers)
library(lme4)
library(broom.mixed)
library("stringr")

eave.data <- eavehelpers::load_recommended_datasets() 
labels(eave.data)


df_pc <- eave.data$serology_primary_care %>% 
  left_join(eave.data$qcovid1) %>%
  left_join(eave.data$c19vaccine) %>%
  left_join(eave.data$demographics)


df_bd <- eave.data$serology_blood_donors %>% 
  left_join(eave.data$qcovid1) %>%
  left_join(eave.data$c19vaccine) %>%
  left_join(eave.data$demographics)

saveRDS(df_pc,"/conf/EAVE/GPanalysis/data/temp/serology_Feb2023.rds")


get_analysis_df <- function(df){
  df %>% filter(Sampledate_iso > d1_datetime) %>%
    mutate_at(c('Sampledate_iso','d1_datetime','d2_datetime','d3_datetime','d4_datetime'),as.Date) %>% 
    mutate(time_since_vaccination = case_when(
      !is.na(d4_datetime) & Sampledate_iso > d4_datetime ~ as.numeric(Sampledate_iso - d4_datetime,units='days'),
      !is.na(d3_datetime) & Sampledate_iso > d3_datetime ~ as.numeric(Sampledate_iso - d3_datetime,units='days'),
      !is.na(d2_datetime) & Sampledate_iso > d2_datetime ~ as.numeric(Sampledate_iso - d2_datetime,units='days'),
      TRUE ~ as.numeric(Sampledate_iso - d1_datetime,units='days'))
    ) %>%
    mutate(dose = case_when(
      !is.na(d4_datetime) & Sampledate_iso > d4_datetime ~ 4,
      !is.na(d3_datetime) & Sampledate_iso > d3_datetime ~ 3,
      !is.na(d2_datetime) & Sampledate_iso > d2_datetime ~ 2,
      TRUE ~ 1))%>%
    mutate(last_product = case_when(
      !is.na(d4_datetime) & Sampledate_iso > d4_datetime ~ d4_product,
      !is.na(d3_datetime) & Sampledate_iso > d3_datetime ~ d3_product,
      !is.na(d2_datetime) & Sampledate_iso > d2_datetime ~ d2_product,
      TRUE ~ d1_product)
    ) %>% 
    mutate_at(vars(contains('_product')), ~ 
                case_when(
                  grepl("AstraZeneca",.) ~ 'AZ',
                  grepl("Pfizer",.) ~ 'Pf',
                  grepl("Moderna",.) ~ 'Md',
                  TRUE ~ .)) %>% 
#                  TRUE ~ as.character(NA))) %>% 
    mutate(course = case_when(
      !is.na(d4_datetime) & Sampledate_iso > d4_datetime ~ paste0(d1_product,"-",d2_product,"-",d3_product),
      !is.na(d3_datetime) & Sampledate_iso > d3_datetime ~ paste0(d1_product,"-",d2_product,"-",d3_product),
      !is.na(d2_datetime) & Sampledate_iso > d2_datetime ~ paste0(d1_product,"-",d2_product),
      TRUE ~ d1_product)
    ) %>% 
    mutate(timeG=cut(time_since_vaccination,
                           breaks=seq(0,160,20),
                           labels=seq(20,160,20)))%>% 
    mutate(ageG=cut(ageYear,
                    breaks=seq(20,100,20))) %>% 
    mutate(BMI=cut(Q_BMI,
                   breaks=c(5,18.5,25,30,60))) %>% 
    mutate(
      n_risks = as.factor(case_when(
        is.na(n_risk_gps) ~ '0',
        n_risk_gps == '0' ~ '0',
        n_risk_gps == '1' ~ '1',
        n_risk_gps == '2' ~ '2',
        TRUE ~ '3+'
      )),
      immuno_supp = as.factor(case_when(
        severely_immuno_supp == 1 ~ 'Severely',
        immuno_supp ==1 ~ 'Yes',
        TRUE ~ 'No'
      ))
    ) %>% 
    mutate(timeG=as.numeric(as.character(timeG))) %>%
    mutate_at(c('dose','last_product','Sex','ch_resident'), as.factor) %>%
    mutate_at(vars(contains('Q_DIAG')),~ifelse(is.na(.),0,.)) %>% 
    mutate_at(vars(contains('Q_DIAG')),as.factor) %>% 
    return 
}


df_ana_pc <- get_analysis_df(df_pc)

valid_courses <- (df_ana_pc %>% group_by(course) %>% summarise(n=n()) %>% filter(n>30))$course
valid_courses

df_ana_pc <- df_ana_pc %>% filter(course %in% valid_courses)
df_ana_pc %>% group_by(course) %>% summarise(n=n())



f_normal <- ~ scale*exp(-0.5*(time-mu)^2/sd^2)
normal <- deriv(f_normal,
                namevec=c("scale","mu","sd"),
                function.arg=c("time","scale","mu","sd"))


max_time <- max(df_ana_pc$time_since_vaccination)
df_ana_pc_scaled <- df_ana_pc %>% #select(time_since_vaccination,IgG,course,last) %>%
                    mutate(IgG = IgG/max(IgG),
                           time_since_vaccination = time_since_vaccination/max(time_since_vaccination))
df_ana_pc_scaled

#startvec = c(start=100,scale=500,mu=3,sd=1)
startvec = c(scale=0.5,mu=0.1,sd=0.3)
m0 <- nls(IgG ~ normal(time_since_vaccination,scale,mu,sd),
          df_ana_pc_scaled,
          #df_ana_pc,
          start=startvec)
startvec <- coef(m0)
startvec

fix_terms <- function(df){
  df %>% mutate(
    term = case_when(
      term == 'mu' ~ 'Mode',
      term == 'scale' ~ 'Scale',
      term == 'sd' ~ 'Skew',
      T ~ term)) %>% return
}

get_nlmer_results <- function(nm){
  fixed <- fixef(nm)
  results <- broom.mixed::tidy(nm, effects = "ran_vals", conf.int = TRUE)
  results <- results %>% rowwise %>% mutate_at(c('estimate','conf.low','conf.high'), 
                                   ~ . + fixed[term][[1]]) %>% ungroup
  
  #print (results %>%
  #  mutate_at(c('estimate','conf.low','conf.high'),
  #  ~ case_when(term=='mu' ~ exp(.)*max_time,
  #              TRUE ~ .
  #              )), n=100)
  
  intercepts <- data.frame(estimate=fixed)
  intercepts$term <- rownames(intercepts)
  intercepts <- intercepts %>% fix_terms
  #print (intercepts)
  
  results %>% mutate(group = case_when(
                                        group == 'course' ~ paste0(str_count(level, "-")+1,' Dose'),
                                        group == 'ageG'~ 'Age', 
                                        group == 'n_risks'~ 'N Risks', 
                                        group == 'last_product' ~ 'Last Vaccine',
                                        TRUE ~ group
                                      ),
                     term = case_when(
                       term == 'mu' ~ 'Mode',
                       term == 'scale' ~ 'Scale',
                       term == 'sd' ~ 'Skew',
                       T ~ term
                     )
    )%>%
    ggplot() +
    geom_pointrange(aes(y=level,x=estimate,xmin=conf.low,xmax=conf.high,color=factor(term))) +
    facet_grid(group ~ term,scales='free',switch='both') +
    geom_vline(aes(xintercept=estimate),data=intercepts,linetype="dashed") +
    labs(x='',y='',color='') + 
    guides(fill="none") + 
    theme(
      legend.position="none",
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.placement = 'outside',
      axis.text.y = element_markdown(angle=0),
      strip.text.y.left = element_text(angle = 0)) %>% return
  
}


qnames <- df_ana_pc_scaled %>% select(contains('DIAG')) %>% summarise(across(where(is.factor), ~ sum(.x == 1, na.rm = TRUE))) %>%
          gather() %>% arrange(desc(value)) %>% head(5)
qnames <- qnames$key

qnames

f <- paste0(
  "IgG ~ normal(log(time_since_vaccination), scale, mu, sd) ~ (sd|last_product) + (scale|course) + (mu|last_product) + ",
  paste0("(mu|",paste0(qnames,sep = ")",collapse=' + (mu|')),
  paste0("+ (scale|",paste0(qnames,sep = ")",collapse=' + (scale|'))
)

nml1 <- nlmer(as.formula(f),
             df_ana_pc_scaled %>% filter(time_since_vaccination<200),
             verbose=1,
             start = startvec)

get_nlmer_results(nml1)






#fixef(nml1)
#time <- seq(0,1,0.01)
#y <- normal(log(time),fixef(nml1)[[1]],fixef(nml1)[[2]],fixef(nml1)[[3]])

#normal(fixef(nml1)[[2]],fixef(nml1)[[1]],fixef(nml1)[[2]],fixef(nml1)[[3]])

#*max(df_ana_pc$time_since_vaccination)
#data.frame(time=time,y=y) %>% ggplot(aes(x=time,y=y)) + geom_line() +
#      geom_vline(xintercept=exp(fixef(nml1)[[2]]),color='red',linetype='dashed') +
#      geom_hline(yintercept=fixef(nml1)[[1]],color='red',linetype='dashed')
      #geom_vline(xintercept=24.9,color='red',linetype='dashed')

#ranef(nml1)

fixed <- broom.mixed::tidy(nml1, effects = "fixed", conf.int = TRUE)
fixed






confint.merMod(nm1,
               method='boot',
               nsim = 500,
               boot.type = "basic")

fixef(nm1)




strat_var <- 'last_product'
strat_var <- 'ageG'
strat_var <- 'dose'
#strat_var <- 'n_risks'
#strat_var <- 'BMI'
#strat_var <- 'Q_DIAG_DIABETES_2'
#strat_var <- 'Sex'


df_2doses_pc <- df_ana_pc %>% filter(as.integer(as.character(dose))>1) %>%
  # & !is.na(Q_BMI) & !is.na(ageG)) %>%
  filter(!is.na(!!sym(strat_var)))

tab1 <- df_2doses_pc %>% group_by(!!sym(strat_var)) %>% summarise(n=n())
fit_2doses_pc <- fit_time(df_2doses_pc,strat_var)
fit_2doses_pc$newdata

f <- fit_2doses_pc$fits[["4"]]
f
summary(f)$coefficients[,1]
confint2(f)

results_2doses_pc <- get_results(fit_2doses_pc$fits) %>% calc_parameters
results_2doses_pc

plot_lognormal(df_2doses_pc,fit_2doses_pc$newdata,strat_var)





var <- 'dose'

fit <- gam(IgG ~ dose + s(time,by=dose,k=5), data=df_ana_pc, family=gaussian())
summary(fit)

term_list = list(dose=c(1,2,3,4),time=seq(0,200,0.1))
new_data <- expand.grid(term_list)
pred <- predict.gam(fit,new_data,se.fit = TRUE)
pred <- cbind(new_data, pred) %>% as_tibble

pred %>% ggplot(aes(x=time,y=fit)) +
         geom_ribbon(aes(ymax=fit+1.96*se.fit,
                         ymin=fit-1.96*se.fit,
                         fill=factor(dose),
                         color=factor(dose)),alpha=0.2) +
         geom_pointrange(
           aes(x=time,y=IgG,ymin=IgG-err,ymax=IgG+err,color=factor(dose)),
           data = df_ana_pc %>% group_by(dose,timeG) %>% 
             summarise(err=sd(IgG)/sqrt(n()),IgG=mean(IgG)) %>%
             mutate(time=as.integer(as.character(timeG)))
         )



df_2doses_pc <- df_ana_pc %>% filter(dose!='1') %>%
                mutate(dp = factor(paste0(as.character(dose)," - ",as.character(last_product)))) %>%
                droplevels()


fit <- gam(IgG ~ dp + ageG + Sex + s(time,by=dp,k=5), data=df_2doses_pc, family=gaussian())
summary(fit)

#term_list = list(
#  dose=levels(df_2doses_pc[['dose']]),
#  last_product=levels(df_2doses_pc[['last_product']]),
#  #dp=levels(df_2doses_pc[['dp']]),
#  Sex=levels(df_2doses_pc[['Sex']])[[1]],
#  ageG=levels(df_2doses_pc[['ageG']])[[1]],
#  time=seq(0,200,0.1)
#  )

term_list <- list()
for (term in labels(fit$terms)){
  new_term <- fit[["var.summary"]][[term]][[1]]
  term_list <- append(term_list, list(new_term))
}
names(term_list) <-  labels(fit$terms)
term_list[["time"]] <-  seq(0,200,0.1)
term_list[["dp"]] <-  levels(fit[["var.summary"]][["dp"]])

new_data <- expand.grid(term_list) %>% as_tibble
#%>% as_tibble %>%
#            mutate(dp = factor(paste0(as.character(dose)," - ",as.character(last_product))))
pred <- predict.gam(fit,new_data,se.fit = TRUE)
pred <- cbind(new_data, pred) %>% as_tibble

pred %>% group_by(dp) %>% summarise(n=n())

pred %>% ggplot(aes(x=time,y=fit)) +
  geom_ribbon(aes(ymax=fit+1.96*se.fit,
                  ymin=fit-1.96*se.fit,
                  fill=factor(dp),
                  color=factor(dp)),alpha=0.2) +
  #geom_point(aes(x=time_since_vaccination,y=IgG,color=factor(dp)),alpha=0.1,data=df_2doses_pc %>% 
  #             filter(time_since_vaccination<200)) +
  #geom_pointrange(
  #  aes(x=time,y=IgG,ymin=IgG-err,ymax=IgG+err,color=factor(dose)),
  #  data = df_2doses_pc %>% group_by(dose,timeG) %>% 
  #    summarise(err=sd(IgG)/sqrt(n()),IgG=mean(IgG)) %>%
  #    mutate(time=as.integer(as.character(timeG)))
  #) +
  #facet_grid(dose ~ last_product ,scales='free') +
  labs(y='IgG [BAU/ml]',x='Days since vaccination') +
  ylim(c(0,2100))




vars <- c(
  'dose',
  'last_product',
  'ageG',##,
  'Sex',
  'BMI',
  'n_risks'#,
  #df_ana_pc %>% select(contains('Q_DIAG')) %>% colnames
)
vars

vars <- c(
  df_ana_pc %>% select(contains('Q_DIAG')) %>% colnames
)
vars



results <- list()
plots <- list()

for (strat_var in vars){
  
  df_2doses_pc <- df_ana_pc %>% filter(as.integer(as.character(dose))>1) %>%
    # & !is.na(Q_BMI) & !is.na(ageG)) %>%
    filter(!is.na(!!sym(strat_var)))
  
  tab1 <- df_2doses_pc %>% group_by(!!sym(strat_var)) %>% summarise(n=n())
  print (tab1)
  print (strat_var)


  fit_2doses_pc <- tryCatch({
    fit_time(df_2doses_pc,strat_var)
  },
  error=function(e){
    message(e)
    return (NULL)
  }
  )
  if(is.null(fit_2doses_pc)){
    next
  }
  
  results_2doses_pc <- get_results(fit_2doses_pc$fits) %>% calc_parameters
  
  results[[strat_var]] <- results_2doses_pc
  plots[[strat_var]] <- plot_lognormal(df_2doses_pc,fit_2doses_pc$newdata,strat_var)
   
}



pdata <- NULL
for (variable in labels(results)){
  
  res <- results[[variable]]
  
  temp <- res %>% rename(y=fit,x=estimate_a,xmax=`97.5 %_a`,xmin=`2.5 %_a`) %>%
        mutate(var='scale') %>%
        select(var,y,x,xmin,xmax) %>%
        rbind(
         res %>% 
         rename(y=fit,x=mode,xmax=mode_ucl,xmin=mode_lcl) %>%
         mutate(var='mode',x=x*100,xmin=100*xmin,xmax=100*xmax) %>%
         select(var,y,x,xmin,xmax)
        ) %>%
        rbind(
          res %>% 
            rename(y=fit,x=skew,xmax=skew_ucl,xmin=skew_lcl) %>%
            mutate(var='skew',x=x,xmin=xmin,xmax=xmax) %>%
            select(var,y,x,xmin,xmax)
        )
  
  print (variable)
  temp <- temp %>% mutate(variable=variable)
  #      mutate(variable=variable)
  pdata <- rbind(pdata,temp)

}
  

df_2doses_pc <- df_ana_pc %>% filter(as.integer(as.character(dose))>1)
fit <-nls(IgG ~ a*(dlnorm(time/100, mean,sd)), 
          data = df_2doses_pc, 
          start = list(a=700,mean=1,sd=1.1),
          lower = list(a=0,mean=0,sd=0),
          upper = list(a=100000,mean=5,sd=5),
          algorithm = 'port')
nominal <- summary(fit)$coefficients[,1]


nominal_scale <- nominal[[1]]
mean <- nominal[[2]]
sd <- nominal[[3]]

nominal_mode = exp(mean - sd^2)
nominal_skew = (exp(sd^2)+2)*sqrt(exp(sd^2)-1)
nominal_skew


#plots[[strat_var]] <- plot_lognormal(df_2doses_pc,fit_2doses_pc$newdata,strat_var)






pdata <- pdata  %>% 
  mutate(
    variable = case_when(
      variable == 'dose' ~ 'Dose',
      variable == 'ageG' ~ 'Age',
      variable == 'last_product' ~ 'Last Vaccine',
      variable == 'n_risks' ~ 'Number of Risk Groups',
      TRUE ~ variable)
    )

pdataQ <- pdata %>% 
  mutate(
    y = ifelse(y==1,eavehelpers::get_label(variable),NA),
    variable = 'QCOVID'
  ) %>% filter(!is.na(y)) %>% group_by(y) %>% 
  filter(!any(xmax/x > 100))
pdataQ


scales_x <- list(
  mode = scale_x_continuous(lim=c(10,70)),
  scale = scale_x_continuous(lim=c(0,10000)),
  #scale = scale_x_log10(lim=c(nominal_scale*0.1,nominal_scale*10)),
  skew = scale_x_log10(lim=c(1,500))
)

p <- pdataQ %>%
  ggplot() +
  geom_pointrange(aes(y=y,x=x,xmin=xmin,xmax=xmax,
                      color=factor(var))) +
  geom_vline(data=filter(pdataQ, var=="mode"), aes(xintercept=nominal_mode*100),linetype="dotted") + 
  geom_vline(data=filter(pdataQ, var=="scale"), aes(xintercept=nominal_scale),linetype="dotted") + 
  geom_vline(data=filter(pdataQ, var=="skew"), aes(xintercept=nominal_skew),linetype="dotted") + 
  guides(fill="none",color="none") +
  labs(x='',y='',color='') +
  facet_grid_sc(variable ~ var,
                #scale='free',
                switch='both',
                scales=list(x=scales_x,y='free')) + #space='free'
  theme(
    #panel.background = element_rect(fill = NA, color = "black"),
    strip.background = element_blank(),
    strip.placement = 'outside',
    axis.text.y = element_markdown(angle=0),
    strip.text.y.left = element_text(angle = 0)) +
  scale_colour_discrete_phs(palette='all') 

p

plots[['Q_DIAG_CKD3']] + labs(y='IgG [BAU/ml]',x='Days since vaccination', color='CKD Level 3') +
  theme_classic() + xlim(5,160)

#p %>% 
#  ggsave(filename = "test.png",width = 10., height = 15., dpi = 300)

p %>% 
  ggsave(filename = "test.png",width = 10., height = 7., dpi = 300)


var <- 'dose'

plots[[var]] +
  labs(x='Days since vaccination',
       y='Average IgG [BAU/ml]',
       color='Last Vaccine Dose'
       ) +
  scale_x_log10()

  ggsave(filename = paste0(var,".png"),width = 10., height = 3., dpi = 300)


  
plots[['last_product']] +
    labs(x='Days since vaccination',
         y='Average IgG [BAU/ml]',
         color='Last Vaccine Product'
    ) 

plots[['n_risks']] +
  labs(x='Days since vaccination',
       y='Average IgG [BAU/ml]',
       color='Number of Risks'
  ) 

plots[['BMI']] +
  labs(x='Days since vaccination',
       y='Average IgG [BAU/ml]',
       color='Body Mass Index'
  ) 


plots[['ageG']] +
  labs(x='Days since vaccination',
       y='Average IgG [BAU/ml]',
       color='Age [Years]'
  ) 



plots[['Sex']] +
  labs(x='Days since vaccination',
       y='Average IgG [BAU/ml]',
       color='Sex'
  ) 



plot_results(results_2doses_pc,'scale') +
  scale_y_log10() +
  scale_x_log10()



df_2doses_bd <- df_ana_bd %>% filter(as.integer(as.character(dose))>1)  # & !is.na(Q_BMI) & !is.na(ageG)) %>%
fit_2doses_bd <- fit_time(df_2doses_bd,'last_product')
results_2doses_bd <- get_results(fit_2doses_bd$fits) %>% calc_parameters

results <- results_2doses_pc %>% mutate(cohort='pc') %>% 
  rbind(results_2doses_bd %>% mutate(cohort='bd'))

p <- plot_results(results) +
  scale_y_log10() +
  scale_x_log10()
p



df_2doses_bd %>% 
  group_by(time,!!sym(strat_var)) %>% 
  summarise(igg=mean(IgG),err=sd(IgG)/sqrt(n())) %>% ungroup %>%
  ggplot(aes(x=time,y=igg)) +
  geom_pointrange(aes(ymin=igg-err,ymax=igg+err,color=!!sym(strat_var))) +
  geom_line(aes(x=time,y=IgG,color=strat_var),
            data=fit_2doses_bd$newdata) +
  #ylim(c(0,2100)) +
  xlim(c(0,150))



df_temp <- df_2doses_bd %>% 
  group_by(time,!!sym(strat_var)) %>% 
  summarise(igg=mean(IgG),err=sd(IgG)/sqrt(n())) %>% ungroup %>%
  mutate(cohort='bd') %>% rbind(
    df_2doses_pc %>% 
      group_by(time,!!sym(strat_var)) %>% 
      summarise(igg=mean(IgG),err=sd(IgG)/sqrt(n())) %>% ungroup %>%
      mutate(cohort='pc')
  )

df_temp %>% 
  ggplot(aes(x=time,y=igg)) +
  geom_pointrange(aes(ymin=igg-err,ymax=igg+err,color=!!sym(strat_var))) +
  #geom_line(aes(x=time,y=IgG,color=strat_var),
  #          data=fit_2doses_bd$newdata) +
  #ylim(c(0,2100)) +
  xlim(c(0,150)) +
  facet_grid(cohort ~ .,scales = "free")#,space = 'free')





