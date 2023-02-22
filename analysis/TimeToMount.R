
library(ggplot2)
library(mgcv)
library(phsstyles)
library(facetscales)
library(eavehelpers)
library(lme4)
library(broom.mixed)
library("stringr")
library(ggtext)
devtools::document("SerologyAnalysis/eavehelpers")
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
                    breaks=seq(10,100,20))) %>% 
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

fit_time <- function(df,strat_var) {
  time_max <- 160
  df <- df %>% filter(time<time_max) %>% droplevels()
  newdata <- NULL 
  fits <- list()
  for(i in 1:length(levels(df[[strat_var]]))){
    
    temp <- df %>% filter(as.numeric(!!sym(strat_var))==i)
    
    fit <-nls(IgG ~ a*(dlnorm(time/100, mean,sd)), 
              data = temp %>% as.data.frame, 
              start = list(a=700,mean=1,sd=1.1),
              lower = list(a=0,mean=0,sd=0),
              upper = list(a=100000,mean=5,sd=5),
              algorithm = 'port')
    
    #fit <-nls(IgG ~ init + a*(dlnorm(time/100, mean,sd)), 
    #          data = temp, 
    #          start = list(init=100,a=700,mean=1,sd=1.1),
    #          lower = list(init=0,a=0,mean=0,sd=0),
    #          upper = list(init=100000, a=100000,mean=5,sd=5),
    #          algorithm = 'port')
    
    
    
    fits[[levels(df[[strat_var]])[[i]]]] = fit
    
    temp <- data.frame(time=seq(0,time_max,0.5),
                       strat_var=levels(df[[strat_var]])[[i]])
    temp$IgG <- predict(fit,newdata=temp)
    boot <- nlsBoot(fit, niter = 10000)
    #nlstools::nlsBootPredict(boot,newdata = temp, interval = "confidence") %>% print
    newdata <- rbind(newdata,temp)
  }
  return (list(newdata=newdata,fits=fits));
}

get_results <- function(fits){
  results <- NULL
  for(name in names(fits)){
    f <- fits[[name]]
    err <- nlstools::confint2(f)
    pred <- data.frame(estimate=summary(f)$coefficients[,1])
    res <- pred %>% cbind(err) 
    res <- res %>% mutate(var=rownames(res),fit=name)
    results <- rbind(results,res)
    
  }
  results <- results %>% 
    pivot_wider(names_from=var,
                values_from=c('estimate','2.5 %','97.5 %')
    ) %>% return 
}

calc_parameters <- function(results){
  results <- results %>% mutate(
    mode = exp(estimate_mean - estimate_sd^2),
    #mode_cl1 = exp(`97.5 %_mean` - `2.5 %_sd`^2),
    #mode_cl2 = exp(`2.5 %_mean` - `97.5 %_sd`^2),
    mode_cl3 = exp(`97.5 %_mean` - `97.5 %_sd`^2),
    mode_cl4 = exp(`2.5 %_mean` - `2.5 %_sd`^2),
    skew = (exp(estimate_sd^2)+2)*sqrt(exp(estimate_sd^2)-1),
    skew_ucl = (exp(`97.5 %_sd`^2)+2)*sqrt(exp(`97.5 %_sd`^2)-1),
    skew_lcl = (exp(`2.5 %_sd`^2)+2)*sqrt(exp(`2.5 %_sd`^2)-1),
    #mode_lcl = exp(`97.5 %_mean` - `97.5 %_sd`^2),
    #mode_ucl = exp(`2.5 %_mean` - `2.5 %_sd`^2),
    scale = estimate_a,
    scale_lcl =`2.5 %_a`,
    scale_ucl=`97.5 %_a`
  ) %>%
    mutate(mode_ucl = pmax(mode_cl3,mode_cl4)) %>%
    mutate(mode_lcl = mode - (mode_ucl - mode)) %>%
    return 
  
}

plot_results <- function(results,y='skew'){
  
  results %>% ggplot(aes(
    x=100*mode,
    y=!!sym(y),
    color=as.factor(fit))) +
    geom_pointrange(aes(
      ymin=!!sym(paste0(y,"_lcl")),
      ymax=!!sym(paste0(y,"_ucl"))
    )) +
    geom_errorbarh(
      aes(
        xmin=100*mode_ucl,
        xmax=100*mode_lcl,
        height=0
      )) +
    #scale_x_log10(limits=c(70,100000)) +
    #scale_y_log10() +
    labs(
      x='Mode',
      y='Skewness',
      color=''
    ) +
    scale_color_brewer(palette="Set1") +
    theme_classic() %>% return
}

plot_lognormal <- function(df,newdata,strat_var){
  df %>% 
    group_by(time,!!sym(strat_var)) %>% 
    summarise(igg=mean(IgG),err=sd(IgG)/sqrt(n())) %>% ungroup %>%
    ggplot(aes(x=time,y=igg)) +
    geom_pointrange(aes(ymin=igg-err,ymax=igg+err,color=!!sym(strat_var))) +
    geom_line(aes(x=time,y=IgG,color=strat_var),
              data=newdata) +
    #ylim(c(0,2100)) +
    xlim(c(0,160)) +
    scale_colour_discrete_phs(palette = 'all') %>% 
    return
}

fix_terms <- function(df){
  df %>% mutate(
    term = case_when(
      term == 'mu' ~ 'Mode [days]',
      term == 'scale' ~ 'Scale [BAU/ml]',
      term == 'sd' ~ 'Skew [days]',
      T ~ term)) %>% return
}

get_nlmer_results <- function(nm,modify=F){
  fixed <- fixef(nm)
  results <- broom.mixed::tidy(nm, effects = "ran_vals", conf.int = TRUE)
  
  if (modify){
    results <- results %>% rowwise %>% mutate_at(c('estimate','conf.low','conf.high'), 
                                                 ~ . + fixed[term][[1]]) %>% ungroup
  }
  
  #results <- results %>%
  #  mutate_at(c('estimate','conf.low','conf.high'),
  #      ~ case_when(term=='mu' ~ exp(.)*max_time,
  #                  TRUE ~ .
  #      )
  #)
  
  intercepts <- data.frame(estimate=fixed)
  intercepts$term <- rownames(intercepts)
  intercepts <- intercepts #%>% fix_terms
  
  return (list(results=results,intercepts=intercepts))
}

fix_group_names <- function(df){
  df %>% mutate(group = case_when(
    group == 'course' ~ paste0(str_count(level, "-")+1,' Dose'),
    group == 'ageG'~ 'Age', 
    group == 'n_risks'~ 'N Risks', 
    group == 'last_product' ~ 'Last Vaccine',
    TRUE ~ group
  )) %>% return
}

plot_nlmer_results <- function(results,intercepts=NULL){
  df <- results %>% fix_terms
  
  p <- df %>% ggplot() +
    geom_pointrange(aes(y=level,x=estimate,xmin=conf.low,xmax=conf.high,color=factor(term))) +
    facet_grid( group ~ term ,scales='free',switch='both')
  
  if(!is.null(intercepts)){
    p <- p + geom_vline(aes(xintercept=estimate),
                        data=intercepts %>% fix_terms,linetype="dashed")
  }
  
  p <- p + labs(x='',y='',color='') + 
    guides(fill="none") + 
    scale_colour_discrete_phs(palette = 'all') +
    theme(
      legend.position="none",
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.placement = 'outside',
      axis.text.y = element_markdown(angle=0),
      strip.text.y.left = element_text(angle = 0))
  
  return (p);
}



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


df_pc <- readRDS("/conf/EAVE/GPanalysis/data/temp/serology_Feb2023.rds")



df_ana_pc <- get_analysis_df(df_pc)

valid_courses <- (df_ana_pc %>% group_by(course) %>% summarise(n=n()) %>% filter(n>30))$course
valid_courses

df_ana_pc <- df_ana_pc %>% filter(course %in% valid_courses)
df_ana_pc %>% group_by(course) %>% summarise(n=n())
df_ana_pc %>% group_by(ageG) %>% summarise(n=n())


df_ana_bd <- get_analysis_df(df_bd) %>% filter(!is.na(IgG))
df_ana_bd <- df_ana_bd %>% filter(course %in% valid_courses)
df_ana_bd %>% group_by(course) %>% summarise(n=n())


#startvec <- c(scale = 1000, mu = log(50), sd=log(3))
#myfunc <- ~ scale * (0.5*(1+tanh(1e6*(time-mu)))*exp(lambda*(time-mu))*(m*mu+c) + 0.5*(1-tanh(1e6*(time-mu)))*(m*time + c))
#nfun <- deriv(myfunc,
#              namevec=c("scale","m","c","mu","lambda"),
#              function.arg=c("time","scale","m","c","mu","lambda"))
#startvec <- c(scale=700, m = 2, c=10, mu=40, lambda=-0.01)

f_normal <- ~  scale*(exp(-0.5*(time-mu)^2/sd^2))
normal <- deriv(f_normal,
                namevec=c("scale","mu","sd"),
                function.arg=c("time","scale","mu","sd"))


max_time_pc <- max(df_ana_pc$time_since_vaccination)
max_igg_pc <- max(df_ana_pc$IgG)

df_ana_pc_scaled <- df_ana_pc %>% ungroup %>% filter(time_since_vaccination<200 & dose!=1) %>% 
  mutate(IgG = IgG/max_igg_pc,
         time_since_vaccination = time_since_vaccination/max_time_pc)


#max_time_bd <- max(df_ana_bd$time_since_vaccination)
#max_time_bd
#max_igg_bd <- max(df_ana_bd$IgG)
#max_igg_bd

#df_ana_bd_scaled <- df_ana_bd %>% ungroup %>% filter(time_since_vaccination<200) %>% 
#  mutate(IgG = IgG/max_igg_bd,
#         time_since_vaccination = time_since_vaccination/max_time_bd)


startvec = c(scale=0.25,mu=-3,sd=1)
m0 <- nls(IgG ~ normal(log(time_since_vaccination),scale,mu,sd),
          df_ana_pc_scaled,
          start=startvec
          #lower=c(start=0,scale=0,mu=-10,sd=-10),
          )
startvec <- coef(m0)
startvec


nlme_pc_1 <- nlmer(IgG ~ normal(log(time_since_vaccination), scale, mu, sd) ~ 
                (sd|last_product) + (scale|course) + (mu|last_product)
              + (scale|ageG) + (sd|ageG) + (mu|ageG) 
              #+ (scale|Sex) 
              + (sd|Sex) + (mu|Sex)
              + (scale|n_risks) + (sd|n_risks) + (mu|n_risks),
              #+ (scale|BMI) + (sd|BMI) + (mu|BMI),#  + (sd|course), #(mu|course)
             df_ana_pc_scaled,
             #df_ana_pc %>% filter(time_since_vaccination<200),# %>% filter(dose!=1) %>% droplevels,
             #control = nlmerControl(optCtrl=list(maxfun=100000)),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             verbose=1,
             start = startvec)

summary(nlme_pc_1)



nlme_pc_2 <- nlmer(IgG ~ normal(log(time_since_vaccination), scale, mu, sd) ~ 
                    # (sd|last_product) + (scale|course) + (mu|last_product)
                   + (ageG + Sex | scale) + (ageG + Sex + n_risks + last_product | mu),
                   #+ (scale|Sex) 
                   #+ (sd|Sex) + (mu|Sex)
                   #+ (scale|n_risks) + (sd|n_risks) + (mu|n_risks),
                   #+ (scale|BMI) + (sd|BMI) + (mu|BMI),#  + (sd|course), #(mu|course)
                   df_ana_pc_scaled,
                   #df_ana_pc %>% filter(time_since_vaccination<200),# %>% filter(dose!=1) %>% droplevels,
                   #control = nlmerControl(optCtrl=list(maxfun=100000)),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                   verbose=1,
                   start = startvec)

summary(nlme_pc_2)
nm <- nlme_pc_2

fixed <- fixef(nm)
results <- broom.mixed::tidy(nm, effects = "ran_vals", conf.int = TRUE)  %>% 
           rowwise %>%
           mutate_at(c('estimate','conf.low','conf.high'), 
                                               ~ . + fixed[group][[1]]) %>% 
           ungroup %>% filter(group=='mu') %>%
           mutate_at(c('estimate','conf.low','conf.high'), ~ max_time_pc*exp(.))
           
results %>% ggplot(aes(y=term,x=estimate)) +
            geom_pointrange(aes(xmin=conf.low,xmax=conf.high,color=level))




results <- get_nlmer_results(nlme_pc_1,modify=T)
results

res <- results$results


levels <- unique(res$level)
levels

df <- NULL
for (lev in levels){
  
  mu <- res %>% filter(level==lev & term=='mu') #%>% as.list()
  if(nrow(mu)==0){
    temp <- fixef(nlme_pc_1)['mu'][1]
    mu <- list(estimate=temp,conf.low=temp,conf.high=temp)
  }
  else{
    mu <- mu %>% as.list()
  }
  
  sd <- res %>% filter(level==lev & term=='sd') #%>% as.list()
  if(nrow(sd)==0){
    temp <- fixef(nlme_pc_1)['sd'][1]
    sd <- list(estimate=temp,conf.low=temp,conf.high=temp)
  }
  else{
    sd <- sd %>% as.list()
  }
  
  scale <- res %>% filter(level==lev & term=='scale') #%>% as.list()
  if(nrow(scale)==0){
    temp <- fixef(nlme_pc_1)['scale'][1]
    scale <- list(estimate=temp,conf.low=temp,conf.high=temp)
  }
  else{
    scale <- scale %>% as.list()
  }
  

  print (lev)
  group <- scale$group
  if(is.null(group)){
    group <- mu$group
  }
  
  
  t <- seq(-10,10,0.01)
  y <- normal(t,scale$estimate,mu$estimate,sd$estimate)
  ymin <- normal(t,scale$conf.low,mu$conf.low,sd$conf.low)
  ymax <- normal(t,scale$conf.high,mu$conf.high,sd$conf.high)
  
  temp <- data.frame(t=max_time_pc*exp(t),
                     y=max_igg_pc*y,
                     ymin=max_igg_pc*ymin,
                     ymax=max_igg_pc*ymax,
                     level=lev,
                     group=group) %>% as_tibble
  #temp <- data.frame(t=t,y=y,ymin=ymin,ymax=ymax,color=lev) %>% as_tibble
  
  df <- df %>% rbind(temp)
}

df %>%
  ggplot(aes(y=y,x=t,color=level,fill=level)) + 
  geom_line(linetype='dashed') +
  geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.5) + 
  labs(x='Days Since Vaccine',y='IgG [BAU/ml]') +
  xlim(5,200) +
  facet_wrap(group ~ level) +
  guides(fill="none",color="none") 



scaled_to_normal <- function(term,.){
  case_when(
    term=='mu' ~ max_time_pc*exp(.),
    term=='sd' ~ exp(.),
    term=='scale' ~ max_igg_pc*.,
    TRUE ~ .
  )
}
results$results <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),
                                                 ~ scaled_to_normal(term,.)
)
results$intercepts <- results$intercepts %>% mutate_at(c('estimate'),~scaled_to_normal(term,.))

res <- results$results %>% fix_group_names
res$group <- factor(res$group,levels=c("2 Dose","3 Dose","Last Vaccine","N Risks","Age","Sex"))

plot_nlmer_results(res,results$intercepts) 









results <- get_nlmer_results(nlme_pc_1,modify=F)
results

t0 <- exp(fixef(nlme_pc_1)[['mu']])*max_time_pc
t0

scaled_to_normal <- function(term,.){
  case_when(
    term=='mu' ~ t0*(exp(.) - 1) ,
    term=='sd' ~ t0*(exp(.) - 1),
    term=='scale' ~ max_igg_pc*(.),
    TRUE ~ .
  )
}

results$results <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),
                                                 ~ scaled_to_normal(term,.)
)

plot_nlmer_results(results$results) + 
  geom_vline(aes(xintercept=0),linetype="dashed")














scaled_to_normal <- function(term,.){
  case_when(
    term=='mu' ~ max_time_pc*exp(.),
    term=='sd' ~ exp(.),
    #term=='scale' ~ max_igg_pc*.,
    TRUE ~ .
  )
}

results$results <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),
                                                 ~ scaled_to_normal(term,.)
)
results$intercepts <- results$intercepts %>% mutate_at(c('estimate'),~scaled_to_normal(term,.))

results$intercept


plot_nlmer_results(results)




f_lognormal <- ~ scale*exp(-0.5*(log(time)-mu)^2/sd^2)
lognormal <- deriv(f_lognormal,
                   namevec=c("scale","mu","sd"),
                   function.arg=c("time","scale","mu","sd"))


startvec = c(scale=1000,mu=3,sd=-2)
m0 <- nls(IgG ~ lognormal(time_since_vaccination,scale,mu,sd),
          df_ana_pc %>% ungroup %>% filter(time_since_vaccination<200),
          start=startvec
)
startvec <- coef(m0)
startvec


nlme_pc_1 <- nlmer(IgG ~ lognormal(time_since_vaccination, scale, mu, sd) ~ 
                     (sd|last_product) + (scale|course) + (mu|last_product)
                   + (scale|ageG) + (sd|ageG) + (mu|ageG) 
                   + (scale|Sex) + (sd|Sex) + (mu|Sex)
                   + (scale|n_risks) + (sd|n_risks) + (mu|n_risks),
                   #+ (scale|BMI) + (sd|BMI) + (mu|BMI),#  + (sd|course), #(mu|course)
                   df_ana_pc %>% ungroup %>% filter(time_since_vaccination<200),
                   control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
                   verbose=1,
                   start = startvec)



results <- get_nlmer_results(nlme_pc_1,modify=F)

plot_nlmer_results(results$results)





startvec = c(scale=0.25,mu=-2,sd=1)
m0 <- nls(IgG ~ normal(log(time_since_vaccination),scale,mu,sd),
          df_ana_bd_scaled,
          start=startvec
)
startvec <- coef(m0)
startvec


nlme_bd_1 <- nlmer(IgG ~ normal(log(time_since_vaccination), scale, mu, sd) ~ 
                (sd|last_product) + (scale|course) + (mu|last_product)
              + (scale|ageG) + (sd|ageG) + (mu|ageG)
              + (scale|Sex) + (sd|Sex) + (mu|Sex)
              + (scale|n_risks) + (sd|n_risks) + (mu|n_risks),
              #+ (scale|BMI) + (sd|BMI) + (mu|BMI),#  + (sd|course), #(mu|course)
              df_ana_bd_scaled,
              #df_ana_pc %>% filter(time_since_vaccination<200),# %>% filter(dose!=1) %>% droplevels,
              #control = nlmerControl(optCtrl=list(maxfun=100000)),
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              verbose=1,
              start = startvec)

results <- get_nlmer_results(nlme_bd_1)

results$results <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),
                                                 ~ scaled_to_normal(term,.)
)
results$intercepts <- results$intercepts %>% mutate_at(c('estimate'),~scaled_to_normal(term,.))

results$intercept

plot_nlmer_results(results)

results_bd <- results


results_full <- results_pc$results %>% mutate(cohort='Primary Care') %>% rbind(
  results_bd$results %>% mutate(cohort='Blood Donors')
)
results_full


plot_nlmer_results(list(results=results_full,intercepts=results_pc$intercepts))


startvec <- fixef(nml1)
startvec

df_ana_pc_scaled %>% colnames

nml2 <- nlmer(IgG ~ normal(log(time_since_vaccination), scale, mu, sd) ~ 
                (sd|last_product) + (scale|course) + (mu|last_product)
              + (scale|ageG) + (sd|ageG) + (mu|ageG)
              + (scale|Sex) + (sd|Sex) + (mu|Sex)
              + (scale|immuno_supp) + (sd|immuno_supp) + (mu|immuno_supp)
              + (scale|n_risks),# + (sd|n_risks) + (mu|n_risks)
              #+ (scale|BMI) + (sd|BMI) + (mu|BMI),
              df_ana_pc_scaled,
              #control = nlmerControl(optCtrl=list(maxfun=100000)),
              #control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              verbose=1,
              start = startvec)

get_nlmer_results(nml2)


startvec <- fixef(nml2)
startvec['start'] = 100

nlme_1dose <- nlmer(IgG ~ start + normal(log(time_since_vaccination), scale, mu, sd) ~ 
                (start|course) + 
                (sd|last_product) + (scale|course) + (mu|last_product)
              + (scale|ageG) + (sd|ageG) + (mu|ageG)
              + (scale|Sex) + (sd|Sex) + (mu|Sex)
              + (scale|immuno_supp) + (sd|immuno_supp) + (mu|immuno_supp)
              + (scale|n_risks),# + (sd|n_risks) + (mu|n_risks)
              #+ (scale|BMI) + (sd|BMI) + (mu|BMI),
              df_ana_pc_scaled,# %>% filter(dose!=1) ,
              control = nlmerControl(optimizer='bobyqa'),#optCtrl=list(maxfun=100000)),
              #control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
              verbose=1,
              start = startvec)


results <- get_nlmer_results(nlme_1dose)

levels <- unique((results$results %>% filter(group=='ageG'))$level)

df <- NULL
for (lev in levels){
  mu <- results$results %>% filter(level==lev & term=='mu') %>% as.list()
  sd <- results$results %>% filter(level==lev & term=='sd') %>% as.list()
  scale <- results$results %>% filter(level==lev & term=='scale') %>% as.list()
  
  t <- seq(-10,10,0.01)
  y <- normal(t,scale$estimate,mu$estimate,sd$estimate)
  ymin <- normal(t,scale$conf.low,mu$conf.low,sd$conf.low)
  ymax <- normal(t,scale$conf.high,mu$conf.high,sd$conf.high)
  

  temp <- data.frame(t=max_time*exp(t),y=max_igg*y,ymin=max_igg*ymin,ymax=max_igg*ymax,color=lev) %>% as_tibble
  
  df <- df %>% rbind(temp)
}

df %>%
  ggplot(aes(y=y,x=t,color=color,fill=color)) + 
  geom_line(linetype='dashed') +
  geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.5) + 
  xlim(0,200)




scaled_to_normal <- function(term,.){
  case_when(
      term=='mu' ~ max_time*exp(.),
      term=='sd' ~ exp(.),
      term=='scale' ~ max_igg*.,
      TRUE ~ .
  )
}

results$results <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),
  ~ scaled_to_normal(term,.)
)
results$intercepts <- results$intercepts %>% mutate_at(c('estimate'),~scaled_to_normal(term,.))

results$intercept

plot_nlmer_results(results)



fixed <- fixef(nlme_1dose)
t <- seq(-10,10,0.01)
data.frame(t=t,y=normal(t,fixed['scale'][1],fixed['mu'][1],fixed['sd'][1])) %>% as_tibble %>%
  ggplot(aes(y=y,x=max_time*exp(t))) + geom_line() + xlim(0,200)



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





