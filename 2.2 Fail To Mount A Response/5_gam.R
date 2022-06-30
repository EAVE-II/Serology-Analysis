
library(mgcv)
library("oddsratio")
library(ggplot2)
library(tidymv)
#install.packages("tidymv")

perform_gam <- function(model,formula){
  
  fit <- gam(formula,family=binomial,data=model)

  return(fit)
  tab <- exp(cbind(coef(glmFit), coef(glmFit))) %>% as.data.frame
  tab <- tab[-1,]
  tab$names <- rownames(tab)
  
  
  names(tab)[1] <- 'OR'
  
  return ()
}

get_pterm_or_from_gam <- function (fit,model) {
  pterms <- labels(fit$pterms)
  or <- NULL
  for(term in pterms){
    levs <- levels(model[[term]])
    if (is.null(levs)){
      levs <- unique(sort(model[[term]]))
    }
    else{
      levs <- levels(droplevels(model[[term]]))
    }
    ref <- levs[1]
    values <- levs[2:length(levs)]
    for (value in values){
      tor <- or_gam(data = model, model = fit,pred = term, values=c(ref,value))
      or <- or %>% rbind (tor)
    }
  }
  or <- or %>% unite('names',sep='',c(predictor,value2)) %>% select(-value1)
  names(or)[2] <- 'OR'
  names(or)[3] <- 'LCL'
  names(or)[4] <- 'UCL'
  or <- or %>% mutate( names=ifelse(is.na(as.character(lookup[names])),
                                    ifelse(is.na(as.character(lookup[gsub('.{1}$', '', names)])),names,lookup[gsub('.{1}$', '', names)]),
                                    lookup[names]))
  return (or)
}


modelA <- modelA %>% filter(stage!=1)

formula <- as.formula(paste0("outcome ~ s(ageYear) + Sex + n_risk_gps + product_binary + stage"))
formula <- as.formula(paste0("outcome ~ s(age) + s(bmi) + Sex + simd2020v2_sc_quintile + n_risk_gps  + s(days_since_first_measurement) + s(days_since_last_vac) + s(days_between_d1_d2)  + product_binary"))
formula <- as.formula(paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + product_binary"))

modelB.formula

model <- modelB %>% droplevels() #%>% filter(product_binary == 0) #%>% filter(stage!=1)



model <- model %>% mutate(vacc=as.factor(ifelse(product_binary==1,'AZ','Pfizer')))
model$cat <- factor(model$cat, levels=c("Pf-Pf","AZ","Pf","AZ-AZ","mixed 2 doses","AZ-AZ-AZ","Pf-Pf-Pf","mixed 3 doses"))
model %>% group_by(cat) %>% summarise(n=n())
#,
#                          ch_resident=as.factor(ch_resident),
#                          Q_DIAG_BLOOD_CANCER=as.factor(Q_DIAG_BLOOD_CANCER))


#vars <- model %>% select(-age,-days_since_first_measurement,-days_since_last_vac,-outcome) %>% colnames()
#model %>% group_by(Q_DIAG_BLOOD_CANCER) %>% summarise(n=n())


model$days_since_last_vac = factor(cut(model$days_since_last_vac, breaks = c(0,30,50,100,150,200,10000), right = T, 
                                               labels = c("15-30","30-50","50-100","100-150","150-200","200+")),
                                               levels=c("50-100","15-30","30-50","100-150","150-200","200+"))

model$days_since_first_measurement = factor(cut(model$days_since_first_measurement, breaks = c(0,100,200,300,1000), right = T, 
                                       labels = c("0-100","100-200","100-300","300-400")),
                                       levels=c("100-200","0-100","100-300","300-400"))



formula <- "outcome ~ ageYear + Sex  +  Q_BMI + simd2020v2_sc_quintile + n_risk_gps + 
                      Q_DIAG_CKD_LEVEL + Q_DIAG_AF + Q_DIAG_BLOOD_CANCER + Q_DIAG_CCF + Q_DIAG_CHD +
                      Q_DIAG_CIRRHOSIS + Q_DIAG_CONGEN_HD + Q_DIAG_COPD + Q_DIAG_DEMENTIA +
                      Q_DIAG_DIABETES_1 + Q_DIAG_FRACTURE + Q_DIAG_IMMU + Q_DIAG_NEURO +
                      Q_DIAG_PARKINSONS + Q_DIAG_PULM_RARE + Q_DIAG_PVD + Q_DIAG_RESP_CANCER + 
                      Q_DIAG_SEV_MENT_ILL + Q_DIAG_STROKE + Q_DIAG_VTE +
                      shielding + ch_resident + immuno_supp + 
                      prior_infection + days_since_last_vac + days_since_first_measurement +
                      cat " #+product_binary  + 

formula
#formula <- "outcome ~ ageYear + Q_DIAG_BLOOD_CANCER"
or <- perform_glm(model,formula,do_adjusted=FALSE,do_unadjusted = TRUE)
or

model %>% group_by(simd2020v2_sc_quintile,outcome) %>% summarise(n=n()) 




pglm <- plot_fit(or,"test",unadjusted=TRUE,adjusted=FALSE,xmax=150)
pglm

formula <- as.formula("outcome ~ s(age,by=cat) + Sex  +  Q_BMI + simd2020v2_sc_quintile + Q_DIAG_CKD_LEVEL + Q_DIAG_AF + Q_DIAG_BLOOD_CANCER + Q_DIAG_CCF + Q_DIAG_CHD +
                      Q_DIAG_CIRRHOSIS + Q_DIAG_CONGEN_HD + Q_DIAG_COPD + Q_DIAG_DEMENTIA +
                      Q_DIAG_DIABETES_1 + Q_DIAG_FRACTURE + Q_DIAG_IMMU + Q_DIAG_NEURO +
                      Q_DIAG_PARKINSONS + Q_DIAG_PULM_RARE + Q_DIAG_PVD + Q_DIAG_RESP_CANCER + 
                      Q_DIAG_SEV_MENT_ILL + Q_DIAG_STROKE + Q_DIAG_VTE +
                      shielding + ch_resident + immuno_supp + 
                      s(days_since_first_measurement) + prior_infection + s(days_since_last_vac) +
                      cat ") #+product_binary
#model <- model %>% mutate(vacc=as.factor(ifelse(product_binary==1,'AZ','Pfizer')))
#formula <- as.formula("outcome ~ s(age,by=cat) + s(days_since_first_measurement) + cat + Q_BMI + Sex + Q_DIAG_BLOOD_CANCER  + ch_resident ") 
gamFit <- perform_gam(model,formula)


or <- get_pterm_or_from_gam(gamFit,model)
or
pgam <- plot_simple(or,xmax=150)
pgam



term_list <- list()
for (term in labels(gamFit$terms)){
    new_term <- gamFit[["var.summary"]][[term]][[1]]
    term_list <- append(term_list, list(new_term))
}

names(term_list) <-  labels(gamFit$terms)

term_list[['age']] <- seq(0,100,1)
term_list[['cat']] <- levels(gamFit[["var.summary"]][['cat']][[1]])

new_data <- expand.grid(term_list)
new_data

pred <- predict.gam(gamFit,new_data,se.fit = TRUE)
pred <- cbind(new_data, pred)
pred <- pred %>% mutate(cat=lookup[paste0('cat',cat)])


func <- function(v) {
  #return(exp(v));
  return (v);
}

ggplot(pred, aes(x=age, y=func(fit), color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(x = age, ymin = func(fit-se.fit), ymax = func(fit+se.fit), color=cat, fill=cat), alpha = 0.3) +
  labs(y='s(age)',x=lookup[['age']],color='',fill='') + 
  scale_colour_discrete_phs(palette='all') +
  scale_fill_discrete_phs(palette='all') +
  theme_classic()

pred %>% filter(cat=='1 Dose AZ') %>% filter(age==50 | age==75 | age==100) %>% select(age,fit,se.fit)

categories <- unique(pred$cat)
data <- NULL
for (cat in categories){
  temp <- pred %>% filter(cat==cat)
  ref.fit <- (temp %>% filter(age==50))$fit
  ref.se.fit <- (temp %>% filter(age==50))$se.fit
  ref.ci_low <- ref.fit - 2*ref.se.fit
  ref.ci_high <- ref.fit + 2*ref.se.fit

  
  temp <- temp %>%  mutate(ci_low = fit-2*se.fit,ci_high=fit+2*se.fit ) %>%
          mutate(or=exp(fit-ref.fit),
                lor=exp(ci_low - ref.ci_low),
                uor=exp(ci_high - ref.ci_high))
  
  data <- rbind(data,temp)
}



ggplot(data%>% filter(cat=='1 Dose AZ'),aes(x=age, y=or, color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  #geom_line(aes(y=lor,x=age),size = 0.7) +
  geom_ribbon(aes(x = age, ymin = lor, ymax = uor, color=cat, fill=cat), alpha = 0.3) +
  labs(y='Odds Ratio',x=lookup[['age']],color='',fill='') + 
  scale_colour_discrete_phs(palette='all') +
  scale_fill_discrete_phs(palette='all') +
  geom_hline(yintercept=1,linetype='dashed') +
  geom_vline(xintercept=50,linetype='dashed') +
  #scale_y_log10() +
  theme_classic()



plot_gam(gamFit, pred = "days_since_first_measurement")

plot_df <- no_plot(gamFit)

set_pred <- grep(paste0("\\b", "age", "\\b"), plot_df)
set_pred

plot_df[[2]]

set_pred <- 6

df <- data.frame(
  x = plot_df[[set_pred]]$x,
  se_upr = plot_df[[set_pred]]$fit + plot_df[[set_pred]]$se,
  se_lwr = plot_df[[set_pred]]$fit - plot_df[[set_pred]]$se,
  y = plot_df[[set_pred]]$fit
)
df

col_line <- "blue"
ci_line_col <- "black"
ci_line_type <- "dashed"
ci_fill <- "grey"
ci_alpha <- 0.4
ci_line_size <- 0.8
sm_fun_size <- 1.1
title <- NULL
xlab <- NULL
ylab <- NULL

ggplot(df, aes_(~x, ~y)) +
  geom_line(colour = col_line, size = sm_fun_size) +
  geom_line(aes_(~x, ~se_upr),
            linetype = ci_line_type,
            colour = ci_line_col, size = ci_line_size
  ) +
  geom_line(aes_(~x, ~se_lwr),
            linetype = ci_line_type,
            colour = ci_line_col, size = ci_line_size
  ) +
  geom_ribbon(aes_(x = ~x, ymin = ~se_lwr, ymax = ~se_upr),
              fill = ci_fill, alpha = ci_alpha
  ) +
  ylab(ylab) +
  xlab(xlab)



plot_gam(gamFit, pred = "age")


formula <- as.formula("outcome ~ s(age,by=cat) + cat ") 
gamFit <- perform_gam(model,formula)
plot_smooths(
  model = gamFit,
  series = age,
  comparison = cat)+ theme_classic()



formula <- as.formula("outcome ~ s(age) + cat + Q_DIAG_BLOOD_CANCER") 
gamFit <- perform_gam(model,formula)
gamFit

values <- list(cat='Pf-Pf',Q_DIAG_BLOOD_CANCER=0)
temp <- tidymv::predict_gam(gamFit,values=values)
temp

temp %>%
  ggplot(aes(age, fit)) +
  geom_smooth_ci(cat)



plot_smooths(
  model = gamFit,
  series = age,
  facet_terms = Q_DIAG_BLOOD_CANCER,
  exclude_terms = " Q_BMI ",
  comparison = cat
  )+ theme_classic()

plot_gam(gamFit, pred = "age")

par(mfrow = c(2,2)); plot(gamFit)



agelims = model %>%
  select(age) %>%
  range
# Generate a sequence of age values spanning the range
age_grid = seq(from = min(agelims), to = max(agelims))


# Predict the value of the generated ages, 
# returning the standard error using se = TRUE
temp <- model

#formula <- as.formula("outcome ~ s(age)")
#gamFit <- perform_gam(model,formula)
pred = predict(gamFit, list(age = age_grid), se = TRUE, type='response')

# Compute error bands (2*SE)
se_bands = with(pred, cbind("upper" = fit+2*se.fit, 
                            "lower" = fit-2*se.fit))

pred
# Plot the spline and error bands
ggplot() +
  geom_line(aes(x = age_grid, y = pred$fit), color = "#0000FF") + 
  geom_ribbon(aes(x = age_grid, 
                  ymin = se_bands[,"lower"], 
                  ymax = se_bands[,"upper"]), 
              alpha = 0.3) +
  xlim(agelims)


summary(gamFit)


plot_difference(
  gamFit,
  series = age,
  difference = list(cat=c("Pf-Pf", "AZ-AZ"))
)



summary(gamFit)
plot(gamFit)

#plot_gam(gamFit, pred = "age")
summary(gamFit)

nrow(model)
model %>% group_by(stage) %>% summarise(n=n())





plot_object <- plot_gam(gamFit, pred = "days_since_last_vac")
plot_object + theme_classic() + labs(x='Days Since Last Vaccination Measurement',y='fx()',title='Pfizer or Moderna')



min(df_ana$Sampledate_iso)

plot_object <- plot_gam(gamFit, pred = "days_since_first_measurement")
or_object <- or_gam(
  data = model, model = gamFit, pred = "days_since_first_measurement",
  values = c(50, 200)
)

plot <- insert_or(plot_object, or_object,
                  or_yloc = 1.5+0.4*i,
                  line_col = 'black',
                  line_size = 0.3,
                  line_type = 'dashed',
                  arrow = TRUE,
                  arrow_length = 4,
                  arrow_col = 'orange',
                  arrow_yloc = 1.5+0.4*i,
                  values = TRUE
)
plot + theme_classic() + labs(x='Days Since The First Serology Measurement',y='s(days_since_first_measurement)')


values <- c(0,20,40,60,100)
plot <-  plot_gam(gamFit, pred = "age")
plot

for(i in 1:length(values)-1){

  
  v1 <- values[i]
  v2 <- values[i+1]
  
  
  or_object <- or_gam(
    data = model, model = gamFit, pred = "age",
    values = c(v1, v2)
  )
  
  plot <- insert_or(plot, or_object,
                    or_yloc = 1.5+0.4*i,
                    line_col = 'black',
                    line_size = 0.3,
                    line_type = 'dashed',
                    arrow = TRUE,
                    arrow_length = 4,
                    arrow_col = 'orange',
                    arrow_yloc = 1.5+0.4*i,
                    values = TRUE
  )
  
}

plot + theme_classic() + labs(x='Age',y='s(ageYear)')


formula <- paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + Q_BMI + days_since_last_vac + days_between_d1_d2 + product_binary")
tab <- perform_glm(modelA,formula)
plot_simple(tab)
tab



formula <- as.formula(paste0("outcome ~ s(days_since_first_measurement)"))# + Sex + n_risk_gps + product_binary + stage"))
gamFit <- perform_gam(df_ana,formula)
plot(gamFit)

formula <- as.formula(paste0("outcome ~ s(as.integer(days_between_d1_d2))"))# + Sex + n_risk_gps + product_binary + stage"))
gamFit <- perform_gam(df_ana %>% filter(stage!=1),formula)
plot(gamFit)

exp(coef(tab))

Vb <- vcov(tab, unconditional = TRUE)
se <- sqrt(diag(Vb))

confint(tab, parm = "age", type = "confidence")

confint(tab)
