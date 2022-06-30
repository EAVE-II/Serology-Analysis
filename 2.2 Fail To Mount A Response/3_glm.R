
library(phsstyles)
library(splines)
library(xtable)
require(gridExtra)

perform_glm_fit <- function(model,formula) {
  
  glmFit <- glm(formula,family=binomial,data=model)
  return (glmFit)
}

perform_glm <- function(model,formula,do_adjusted=TRUE,do_unadjusted=TRUE){
  formula <-  gsub("[\r\n]", " ", formula)
  
  tab <- NULL
  if (do_adjusted){
    glmFit <- perform_glm_fit(model,formula)
    
    print (with(summary(glmFit), 1 - deviance/null.deviance))
    print (glmFit)
    
    tab <- exp(cbind(coef(glmFit), confint(glmFit))) %>% as.data.frame
    tab <- tab[-1,]
    tab$names <- rownames(tab)
    
    
    names(tab)[1] <- 'OR'
    names(tab)[2] <- 'LCL'
    names(tab)[3] <- 'UCL'
  }
 
  
  if (do_unadjusted){
    vars <- strsplit(formula,split="~")[[1]][2]
    vars <- strsplit(unlist(gsub(" ","",vars)),split="+",fixed=TRUE)[[1]]
    tab_unadjusted <- NULL
    for (name in vars){
      formula <- paste0("outcome ~ ",name)
      glmFit <- glm(formula,family=binomial,data=model)
      glmFit.summary <- summary(glmFit)
      p = coef(glmFit.summary)[,"Pr(>|z|)"][[2]]
      print ( p)
      t <- exp(cbind(coef(glmFit), confint(glmFit))) %>% as.data.frame
      t <- t[-1,]
      
      names(t)[1] <- 'uOR'
      names(t)[2] <- 'uLCL'
      names(t)[3] <- 'uUCL'
      
      t$p <- p
      t<- t %>% mutate(p=ifelse(p<0.05,'<0.05',round(p,2)))
      
      
      n <- as.data.frame.matrix(table(model[[name]],model$outcome)) %>% 
           mutate(percentage=round(100*`1`/`0`,2)) %>% 
           rownames_to_column('rn') %>%
           mutate(rn=ifelse(rn==1,name,paste0(name,rn))) %>% 
           column_to_rownames('rn')
      
      print (nrow(n))
      print (n)
      
      
      t <- merge(n,t,by=0,all=T) %>%  column_to_rownames('Row.names')
      
      
      tab_unadjusted = rbind(tab_unadjusted,t)
      
    }
    
    if (do_adjusted) {
      tab <- cbind(tab,tab_unadjusted)
    }
    else {
      tab <- tab_unadjusted
      tab$names <- rownames(tab)
      
    }
  }
  
  print (tab)
  print (tab$names)
  tab <- tab %>% mutate( names=ifelse(is.na(as.character(lookup[names])),names,lookup[names]))
  rownames(tab) <- tab$names
  print (tab)
  
  name <- deparse(substitute(model))
  
  #xtable(tab)
  #xtable(anova(glmFit))
  

  #print(xtable(tab %>% select(-names), type = "latex"),file=paste0(results_folder,name,"_fit.tex"))
  
  #print(xtable(anova(glmFit), type = "latex"),file=paste0(results_folder,name,"_anova.tex"))
  
  return (tab)
}
  
plot_fit <- function(tab,name,xmin=0.01,xmax=50,adjusted=TRUE,unadjusted=TRUE){
  
  tab$names <- factor(tab$names,levels=tab$names)
  
  p<- ggplot(tab, aes(x=names))# reorder(names,OR))) 
    
  if(adjusted){
    p<- p + geom_pointrange(aes(y=OR, ymin=LCL, ymax=UCL, fill='Adjusted'), width=.2,alpha=0.7,shape=21,
                    position=position_dodge(.9)) 
  }
  if(unadjusted){
    p<- p + geom_pointrange(aes(y=uOR, ymin=uLCL, ymax=uUCL, fill='Unadjusted'), width=.2,alpha=0.6,shape=21,
                    position=position_dodge(.9))
    #if(!adjusted){
    #  p<- p + geom_label(aes(x=names,y=xmax,label=round(uOR,2)))
    #}
  }
  p<- p + 
    scale_fill_manual(values=c(phs_colours("phs-magenta"), phs_colours("phs-green"))) +
    labs(title='',x='',y='OR (95% CI)',fill='Fit Type') +
    geom_hline(yintercept=c(1), linetype="dashed") +
    geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
    #ylim(0., 10) +
    scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
    coord_flip(ylim=c(xmin,xmax)) +
    scale_y_log10() +
    theme_classic()
  ggsave(paste0(results_folder,name,".png"),height=7)
  return (p)
}

plot_simple <- function(tab,xmin=0.01,xmax=50){
  tab$names <- factor(tab$names,levels=tab$names)
  
  p<- ggplot(tab, aes(x=names))+# reorder(names,OR))) + 
        geom_pointrange(aes(y=OR, ymin=LCL, ymax=UCL, fill='Adjusted'), width=.2,alpha=0.7,shape=21,
                        position=position_dodge(.9)) +
        scale_fill_manual(values=c(phs_colours("phs-magenta"), phs_colours("phs-green"))) +
        ylim(0., 16.0) +
        labs(title='',x='',y='OR (95% CI)') +
        geom_hline(yintercept=c(1), linetype="dashed") +
        geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
        #ylim(0., 10) +
        scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
        coord_flip(ylim=c(xmin,xmax)) +
        theme_classic() + guides(fill="none")
  return (p)
}

plot_strat <- function(tab){
  tab$names <- factor(tab$names,levels=unique(tab$names))
  p<- ggplot(tab ,aes(x=names))+#reorder(names,OR))) + 
    geom_pointrange(aes(y=OR, ymin=LCL, ymax=UCL, fill=strat), width=.2,alpha=0.7,shape=21,
                    position=position_dodge(.9)) +
    labs(title='',x='',y='OR (95% CI)',fill='') +
    geom_hline(yintercept=c(1), linetype="dashed") +
    #geom_hline(yintercept=c(3) linetype="dotted", alpha=0.5) +
    geom_vline(xintercept=c(2:100)-0.5, linetype='dotted',alpha=0.8) +
    scale_y_log10() + #limits=c(0.1,50)) +
    coord_flip(ylim=c(0.02,200)) +
    theme_classic()
  return (p)
}

#modelA %>% group_by(outcome) %>% summarise(n=n())
#modelA.formula <- paste0("outcome ~ ageYear + Q_BMI + Sex + n_risk_gps  + product_binary + stage ")
#modelA.fit <- perform_glm_fit(modelA,modelA.formula)

modelA.formula <- as.formula(paste0("outcome ~ ageYear + Sex + simd2020v2_sc_quintile + n_risk_gps + product_binary"))
modelA.tab <- perform_glm(modelA %>% filter(stage==2),modelA.formula,do_unadjusted = FALSE)

modelA.fit <- perform_glm_fit(modelA,modelA.formula)

exp(cbind(coef(modelA.fit), confint(modelA.fit))) %>% as.data.frame

modelA.tab2 <- or_glm(data = modelA, model = modelA.fit, incr=list(product_binary=1))
modelA.tab2
tab <- modelA.tab2 %>% as
names(tab)[1] <- 'names'
names(tab)[2] <- 'OR'
names(tab)[3] <- 'LCL'
names(tab)[4] <- 'UCL'
tab
#tab <- tab %>% mutate( names=lookup[names])
tab
p <- plot_simple(tab)
p

modelA.tab
modelA.p <- plot_simple(modelA.tab)#,'modelA')#xmin=0.1,xmax=10)#,'modelA') 
modelA.p  + labs(title='Stage 2')

grid.arrange(p,modelA.p,pgam,nrow=1)


modelA.tab_s2 <- perform_glm(modelA %>% filter(stage==2),modelA.formula,do_unadjusted = TRUE)
modelA.tab_s2 %>% drop_na()
modelA.p_s2 <- plot_fit(modelA.tab_s2 %>% drop_na(),'modelA')#xmin=0.1,xmax=10)#,'modelA') 
modelA.p_s2

modelA.tab_s3 <- perform_glm(modelA %>%  mutate(cat = relevel(cat,"Pf-Pf-Pf"))  %>% filter(stage>2),modelA.formula,do_unadjusted = TRUE)
modelA.tab_s3 %>% drop_na()
modelA.p_s3 <- plot_fit(modelA.tab_s3 %>% drop_na(),'modelA',xmin=0.1,xmax=300)#,'modelA') 
modelA.p_s3

modelA.tab_s1 <- perform_glm(modelA %>%  mutate(cat = relevel(cat,"Pf"))  %>% filter(stage==1),modelA.formula,do_unadjusted = TRUE)
modelA.tab_s1 %>% drop_na()
modelA.p_s1 <- plot_fit(modelA.tab_s1 %>% drop_na(),'modelA')#xmin=0.1,xmax=10)#,'modelA') 
modelA.p_s1

grid.arrange(modelA.p_s1,modelA.p_s2,modelA.p_s3,ncol=3)



modelC %>% group_by(d1_product) %>% summarise(n=n())

modelC.tab <- perform_glm(modelC,modelC.formula)
modelC.tab

modelC.p <- plot_fit(modelC.tab,'modelC')
modelC.p


df_ana %>% select(outcome,age)


modelB.tab <- perform_glm(modelB,modelB.formula,do_unadjusted = TRUE)
modelB.tab

modelB.p <- plot_simple(modelB.tab)
modelB.p + labs(title='Adjusted Odds')

modelB.p <- plot_fit(modelB.tab,'modelB')
modelB.p



f <-  "outcome ~ ageYear + Sex + Q_DIAG_BLOOD_CANCER  + product_binary + stage"
modelB.tab <- perform_glm(modelB,f,do_unadjusted = FALSE)
modelB.p <- plot_simple(modelB.tab)
modelB.p + labs(title='Adjusted Odds')

f <-  "outcome ~ ageYear + Sex + Q_DIAG_BLOOD_CANCER  + product_binary + stage*product_binary"
modelB.tab_int <- perform_glm(modelB,f,do_unadjusted = FALSE)
modelB.tab
modelB.tab_int

modelB.p <- plot_simple(modelB.tab_int)
modelB.p + labs(title='Adjusted Odds')


f <-  "outcome ~ Q_DIAG_BLOOD_CANCER  + product_binary "
modelB.tab <- perform_glm(modelB,f,do_unadjusted = FALSE)
modelB.p <- plot_simple(modelB.tab)
modelB.p + labs(title='Adjusted Odds')

f <-  "outcome ~  product_binary*Q_DIAG_BLOOD_CANCER"
modelB.tab_int <- perform_glm(modelB,f,do_unadjusted = FALSE)
modelB.tab
modelB.tab_int

modelB.p <- plot_simple(modelB.tab_int)
modelB.p + labs(title='Adjusted Odds')

modelB %>% filter(Q_DIAG_BLOOD_CANCER==1) %>% group_by(product_binary,outcome) %>% summarise(n=n())


modelB.tab_az <- perform_glm(modelB %>% filter(product_binary==1),modelB.formula2)
modelB.tab_az
modelB.p <- plot_simple(modelB.tab_az)
modelB.p + labs(title='Last Vaccine AZ')

f <- gsub("\\+ Q_DIAG_CONGEN_HD ","",modelB.formula2)
modelB.tab_other <- perform_glm(modelB %>% filter(product_binary==0), gsub("\\+ Q_DIAG_DEMENTIA ","",f))
modelB.p <- plot_simple(modelB.tab_other)
modelB.p + labs(title='Last Vaccine Pfizer')

modelB.tab_other['strat'] = 'Pfizer or Moderna'
modelB.tab_other
modelB.tab_az['strat'] = 'AstraZeneca'
modelB.tab_az

temp <- modelB.tab_az %>% select(names,OR,LCL,UCL,strat) %>% rbind(modelB.tab_other %>% select(names,OR,LCL,UCL,strat) %>% drop_na() )
temp

plot_strat(temp)


nrow(modelB)
nrow(modelB %>% filter(stage==1))
nrow(modelB %>% filter(stage==2))
nrow(modelB %>% filter(stage=='2+'))

f <- modelB.formula2
f
f <- "outcome ~ ageYear + Sex + Q_BMI+  Q_DIAG_CKD_LEVEL + Q_DIAG_AF +
      Q_DIAG_BLOOD_CANCER + Q_DIAG_CCF + Q_DIAG_CHD + Q_DIAG_CIRRHOSIS + 
      Q_DIAG_COPD + Q_DIAG_DEMENTIA + Q_DIAG_DIABETES_1 + Q_DIAG_FRACTURE +
      Q_DIAG_NEURO + Q_DIAG_PULM_RARE + Q_DIAG_PVD + Q_DIAG_RESP_CANCER +
      Q_DIAG_STROKE + Q_DIAG_VTE + prior_infection  + stage + product_binary"
modelB.tab_d1 <- perform_glm(modelB %>% filter(stage==1),gsub("\\+ stage","",f))
modelB.tab_d2 <- perform_glm(modelB %>% filter(stage==2),gsub("\\+ stage","",f))
modelB.tab_d3 <- perform_glm(modelB %>% filter(stage=='2+'),gsub("\\+ stage","",f))

modelB.tab_d1['strat']  = '1 Dose'
modelB.tab_d2['strat'] = '2 Doses'
modelB.tab_d3['strat'] = '3+ Doses'
temp <- modelB.tab_d1 %>% select(names,OR,LCL,UCL,strat) %>% rbind(modelB.tab_d2 %>% select(names,OR,LCL,UCL,strat) ) %>%
        rbind(modelB.tab_d3 %>% select(names,OR,LCL,UCL,strat) %>% filter(LCL>0.01 & UCL<100) %>% drop_na()) 
plot_strat(temp)

p1 <- plot_simple(modelB.tab_d1)
p1 + labs(title='Dose 1')
p2 <- plot_simple(modelB.tab_d2)
p2 + labs(title='Dose 2')
p3 <- plot_simple(modelB.tab_d3)
p3 + labs(title='Dose 3+')


f <- "outcome ~ ageYear + Sex + Q_DIAG_CKD_LEVEL + Q_DIAG_AF +
      Q_DIAG_BLOOD_CANCER + Q_DIAG_CCF + Q_DIAG_CHD + Q_DIAG_CIRRHOSIS + 
      Q_DIAG_COPD + Q_DIAG_DEMENTIA + Q_DIAG_DIABETES_1 + Q_DIAG_FRACTURE +
      Q_DIAG_NEURO + Q_DIAG_PULM_RARE + Q_DIAG_PVD + Q_DIAG_RESP_CANCER +
      Q_DIAG_STROKE + Q_DIAG_VTE + prior_infection  + stage "
modelB.tab_d1_az <- perform_glm(modelB %>% filter(stage==1 & product_binary==1),gsub("\\+ stage","",f))
modelB.tab_d2_az <- perform_glm(modelB %>% filter(stage==2 & product_binary==1),gsub("\\+ stage","",f))
modelB.tab_d3_az <- perform_glm(modelB %>% filter(stage=='2+' & product_binary==1),gsub("\\+ stage","",f))

modelB.tab_d1_other <- perform_glm(modelB %>% filter(stage==1 & product_binary==0),gsub("\\+ stage","",f))
modelB.tab_d2_other <- perform_glm(modelB %>% filter(stage==2 & product_binary==0),gsub("\\+ stage","",f))
modelB.tab_d3_other <- perform_glm(modelB %>% filter(stage=='2+' & product_binary==0),gsub("\\+ stage","",f))


modelB.tab_d1_az['strat']  = 'AstraZeneca'
modelB.tab_d2_az['strat'] = 'AstraZeneca'
modelB.tab_d3_az['strat'] = 'AstraZeneca'

modelB.tab_d1_other['strat']  = 'Pfizer or Moderna'
modelB.tab_d2_other['strat'] = 'Pfizer or Moderna'
modelB.tab_d3_other['strat'] = 'Pfizer or Moderna'

modelB.tab_d1_other %>% select(names,OR)
temp <- modelB.tab_d1_az %>% select(names,OR,LCL,UCL,strat) %>% rbind(modelB.tab_d1_other %>% select(names,OR,LCL,UCL,strat) )
plot_strat(temp)



grid.arrange(p1,p2)




modelB.p <- plot_fit(modelB.tab_other%>% drop_na(),'modelB_other')
modelB.p 


modelA.formula <- paste0("outcome ~ ageYear + Q_BMI + Sex + n_risk_gps  + product_binary + stage ")
modelA.fit <- perform_glm_fit(modelA,modelA.formula)
modelA.tab <- perform_glm(modelA,modelA.formula,do_unadjusted = TRUE)
modelA.p <- plot_fit(modelA.tab,'modelA')#xmin=0.1,xmax=10)#,'modelA') 
modelA.p


modelA.formula <- paste0("outcome ~ ageYear + Sex  + Q_BMI + n_risk_gps*stage + product_binary  ")
modelA.fit_int <- perform_glm_fit(modelA,modelA.formula)
modelA.tab_int <- perform_glm(modelA,modelA.formula,do_unadjusted = FALSE)
modelA.p_int <- plot_simple(modelA.tab_int,xmin=0.1,xmax=10)#,'modelA') 
modelA.p_int

grid.arrange(modelA.p,modelA.p_int,nrow=1)


temp <- strsplit(gsub("\\+ product_binary","",modelA.formula),"\\+")
for (i in 1:lengths(temp)){
  f <- temp[[1]]
  f[[i]] <- paste0(f[[i]],"* product_binary ")
  f <- paste0(f,collapse="+")
  print (f)
  modelA.fit_int <- perform_glm_fit(modelA,f)
  comp <- anova(modelA.fit,modelA.fit_int,test='Chisq')
  p <- comp$"Pr(>Chi)"[[2]]
  print (paste0("Pr(>Chi) = ",p))
  
}

f <- "outcome ~ ageYear * product_binary + Sex + n_risk_gps + Q_BMI  + stage * product_binary"
modelA.fit_int <- perform_glm_fit(modelA,f)
comp <- anova(modelA.fit,modelA.fit_int,test='Chisq')
p <- comp$"Pr(>Chi)"[[2]]
print (paste0("Pr(>Chi) = ",p))
modelA.tab_int <- perform_glm(modelA,f,do_unadjusted = FALSE)
modelA.p_int <- plot_simple(modelA.tab_int,xmin=0.1,xmax=10)#,'modelA') 
modelA.p_int
grid.arrange(modelA.p,modelA.p_int,nrow=1)

comp <- anova(modelA.fit,modelA.fit_int,test='Chisq')
p <- comp$"Pr(>Chi)"[[2]]
p

modelA.formula <- paste0("outcome ~ ageYear + Sex + n_risk_gps + Q_BMI + stage*product_binary  ")
modelA.fit_int <- perform_glm_fit(modelA,modelA.formula)
comp <- anova(modelA.fit,modelA.fit_int,test='Chisq')
p <- comp$"Pr(>Chi)"[[2]]
p

modelA.formula3 <- paste0("outcome ~ ageYear")

modelA.formula4 <- paste0("outcome ~ ns(ageYear,df=5)")

modelA.formula2 <- paste0("outcome ~ pspline(ageYear,df=0,intercept=FALSE) + Sex + n_risk_gps + ",paste0(modelA.vars,collapse = " + "))
modelA.formula2 
modelA.p2 <- perform_glm(modelA,modelA.formula2)
modelA.p2




modelB.p2 <- perform_glm(modelB,modelB.formula2)
modelB.p2

modelC_1.p <- perform_glm(modelC_1,modelC_1.formula)
modelC_1.p

modelC_2.p <- perform_glm(modelC_2,modelC_2.formula)
modelC_2.p
ggsave(paste0(results_folder,'modelC_2',".png"),height=6)


modelA_2.p <- perform_glm(modelA_2,modelA_2.formula)
modelA_2.p

modelA_3.p <- perform_glm(modelA_3,modelA_3.formula)
modelA_3.p
