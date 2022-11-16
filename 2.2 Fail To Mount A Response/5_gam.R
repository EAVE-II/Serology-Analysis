
library(mgcv)
library("oddsratio")
library(ggplot2)
library(tidymv)
library(gratia)
library(purrr)
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
      print (value)
      print (tor)
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


get_or_from_gam <- function(mod, level = 0.95) {
  # a method for extracting confidence intervals and returning a tidy data frame
  require(mgcv)
  require(dplyr)
  
  mod.s <- summary(mod)
  
 
  E <- data.frame(OR = mod.s$p.coeff, pv = mod.s$p.pv) %>%
    mutate(names = row.names(.)) %>%
    select(names, OR, pv)
  
  SE <- data.frame(se= mod.s$se) %>%
    mutate(names = row.names(.)) %>%
    select(names, se)
  
  nu <- mod.s$residual.df
  
  
  inner_join(E, SE) %>%
    mutate(
           original_name = names,
           names = ifelse(is.na(lookup[names]),names,lookup[names]),
           LCL = OR +
             se * qt(df = nu,
                     p = (1 - level) / 2),
           UCL = OR +
             se * qt(df = nu,
                     p = 1 - (1 - level) / 2)) %>% slice(-1) %>% mutate_at(vars(-names,-original_name,-pv),exp) %>% 
    return
  
}


fit_gam <- function(gamFit,var,cats = NULL, value = NULL) {
  term_list <- list()
  for (term in labels(gamFit$terms)){
    new_term <- gamFit[["var.summary"]][[term]][[1]]
    term_list <- append(term_list, list(new_term))
  }
  
  names(term_list) <-  labels(gamFit$terms)
  
  if (var =='age'){
    term_list[['age']] <- seq(0,80,1)
  }
  else if(var =='days_since_first_measurement'){
    term_list[['days_since_first_measurement']] <- seq(0,450,1)
  }
  else if(var =='days_since_last_vac'){
    term_list[['days_since_last_vac']] <- seq(0,300,1)
  }
  
  if (!is.null(cats)){
    for (lab in labels(cats)){
      term_list[[lab]] <- cats[[lab]]
    }
  }
  if(!is.null(value)){
    term_list[[var]] <- value
  }
  
  
  
  new_data <- expand.grid(term_list)
  
  pred <- predict.gam(gamFit,new_data,se.fit = TRUE)
  pred <- cbind(new_data, pred)
  pred <- pred %>% mutate(cat=lookup[paste0('cat',cat)])
  pred <- pred %>% mutate(stage=ifelse(grepl('1',cat),'Dose 1',ifelse(grepl('2',cat),'Dose 2','Dose 3')))
  

  return (pred);
}

plot_gam <- function(gamFit,var, var_ref=0, cats = NULL) {
 
 
  palette <- 'none'
  if (!is.null(cats)){
    palette <- 'all' 
  }
  
  pred <- fit_gam(gamFit,var,cats)
  
  
  ref <- pred %>% filter(!!as.name(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
  
  
  pred <- pred %>% left_join(ref, by='cat')
  
  
  
  func <- function(a,b){
    return (exp(a - b));
  }
  
  
  p <- ggplot(pred, aes(x=!!as.name(var), y=func(fit,ref_fit), color=cat)) +
    geom_line(size = 0.5, linetype='dashed') +
    geom_ribbon(aes(x = !!as.name(var), ymin = func(fit-se.fit,ref_fit), ymax = func(fit+se.fit,ref_fit), color=cat, fill=cat), alpha = 0.3) +
    labs(y='Odds Ratio',x=var,color='',fill='') + 
    scale_colour_discrete_phs(palette=palette) +
    geom_vline(xintercept=var_ref,linetype='dotted') +
    geom_hline(yintercept=1,linetype='dotted') +
    scale_fill_discrete_phs(palette=palette) +
    #scale_y_log10() +
    theme_classic()
  
  if(palette=='none'){
    print ('here')
    p <- p + theme(legend.position="none")
  }
  return(p)
}

add_names <- function(df,formula){
  labs <- labels(terms(as.formula(formula)))
  temp <- df %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
    unnest(original) %>%
    filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 
  return (df %>% left_join(temp))
}


# primary care model A
pc.modelA.gamFit <- perform_gam(pc.modelA,as.formula(pc.modelA.gam.formula))

pc.modelA.gamFit.or <- get_or_from_gam(pc.modelA.gamFit) %>% add_names(pc.modelA.gam.formula)
colnames(pc.modelA.gamFit.or)

pc.modelA.gamFit.or %>% plot_simple(labels=F) 



#sensitivity analysis Loose and Tight
pc.modelA.T.gamFit <- perform_gam(pc.modelA.T,as.formula(pc.modelA.T.gam.formula))
pc.modelA.T.gamFit.or <- get_or_from_gam(pc.modelA.T.gamFit) %>% add_names(pc.modelA.T.gam.formula)
pc.modelA.T.gamFit.or %>% plot_simple(labels=F) 

pc.modelA.L.gamFit <- perform_gam(pc.modelA.L,as.formula(pc.modelA.L.gam.formula))
pc.modelA.L.gamFit.or <- get_or_from_gam(pc.modelA.L.gamFit) %>% add_names(pc.modelA.L.gam.formula)
pc.modelA.L.gamFit.or %>% plot_simple(labels=F) 



temp <- pc.modelA.gamFit.or %>% mutate(ftype='Nominal') %>% 
  rbind( pc.modelA.T.gamFit.or %>% mutate(ftype='Tight') ) %>% 
  rbind( pc.modelA.L.gamFit.or %>% mutate(ftype='Loose') ) %>% mutate(ftype=as.factor(ftype))


levels <- (pc.modelA.L.gamFit.or %>% select(original) %>% unique %>% map_dfr(rev))$original
p <- temp %>% ggplot(aes(x=as.factor(names)))+
  geom_pointrange(aes(y=OR, ymin=LCL, ymax=UCL, fill=ftype), width=.2,alpha=0.7,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-green"),phs_colours("phs-magenta"),phs_colours("phs-rust"))) +
  #scale_fill_discrete_phs(palette='main') +
  ylim(0., 16.0) +
  labs(title='',x='',y='OR (95% CI)',fill='f-type') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  ##ylim(0., 10) +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  coord_flip(ylim=c(0.05,50)) +
  theme_classic() + guides() + #fill="none") 
  facet_grid(factor(original,levels=levels) ~ ftype,scales="free", space="free") +                   
  theme(strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'none',
        legend.justification = 'center') +
  theme(panel.background = element_rect(fill = NA, color = "black")) 
p
 
ggsave(paste0(folder,"modelA_pc_sensitivity_ORs.pdf"), p , width=7, height=9, dpi=300, units="in")




# blood donors model A

bd.modelA.gamFit <- perform_gam(bd.modelA,as.formula(bd.modelA.gam.formula))
bd.modelA.gamFit.or <- get_or_from_gam(bd.modelA.gamFit) %>% add_names(bd.modelA.gam.formula)
#  mutate_at(vars(-names,-original_name),~round(.,2)) %>%
#  mutate(pv = paste0('p ',ifelse(pv<0.05,'<0.05',paste0('= ',sprintf('%.2f',pv)))))

#labs <- labels(terms(as.formula(bd.modelA.gam.formula)))
#temp <- bd.modelA.gamFit.or %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
#  unnest(original) %>%
#  filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 

#bd.modelA.gamFit.or <- bd.modelA.gamFit.or %>% left_join(temp)

bd.modelA.gamFit.or %>% plot_simple(labels=F)




#sensitivity analysis Loose and Tight
bd.modelA.T.gamFit <- perform_gam(bd.modelA.T,as.formula(bd.modelA.T.gam.formula))
bd.modelA.T.gamFit.or <- get_or_from_gam(bd.modelA.T.gamFit) %>% add_names(bd.modelA.T.gam.formula)

bd.modelA.L.gamFit <- perform_gam(bd.modelA.L,as.formula(bd.modelA.L.gam.formula))
bd.modelA.L.gamFit.or <- get_or_from_gam(bd.modelA.L.gamFit) %>% add_names(bd.modelA.L.gam.formula)


bd.modelA.T.gamFit.or %>% filter(se<1000)



temp <- bd.modelA.gamFit.or %>% mutate(ftype='Nominal') %>% 
  rbind( bd.modelA.T.gamFit.or  %>% filter(se<15) %>% mutate(ftype='Tight') ) %>% 
  rbind( bd.modelA.L.gamFit.or %>% mutate(ftype='Loose') ) %>% mutate(ftype=as.factor(ftype))


levels <- (pc.modelA.L.gamFit.or %>% select(original) %>% unique %>% map_dfr(rev))$original
p <- temp %>% ggplot(aes(x=as.factor(names)))+
  geom_pointrange(aes(y=OR, ymin=LCL, ymax=UCL, fill=ftype), width=.2,alpha=0.7,shape=21,
                  position=position_dodge(.9)) +
  scale_fill_manual(values=c(phs_colours("phs-green"),phs_colours("phs-magenta"),phs_colours("phs-rust"))) +
  #scale_fill_discrete_phs(palette='main') +
  ylim(0., 16.0) +
  labs(title='',x='',y='OR (95% CI)',fill='f-type') +
  geom_hline(yintercept=c(1), linetype="dashed") +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  ##ylim(0., 10) +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  coord_flip(ylim=c(0.05,50)) +
  theme_classic() + guides() + #fill="none") 
  facet_grid(factor(original,levels=levels) ~ ftype,scales="free", space="free") +                   
  theme(strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        legend.position = 'none',
        legend.justification = 'center') +
  theme(panel.background = element_rect(fill = NA, color = "black")) 
p

ggsave(paste0(folder,"modelA_bd_sensitivity_ORs.pdf"), p , width=7, height=9, dpi=300, units="in")









levels <- (pc.modelA.gamFit.or %>% select(original) %>% unique %>% map_dfr(rev))$original
temp <- pc.modelA.gamFit.or %>% mutate(strat='Primary Care') %>% rbind(bd.modelA.gamFit.or %>% mutate(strat='Blood Donors')) 


p.modelA.ORs <- plot_strat(temp) + facet_grid( factor(original,levels=levels) ~ fct_rev(strat) ,scales="free", space="free") +  
                   theme(strip.background.x = element_blank(),
                         strip.text.y = element_blank(),
                         #legend.position = 'top',
                         legend.position = 'none',
                         legend.justification = 'center') +
  theme(panel.background = element_rect(fill = NA, color = "black")) 

ggsave(paste0(folder,"modelA_ORs.pdf"), p.modelA.ORs, width=7, height=9, dpi=300, units="in")




p.pc.modelA.gam.age <- plot_gam(pc.modelA.gamFit,'age',var_ref=40,list(cat=levels(pc.modelA$cat))) + 
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")
ggsave(paste0(folder,"modelA_pc_gam_age.pdf"), p.pc.modelA.gam.age, width=12, height=5.5, dpi=300, units="in")




p.bd.modelA.gam.age <- plot_gam(bd.modelA.gamFit,'age',var_ref=40,list(cat=levels(bd.modelA$cat))) + 
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

ggsave(paste0(folder,"modelA_bd_gam_age.pdf"), p.bd.modelA.gam.age, width=12, height=5.5, dpi=300, units="in")



var <- 'days_since_first_measurement'
var_ref <- 0

pred1 <- fit_gam(pc.modelA.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit)  %>% mutate(cat='Primary Care')
ref <- pred1 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred1 <- pred1 %>% left_join(ref,by='cat')


pred2 <- fit_gam(bd.modelA.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit) %>% mutate(cat='Blood Donors')
ref <- pred2 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred2 <- pred2 %>% left_join(ref,by='cat')


pred <- pred2 %>% rbind(pred1)

palette <- 'all' 

p.modelA.gam.days_since_first_measurement <- pred %>%
  mutate(date=as.Date(!!sym(var),origin='2020-12-20')) %>% 
  ggplot(aes(x=date, y=exp(fit - ref_fit), color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(x = date, ymin = exp(fit-se.fit - ref_fit ), ymax = exp(fit+se.fit - ref_fit), color=cat, fill=cat), alpha = 0.3) +
  labs(y='Odds Ratio',x='Date',color='',fill='') + 
  scale_colour_discrete_phs(palette=palette) +
  geom_vline(xintercept=var_ref,linetype='dotted') +
  geom_hline(yintercept=1,linetype='dotted') +
  scale_fill_discrete_phs(palette=palette) +
  #scale_y_log10() +
  theme_classic()  +
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

p.modelA.gam.days_since_first_measurement

ggsave(paste0(folder,"modelA_gam_days_since_first_measurement.pdf"), p.modelA.gam.days_since_first_measurement, width=6, height=9, dpi=300, units="in")



var <- 'days_since_last_vac'
var_ref <- 100

pred1 <- fit_gam(pc.modelA.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit)  %>% mutate(cat='Primary Care')
ref <- pred1 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred1 <- pred1 %>% left_join(ref,by='cat')


pred2 <- fit_gam(bd.modelA.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit) %>% mutate(cat='Blood Donors')
ref <- pred2 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred2 <- pred2 %>% left_join(ref,by='cat')


pred <- pred2 %>% rbind(pred1)

palette <- 'all' 

p.modelA.gam.days_since_last_vac <- ggplot(pred, aes(x=!!sym(var), y=exp(fit - ref_fit), color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(x = !!sym(var), ymin = exp(fit-se.fit - ref_fit ), ymax = exp(fit+se.fit - ref_fit), color=cat, fill=cat), alpha = 0.3) +
  labs(y='Odds Ratio',x=lookup[var],color='',fill='') + 
  scale_colour_discrete_phs(palette=palette) +
  geom_vline(xintercept=var_ref,linetype='dotted') +
  geom_hline(yintercept=1,linetype='dotted') +
  scale_fill_discrete_phs(palette=palette) +
  #scale_y_log10() +
  xlim(15,200) +
  theme_classic()  +
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

ggsave(paste0(folder,"modelA_gam_days_since_last_vac.pdf"), p.modelA.gam.days_since_last_vac, width=6, height=9, dpi=300, units="in")





# primary care model B
pc.modelB.gamFit <- perform_gam(pc.modelB,as.formula(pc.modelB.gam.formula))


pc.modelB.gamFit.or <- get_or_from_gam(pc.modelB.gamFit) %>% 
  mutate_at(vars(-names,-original_name),~round(.,2)) %>%
  mutate(pv_orig = pv, pv = paste0('p ',ifelse(pv<0.05,'<0.05',paste0('= ',sprintf('%.2f',pv)))))



labs <- labels(terms(as.formula(pc.modelB.gam.formula)))
temp <- pc.modelB.gamFit.or %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
  unnest(original) %>%
  filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 

pc.modelB.gamFit.or <- pc.modelB.gamFit.or %>% left_join(temp)
#pc.modelB.gamFit.or 

#labs <- unique(pc.modelB.gamFit.or$original)
#pc.modelB.gamFit.uor <- NULL
#for (lab in labs){
#  f <-  perform_gam(pc.modelB,as.formula(paste0('outcome ~ ',lab))) 
#  temp <- get_or_from_gam(f)# %>% select(original_name,OR,LCL,UCL) %>% rename(uOR=OR,uLCL=LCL,uUCL=UCL)
#  pc.modelB.gamFit.uor <- rbind(pc.modelB.gamFit.uor,temp)
#}
#labs <- labels(terms(as.formula(pc.modelB.gam.formula)))
#temp <- pc.modelB.gamFit.uor %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
#  unnest(original) %>%
#  filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 

#pc.modelB.gamFit.uor <- pc.modelB.gamFit.uor %>% left_join(temp)

#pc.modelB.gamFit.uor <- pc.modelB.gamFit.uor %>% mutate(strat2='unadjusted')
#pc.modelB.gamFit.or  <- pc.modelB.gamFit.or %>% mutate(strat2='adjusted')
#pc.modelB.gamFit.or 



#pc.modelB.gamFit.or %>% plot_simple(labels=F,xmax=20000) + 
#  geom_text(aes(x=names,y=150,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
#            hjust=0)

# blood donors model B

bd.modelB.gamFit <- perform_gam(bd.modelB,as.formula(bd.modelB.gam.formula))

bd.modelB.gamFit.or <- get_or_from_gam(bd.modelB.gamFit) %>% 
  mutate_at(vars(-names,-original_name),~round(.,2)) %>%
  mutate(pv_orig = pv, pv = paste0('p ',ifelse(pv<0.05,'<0.05',paste0('= ',sprintf('%.2f',pv)))))

labs <- labels(terms(as.formula(bd.modelB.gam.formula)))
temp <- bd.modelB.gamFit.or %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
  unnest(original) %>%
  filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 

bd.modelB.gamFit.or <- bd.modelB.gamFit.or %>% left_join(temp)

#labs <- unique(bd.modelB.gamFit.or$original)
#bd.modelB.gamFit.uor <- NULL
#for (lab in labs){
#  f <-  perform_gam(bd.modelB,as.formula(paste0('outcome ~ ',lab))) 
#  temp <- get_or_from_gam(f)# %>% select(original_name,OR,LCL,UCL) %>% rename(uOR=OR,uLCL=LCL,uUCL=UCL)
#  bd.modelB.gamFit.uor <- rbind(bd.modelB.gamFit.uor,temp)
#}
#labs <- labels(terms(as.formula(bd.modelB.gam.formula)))
#temp <- bd.modelB.gamFit.uor %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
#  unnest(original) %>%
#  filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 
#bd.modelB.gamFit.uor <- bd.modelB.gamFit.uor %>% left_join(temp)


#bd.modelB.gamFit.uor <- bd.modelB.gamFit.uor %>% mutate(strat2='unadjusted')
#bd.modelB.gamFit.or  <- bd.modelB.gamFit.or %>% mutate(strat2='adjusted')



#bd.modelB.gamFit.or %>% plot_simple(labels=F,xmax=20000) + 
#  geom_text(aes(x=names,y=150,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
#            hjust=0)


#levels <- (bd.modelB.gamFit.or %>% select(original) %>% unique %>% map_df(rev))$original
#pc.temp <- pc.modelB.gamFit.uor %>% mutate(strat2='Unadjusted') %>% rbind(pc.modelB.gamFit.or %>% mutate(strat2='Adjusted')) 
#bd.temp <- bd.modelB.gamFit.uor %>% mutate(strat2='Unadjusted') %>% rbind(bd.modelB.gamFit.or %>% mutate(strat2='Adjusted')) 
#temp <-  pc.temp %>% mutate(strat='Primary Care') %>% rbind(bd.temp %>% mutate(strat='Blood Donors')) 
#colnames(pc.modelB.gamFit.or)
#colnames(bd.modelB.gamFit.or)

temp <-  pc.modelB.gamFit.or %>% mutate(strat='Primary Care') %>% rbind(bd.modelB.gamFit.or %>% mutate(strat='Blood Donors')) 
temp <- temp %>% mutate(original = ifelse(grepl('Q_DIAG',original),'Q_DIAG',original))
temp


levels <- (pc.modelB.gamFit.or %>% select(original) %>% unique %>% map_dfr(rev))$original
p.modelB.ORs <- plot_strat(temp) + facet_grid( factor(original,levels=levels) ~ fct_rev(strat) ,scales="free", space="free") +  
  theme(strip.background.x = element_blank(),
        strip.text.y = element_blank(),
        #legend.position = 'top',
        legend.position = 'none',
        legend.justification = 'center') +
  theme(panel.background = element_rect(fill = NA, color = "black")) 

p.modelB.ORs



ggsave(paste0(folder,"modelB_ORs.pdf"), p.modelB.ORs, width=7, height=12, dpi=300, units="in")




p.pc.modelB.gam.age <- plot_gam(pc.modelB.gamFit,'age',var_ref=40,list(cat=levels(pc.modelB$cat))) + 
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

ggsave(paste0(folder,"modelB_pc_gam_age.pdf"), p.pc.modelB.gam.age, width=12, height=5.5, dpi=300, units="in")


p.bd.modelB.gam.age <- plot_gam(bd.modelB.gamFit,'age',var_ref=40,list(cat=levels(bd.modelB$cat))) + 
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

ggsave(paste0(folder,"modelB_bd_gam_age.pdf"), p.bd.modelB.gam.age, width=12, height=5.5, dpi=300, units="in")





var <- 'days_since_first_measurement'
var_ref <- 0

pred1 <- fit_gam(pc.modelB.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit)  %>% mutate(cat='Primary Care')
ref <- pred1 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred1 <- pred1 %>% left_join(ref,by='cat')


pred2 <- fit_gam(bd.modelB.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit) %>% mutate(cat='Blood Donors')
ref <- pred2 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred2 <- pred2 %>% left_join(ref,by='cat')


pred <- pred2 %>% rbind(pred1)

palette <- 'all' 

p.modelB.gam.days_since_first_measurement <- pred %>%
  mutate(date=as.Date(!!sym(var),origin='2020-12-20')) %>% 
  ggplot(aes(x=date, y=exp(fit - ref_fit), color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(x = date, ymin = exp(fit-se.fit - ref_fit ), ymax = exp(fit+se.fit - ref_fit), color=cat, fill=cat), alpha = 0.3) +
  labs(y='Odds Ratio',x='Date',color='',fill='') + 
  scale_colour_discrete_phs(palette=palette) +
  geom_vline(xintercept=var_ref,linetype='dotted') +
  geom_hline(yintercept=1,linetype='dotted') +
  scale_fill_discrete_phs(palette=palette) +
  #scale_y_log10() +
  theme_classic()  +
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

p.modelB.gam.days_since_first_measurement

ggsave(paste0(folder,"modelB_gam_days_since_first_measurement.pdf"), p.modelB.gam.days_since_first_measurement, width=6, height=9, dpi=300, units="in")



var <- 'days_since_last_vac'
var_ref <- 100

pred1 <- fit_gam(pc.modelB.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit)  %>% mutate(cat='Primary Care')
ref <- pred1 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred1 <- pred1 %>% left_join(ref,by='cat')


pred2 <- fit_gam(bd.modelB.gamFit,var,list(cat='Pf-Pf')) %>% select(!!sym(var),fit,se.fit) %>% mutate(cat='Blood Donors')
ref <- pred2 %>% filter(!!sym(var)==var_ref) %>% select(fit,cat) %>% rename(ref_fit = fit)
pred2 <- pred2 %>% left_join(ref,by='cat')


pred <- pred2 %>% rbind(pred1)

palette <- 'all' 

p.modelB.gam.days_since_last_vac <- ggplot(pred, aes(x=!!sym(var), y=exp(fit - ref_fit), color=cat)) +
  geom_line(size = 0.5, linetype='dashed') +
  geom_ribbon(aes(x = !!sym(var), ymin = exp(fit-se.fit - ref_fit ), ymax = exp(fit+se.fit - ref_fit), color=cat, fill=cat), alpha = 0.3) +
  labs(y='Odds Ratio',x=lookup[var],color='',fill='') + 
  scale_colour_discrete_phs(palette=palette) +
  geom_vline(xintercept=var_ref,linetype='dotted') +
  geom_hline(yintercept=1,linetype='dotted') +
  scale_fill_discrete_phs(palette=palette) +
  #scale_y_log10() +
  xlim(15,200) +
  theme_classic()  +
  facet_grid(~cat,scales="free") +  
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="top")

p.modelB.gam.days_since_last_vac 
ggsave(paste0(folder,"modelB_gam_days_since_last_vac.pdf"), p.modelB.gam.days_since_last_vac, width=6, height=9, dpi=300, units="in")







#grid.arrange(p,p2,ncol=2)



plot_my_gam <- function(gamFit,predictor,cat,ref,value){
  
  r <- fit_gam(gamFit,predictor,cat,ref) %>% select(!!as.name(predictor),cat,fit,stage)
  v <- fit_gam(gamFit,predictor,cat,value) %>% select(!!as.name(predictor),cat,fit,stage)
  
  test <- r %>% left_join(v,by='cat') %>% mutate(OR = exp(fit.y)/exp(fit.x)) %>% mutate(stage=stage.x)
  print (test)
  
  p1 <- plot_gam(gamFit,predictor,cat) + theme(legend.position="top")
  
  
  p1 <- p1 + geom_vline(xintercept = ref, linetype="dashed", color = "red", size = 1.) +
    geom_vline(xintercept = value, linetype="dashed", color = "red", size = 1.) +
    geom_segment(aes(x = !!as.name(paste0(predictor,'.x')), y = fit.x, xend = !!as.name(paste0(predictor,'.y')), yend = fit.y,color=cat), data = test,
                 size=2,show.legend = F,
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    geom_label(aes(x=!!as.name(paste0(predictor,'.y')),y=fit.y,label=round(OR,2)),data=test,
               show.legend=F,
               #position=position_jitter(height=0.5),
               hjust=-0.2
    ) +
    annotate("text",
             x = mean(c(ref,value)),
             y = 1,
             label = 'Odds Ratio(s)') 
  
  return(p1);
}

#plot_my_gam(pc.modelB.gamFit,'age',list(cat=levels(pc.modelB$cat)),20,40) #+ facet_grid(~stage) + ylim(-10,5)




plot_gam(pc.modelA.gamFit,'days_since_first_measurement',var_ref=0) +
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")

plot_gam(pc.modelA.gamFit,'days_since_last_vac',var_ref=15) +
  coord_flip() + scale_y_log10(labels = ~ signif(.x,digits=1)) + # limits=c(0.11,11) ) + 
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none")





temp <-  pc.modelB.gamFit.or %>% select(names) %>% full_join(bd.modelB.gamFit.or) %>% mutate(OR=paste0(round(OR,2)," (",round(LCL,2)," - ", round(UCL,2),")")) %>% 
         select(original,names,OR) %>% left_join(bd.modelB.gamFit %>% 
                                          mutate(OR=paste0(round(OR,2)," (",round(LCL,2)," - ", round(UCL,2),")")) %>% 
                                          rename(uOR=OR) %>% select(original,names,uOR)) %>%
         map_df(rev) %>%
         group_by(original) %>% mutate(original = replace(original,-1,"")) %>% ungroup 
         #rename(origin)

View(temp)

plot_strat(temp) + theme(legend.position = 'top', legend.justification = 'left') +
  coord_flip(ylim=c(0.02,50e3)) +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  geom_text(aes(x=names,y=150,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            data=temp%>%filter(strat=='Primary Care'),hjust=0, color=phs_colours("phs-blue")) +
  geom_text(aes(x=names,y=50e2,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            data=temp%>%filter(strat!='Primary Care'),hjust=0, color=phs_colours("phs-magenta"))





p1 <- plot_my_gam(gamFit,'days_since_first_measurement',NULL,0,200) + ylim(-5,-1) + guides(fill=FALSE, color=FALSE)
p2 <- plot_my_gam(gamFit,'days_since_last_vac',NULL,15,150) + guides(fill=FALSE, color=FALSE)
p2

grid.arrange(p1,p2,nrow=1)


oddsratio::insert_or(p1,test %>% mutate(predictor='age',value1=age.x,value2=age.y,oddsratio=OR) %>% head(1),
                     or_yloc = 5,
                     values_xloc = 0.04, line_size = 0.5,
                     line_type = "dotdash", text_size = 6,
                     values_yloc = 0.01, arrow_col = "red")


plot_gam(gamFit,'age') + theme(legend.position="top")


p1 <- plot_gam(gamFit,'age',list(cat=c('Pf','AZ'))) + theme(legend.position="top")
p1
p2 <- plot_gam(gamFit,'age',list(cat=c('Pf-Pf','AZ-AZ'))) + theme(legend.position="top")
modelA %>% group_by(cat) %>% summarise(n=n())
p3 <- plot_gam(gamFit,'age',list(cat=c('Pf-Pf-Pf','mixed 3 doses'))) + theme(legend.position="top")


grid.arrange(p1,p2,p3,ncol=3)

plot_gam(gamFit,'days_since_first_measurement') 
plot_gam(gamFit,'days_since_last_vac') 



gamFit <- perform_gam(modelB,as.formula(modelB.gam.formula))
or <- get_or_from_gam(gamFit)
or <- or %>% mutate_at(vars(-names),~round(.,2)) %>%
  mutate(pv = paste0('p ',ifelse(pv<0.05,'<0.05',paste0('= ',sprintf('%.2f',pv)))))
plot_simple(or,labels=F,xmax=20000) + 
  geom_text(aes(x=names,y=150,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            hjust=0)


#or_bd <- or
#or_pc <- or
#or_pc_tight <- or %>% filter(se<90)
or_pc_loose <- or %>% filter(!grepl("(Other)",names))%>% 
                    filter(!grepl("Cerebralpalsy",names))

#temp <- or_pc %>% mutate(strat='Primary Care') %>% rbind(or_bd %>% mutate(strat='Blood Donors'))
temp <- or_pc %>% mutate(strat='Impaired Response') %>% rbind(or_pc_tight %>% mutate(strat='Fail-to-Mount')) %>%
  rbind(or_pc_loose %>% mutate(strat='IgG Lowest 10%'))
temp

plot_strat(temp,'phs-rust') + theme(legend.position = 'top', legend.justification = 'left') +
  coord_flip(ylim=c(0.07,5e5)) +
  geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
  scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
  geom_text(aes(x=names,y=100,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            data=temp%>%filter(strat=='IgG Lowest 10%'),hjust=0, color=phs_colours("phs-green")) +
  geom_text(aes(x=names,y=20e2,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            data=temp%>%filter(strat=='Impaired Response'),hjust=0, color=phs_colours("phs-blue")) +
  geom_text(aes(x=names,y=50e3,label=paste0(sprintf('%.2f',OR),' (',sprintf('%.2f',LCL),' - ',sprintf('%.2f',UCL),'); ',pv)),
            data=temp%>%filter(grepl('Fail',strat)),hjust=0, color=phs_colours("phs-rust"))




gamFit <- perform_gam(modelB.tight,as.formula(modelB.tight.gam.formula))
or <- get_or_from_gam(gamFit)
plot_simple(or %>% filter(UCL<100))

  


#modelC <- modelA %>% filter(stage==3 | stage=='4+') %>% droplevels()
#gamFit.C <- perform_gam(modelC,as.formula(modelC.gam.formula))


modelC.B.gam.formula



gamFit.C  <- perform_gam(modelC.A %>% filter(!(Q_BMI=='0-20')),as.formula(modelC.A.gam.formula))
gamFit.C  <- perform_gam(modelC.B %>% filter(!(Q_BMI=='0-20')),as.formula(modelC.B.gam.formula))

or <- get_or_from_gam(gamFit.C)
or
#600x700
plot_simple(or)

plot_gam(gamFit.C)

modelD %>% group_by(Q_DIAG_BLOOD_CANCER) %>% summarise(n=n())

modelD.gam.formula <- "outcome ~ s(age) + Sex + Q_BMI  + Q_DIAG_CHD + Q_DIAG_DIABETES_1"
gamFit.D  <- perform_gam(modelD,as.formula(modelD.gam.formula))

or <- get_or_from_gam(gamFit.D)
or
#600x700
plot_simple(or)
  
  
  pred %>% filter(cat=='1 Dose AZ') %>% filter(age==50 | age==75 | age==100) %>% select(age,fit,se.fit)
  
  categories <- unique(pred$cat)
  categories
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

data


data %>% group_by(cat) %>% summarise(n=n())

#filter(cat=='AZ-AZ')
ggplot(data,aes(x=age, y=or, color=cat)) +
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





plot_gam(gamFit, pred = "days_since_last_vac")



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



plot_simple(or)

write.csv(or%>%select(-names),"/home/calumm09/temp.csv")


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


          

qnames <- get_qnames(modelC,n=5)
qnames

formula <- paste0("outcome ~ s(age,by=cat) + Sex  +  Q_BMI + simd2020v2_sc_quintile + 
                      ",paste0(qnames,collapse = " + ")," +
                      shielding + immuno_supp +
                      prior_infection + s(days_since_last_vac) + s(days_since_first_measurement) +
                      cat ")


formula <- paste0("outcome ~ s(age,by=cat) + Sex + Q_BMI + simd2020v2_sc_quintile +
                      ",paste0(qnames,collapse = " + ")," +
                      cat ")


gamFit 
gamFit <- perform_gam(modelC,as.formula(formula))

modelC %>% group_by(Q_BMI,outcome) %>% summarise(n=n())


modelC %>% group_by(outcome) %>% summarise(n=n())

hist(modelC$IgG)


modelC %>% group_by(simd2020v2_sc_quintile,outcome) %>% summarise(n=n()) %>% ungroup %>%
  group_by(simd2020v2_sc_quintile) %>%
  mutate(per =  100 *n/sum(n)) %>% ungroup 




or_glm(data=modelC, model=gamFit)
c <- coef(gamFit)[['Q_DIAG_CKD_LEVEL']]
exp(c)

term_list <- list()
for (term in labels(gamFit$terms)){
  new_term <- gamFit[["var.summary"]][[term]][[1]]
  term_list <- append(term_list, list(new_term))
}

names(term_list) <-  labels(gamFit$terms)
term_list[['age']] <- 50
term_list[['days_since_first_measurement']] <- 300
term_list[['days_since_last_measurement']] <- 21
term_list

new_data <- expand.grid(term_list)
new_data

pred <- predict.gam(gamFit,new_data,se.fit = TRUE)
pred <- cbind(new_data, pred)
pred


exp(pred %>% select(fit,se.fit))

tor <- or_gam(data = modelC, model = gamFit, pred = 'Q_DIAG_CKD_LEVEL', values=c(0,1))
tor


or <- get_pterm_or_from_gam(gamFit,modelC)
or
pgam <- plot_simple(or,xmax=20,xmin=0.15)
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
