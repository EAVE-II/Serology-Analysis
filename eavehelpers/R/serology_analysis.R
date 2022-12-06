library(tibble)
library(dplyr)
library(tidyr)
library(mgcv)
library(tidymv)
library(gratia)
library(purrr)
library(markdown)
library(ggtext)
library(knitr)
library(kableExtra)


#' @export
get_qcovid_names <- function(df,var=insufficient_response,condition='Yes',n=5){
  var <- rlang::quo_name(rlang::enquo(var))
  qnames <- df %>% filter(!!sym(var)==condition) %>% select(contains("Q_DIA")) %>% 
    summarise(across(where(is.factor), ~ sum(.x == 'Yes', na.rm = TRUE))) %>%
    gather() %>% 
    filter(value>=n)
  qnames <- qnames$key
  return (qnames)
}



#' Work out which variables are good enough to be used for analysis 
#' based upon how good the p-value is for them when used in a univariate GLM fit with the outcome
#' @export
get_vars_to_use <- function(model,outcome_var,thres=0.1){
  
  #get all variables expect the outcome variable
  variables <- model %>% select(-!!sym(outcome_var)) %>% names() #select(-any_of(c("outcome","ageYear","Sex","n_risk_gps")))  %>% names()

  dfp <- NULL
  #loop over all variables
  for (name in variables){
    #perform a quick GLM with the variable
    formula <- paste0(outcome_var," ~ ",name)
    glmFit <- tryCatch({
      glm(formula,family=binomial,data=model)
    },
    error=function(e){
      message(e)
      return (NULL)
    })
    if(is.null(glmFit)){
      next
    }
    
    glmFit.summary <- summary(glmFit)
    #get a p-value
    p = coef(glmFit.summary)[,"Pr(>|z|)"][[2]]
    temp <- data.frame(name,p)
    #names(temp) <- c('name','p')
    dfp = rbind(dfp, temp)
  }
  #return which variables have a p-value less than the specified threshold
  vars_to_use <- as.vector((dfp %>% filter(p < thres))$name)
  return (vars_to_use)
}

#' Hacky bit of code to classify product-stage (vaccine dose) combinations
#' @export
calculate_product_category <- function (stage,d1,d2,d3,d4){
  
  #for up to a possible 4 vaccines, given the 'stage' (defined by serology measurement date)
  # define variables such that:
  # 0 - vaccines not administered (at time of this stage)
  # 1 - AZ administered (at this stage)
  # 2 - Anything else (Pfizer, little Moderna) administered (at this stage)
  
  d1 = ifelse(d1 == 'Covid-19 Vaccine AstraZeneca',1,2)
  d2 = ifelse(is.na(d2),0,ifelse(d2 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d3 = ifelse(is.na(d3),0,ifelse(d3 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d4 = ifelse(is.na(d4),0,ifelse(d4 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d2 <- ifelse(stage<2,0,d2) # d2 has to be 0 as stage is <2 
  d3 <- ifelse(stage<3,0,d3)
  d4 <- ifelse(stage<4,0,d4)

  # e.g. 2211 - two doses AZ followed by 2 doses of non-AZ   
  cat <- as.integer(paste0(d4,d3,d2,d1))
  
  cat <- case_when(
    cat==1 ~ 'One dose AZ',
    cat==2 ~ 'One dose Pfizer/Moderna',
    cat==11 ~ 'Two doses of AZ',
    cat==22 ~ 'Two doses of Pfizer/Moderna',
    cat==222 | cat==2222 ~ 'Mixed 3+ doses (no AZ)',
    TRUE ~ 'Mixed 3+ doses (including AZ)'
  )

  cat <- factor(cat,
                levels = c('Two doses of Pfizer/Moderna',
                           'One dose AZ',
                           'One dose Pfizer/Moderna',
                           'Two doses of AZ',
                           'Mixed 3+ doses (including AZ)',
                           'Mixed 3+ doses (no AZ)'))
  return (cat)
}

#' @export
code_vars <- function(model){
  
  #start coding some variable with references (first factor level)
  model <- model %>% 
            mutate(Sex = as.factor(Sex), #make into a factor
                   age=ageYear, #keep a record of the original age before categorisation
                   n_risk_gps=as.factor(as.integer(n_risk_gps)-1), #coversion from <ord> to <factor>
                   cat = calculate_product_category(stage,d1_product,d2_product,d3_product,d4_product)#,
            )
  #setup SIMD, factor version of age and factor version of BMI
  # BMI levels are defined by standard defintions of under/over weight and obese
  model <- model %>% 
            mutate_at(vars(one_of('simd2020v2_sc_quintile')), ~relevel(as.factor(.x),'3')) %>% #SIMD=3 as reference
            mutate_at(vars(one_of('ageYear')), ~factor(cut(.x, breaks = c(0,20,40,60,150), right = T, 
                                                   labels = c("0-19","20-39","40-59","60+")),
                                               levels=c("40-59","0-19","20-39","60+"))) %>%
            mutate_at(vars(one_of('Q_BMI')), ~relevel(cut(ifelse((is.na(.x) | .x<5 | .x>80 ),-1,.x), breaks = c(-2,0,20,25,30,80), right = T, 
                                                  labels = c("Unknown","0-20","20-25","25-30","30+")),ref='20-25')
    )
  
  #create a better variable for the definition of immuno suppressed 
  #we have a duplication from QCOVID and the vaccine record
  #we decided to combine them into a factor level
  #Note: we dropped the Q_COVID - Q_DIAG_IMMU variable as it double counts, 
  #      but is also very rare (so was never used in any analysis)
  model <- model %>%
          #mutate(immuno = strtoi(paste0(immuno_supp,severely_immuno_supp,0))) %>%
          mutate(immuno = case_when(
                                    (severely_immuno_supp==1) | (Q_DIAG_IMMU>0) ~ 'Severely',
                                    immuno_supp==1 ~ 'Yes',
                                    TRUE ~ 'No')) %>%
          mutate(immuno = factor(immuno,levels=c('No','Yes','Severely'))) %>% ##ifelse(immuno==1,0,immuno))) %>%
          select(-immuno_supp,-severely_immuno_supp,-Q_DIAG_IMMU)
  
  #make these binary variables into a Yes/No
  model <- model %>%
          mutate_at(c("insufficient_response","prior_infection","ch_resident","shielding"),
                    ~as.factor(ifelse(.==1,'Yes','No'))) %>% 
          mutate_at(vars(contains('Q_DIAG')),
                    ~as.factor(ifelse(.==1,'Yes','No')))
          
  return (model)
}

#' @export
get_modelA_df <- function(df){
  require(tidyr)
  return (df %>% code_vars())
}

#' @export
get_modelB_df <- function(df){
  require(tidyr)
  return (df %>% filter(n_risk_gps>0) %>% code_vars())
}


#' @export
calculate_unadjusted_gam <- function(model,formula){
  terms <- labels(terms(formula))
  outcome <- all.vars(formula)[1]
  fits <- NULL
  for (var in terms){
    f <- as.formula(paste0(outcome," ~ ",var))
    fit <- perform_gam(model,f)
    fits[[var]] <- fit
  }
  return (fits);
}



#' @export
perform_gam <- function(model,formula){
  fit <- gam(formula,family=binomial,data=model)
  return(fit)
}

#' @export
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
  #or <- or %>% mutate( names=ifelse(is.na(as.character(lookup[names])),
  #                                  ifelse(is.na(as.character(lookup[gsub('.{1}$', '', names)])),names,lookup[gsub('.{1}$', '', names)]),
  #                                  lookup[names]))
  return (or)
}


#' Get a map between the coefficient name and 
#'  - the original variable it came from
#'  - the level it corresponds to
#' This function is for use with factor levels in GAM Fits
#' @export 
get_coeff_name_map <- function(gamFit){
  vars <- gamFit$var.summary
  retval <- list()
  for (varname in labels(vars)) {
    levs = levels(vars[[varname]])
    for (lev in levs){
      mlev = paste0(varname,lev)
      retval[[mlev]] = list(var=varname,level=lev)
    }
  }
  return (retval);
}


#' @export
get_or_from_gam <- function(mod, level = 0.95) {
  # a method for extracting confidence intervals and returning a tidy data frame
  # - found online via stackexchange 
  require(mgcv)
  require(dplyr)
  
  name_map <- get_coeff_name_map(mod)

  unadjusted <- attributes(mod)$unadjusted

  
  get_var_name <- function(name){
    var <- name_map[name] %>% map(c('var')) %>% as.character
    #var <- ifelse(is.null(var),'test',var)
    return (var)
  }
  get_level_name <- function(name){
    lev <- name_map[name] %>% map(c('level')) %>% as.character
    return (lev)
  }
  
  get_odds <- function(mod){
    mod.s <- summary(mod)
    #get the coefficients for non-splined terms
    #. p.coeff. is an array of estimates of the strictly parametric model coefficients
    # - also get the pvalue of these parametirc coefficients
    E <- data.frame(OR = mod.s$p.coeff, pv = mod.s$p.pv) %>%
      mutate(names = row.names(.)) %>%
      select(names, OR, pv)
    
    #get the se separately 
    SE <- data.frame(se= mod.s$se) %>%
      mutate(names = row.names(.)) %>%
      select(names, se)
    #get the deviance
    nu <- mod.s$residual.df
    
    #join them up so there's an estimate and error for each parametric coefficient
    tab <- inner_join(E, SE) %>%
      mutate(
        LCL = OR + #calculate the lower confidence interval given the level (default is 95%)
          se * qt(df = nu,
                  p = (1 - level) / 2),
        UCL = OR + #same for upper confidence
          se * qt(df = nu,
                  p = 1 - (1 - level) / 2)) %>% 
      slice(-1) %>% #remove the first - this is always the intercept coefficient which we dont want any OR from 
      mutate_at(vars(OR,LCL,UCL),exp) #%>%  #expontentiate the OR (and CIs) - to get the actual OR (exp of coeff) 
    
    tab <- tab %>%
      mutate(var = get_var_name(names),
             level = get_level_name(names),
             label = eavehelpers::get_label(var)) %>%
      #mutate(var = ifelse(is.null(var),names,var)) %>% 
      as_tibble() 
    return (tab)
  }
  
  tab <- get_odds(mod)
  vars <- unique(tab$var)
  
  utab <- NULL
  for (var in vars){
    temp <- get_odds(unadjusted[[var]]) %>% select(names,OR,LCL,UCL) %>%
            rename(uOR=OR,uLCL=LCL,uUCL=UCL)
    utab <- rbind(temp,utab)
  }
  
  tab <- tab %>% left_join(utab)
  
  tab <- tab %>% group_by(var) %>% mutate(pos=row_number()+1) %>% ungroup
  references <- do.call(rbind, lapply(unique(tab$var), function(var){
    tibble(var=var,
           level=levels(mod$var.summary[[var]])[1],
           pos=1,
           OR=1,
           UCL=1,
           LCL=1) %>% 
      mutate(names=paste0(var,level),
             level = paste0("<b>",level, " (ref)", "</b>"),
             label=eavehelpers::get_label(var))
  }))
  
  tab %>% bind_rows(references) %>% arrange(var,pos) %>% return
}


#' @export
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

#' @export
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
    labs(y='Odds Ratio',x=lookup[var],color='',fill='') + 
    scale_colour_discrete_phs(palette=palette) +
    geom_vline(xintercept=var_ref,linetype='dotted') +
    geom_hline(yintercept=1,linetype='dotted') +
    scale_fill_discrete_phs(palette=palette) +
    #scale_y_log10() +
    theme_classic()
  
  if(palette=='none'){
    p <- p + theme(legend.position="none")
  }
  return(p)
}

#' @export
add_names <- function(df,formula){
  labs <- labels(terms(as.formula(formula)))
  temp <- df %>% select(original_name) %>% as_tibble %>% mutate(original=list(labs)) %>%
    unnest(original) %>%
    filter(unlist(Map(function(x, y) grepl(x, y), original, original_name))) 
  return (df %>% left_join(temp))
}


#' @export
plot_ratios <- function(tab,fill=NULL,facet=NULL,
                        xmin=0.01,xmax=100,xtitle='OR (95% CI)',
                        colors=NULL){
  require(markdown)
  require(ggtext)
  
  
  if(is.null(colors)){
    colors = c(phs_colours("phs-blue"), phs_colours("phs-rust"), phs_colours("phs-green"))
  }
  
  p <- ggplot(tab, aes(x=reorder(level,-pos)))
  
  fill_var =  rlang::quo_name(rlang::enquo(fill))

  if (fill_var %in% colnames(tab)){
    p <- p + 
      geom_pointrange(aes(y=Nominal, ymin=LCL, ymax=UCL, fill=!!sym(fill_var), color=!!sym(fill_var)), width=.2,alpha=1,shape=21,
                      position=position_dodge(.9))
  }
  else {
    p <- p + 
      geom_pointrange(aes(y=Nominal, ymin=LCL, ymax=UCL, fill='adjusted', color='adjusted'), width=.2,alpha=1,shape=21,
                      position=position_dodge(.9)) 
  }
  

  p <- p + scale_fill_manual(values=colors)  + scale_color_manual(values=colors) 
  
  p <- p + 
    labs(title='',x='',y=xtitle) +
    geom_hline(yintercept=c(1), linetype="dashed") +
    geom_hline(yintercept=c(0.1,0.333,3,10), linetype="dotted") +
    #ylim(0., 10) +
    scale_y_log10(breaks=c(0.01,0.1,0.33,1,3,10)) +
    coord_flip(ylim=c(xmin,xmax)) +
    theme_classic() + guides(fill="none",color="none") +
    facet_grid(as.factor(label) ~ .,scale='free',space='free', switch='both') +
    theme(
      panel.background = element_rect(fill = NA, color = "black"),
      strip.background = element_blank(),
      strip.placement = 'outside',
      axis.text.y = element_markdown(angle=0),
      strip.text.y.left = element_text(angle = 0)) 
  
  if (!is.null(facet)){
    facet_formula =  as.formula(rlang::quo_name(rlang::enquo(facet)))
    p <- p + facet_grid(facet_formula ,scale='free',space='free',switch='y' ) +
             theme(
                  panel.background = element_rect(fill = NA, color = "black"),
                  strip.background = element_blank(),
                  strip.placement = 'outside',
                  axis.text.y = element_markdown(angle=0),
                  strip.text.y.left = element_text(angle = 0)
            )
  }
  
  return (p)
}

#' @export
display_ratios_table <- function (tab){
  require(knitr)
  tab %>% 
    mutate(Ratio=ifelse(pos==1,"-",paste0(round(Nominal,2)," (",round(LCL,2)," - ", round(UCL,2), " )"))) %>% 
    mutate(label = replace(label, duplicated(label), '')) %>% 
    select(label,level,Ratio) %>% 
    kable() %>%
    kable_classic(full_width = F, html_font = "Cambria")
}


