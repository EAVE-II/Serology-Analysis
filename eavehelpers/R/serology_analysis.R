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
library(data.table)
library(forcats)
library(broom)
library(survival)

library(ggtext)
library(kableExtra)
library(rlang)
library(visreg)

#' Helper function that will extract the QCOVID names from a dataframe
#' after applying filtering on a condition,
#' and requiring there to be >n events with condition
#'
#' It is used to find a QCOVID variables that have >n events after applying a filtering
#'
#' @param df input dataframe object
#' @param var the name of the variable in the input dataframe to be checked
#' @param condition the value of the var to be filtered on (default 'Yes')
#' @param n the number of events required to be present for the filtered condition  (default 5)
#' @return a list of QCOVID names
#' @export
get_qcovid_names <- function(df,var=insufficient_response,condition='Yes',n=5){
  #extract the variable
  var <- rlang::quo_name(rlang::enquo(var))
  qnames <- df %>% filter(!!sym(var)==condition) %>% select(contains("Q_DIA")) %>%
    summarise(across(where(is.factor), ~ sum(.x == 'Yes', na.rm = TRUE))) %>%
    gather() %>%
    filter(value>=n)
  #just return the names
  qnames <- qnames$key
  return (qnames)
}



#' Work out which variables are good enough to be used for analysis
#' based upon how good the p-value is for them when used in a univariate GLM fit with the outcome
#' @param model the input dataframe model
#' @param outcome_var the name of the outcome variable (where outcome_var ~ x + y + z)
#' @param thres the threshold value for the p-value (default 0.1)
#' @return a list of variable names
#' @export
get_vars_to_use <- function(model,outcome_var,thres=0.1){

  #get all variables expect the outcome variable
  variables <- model %>% select(-!!sym(outcome_var)) %>% names()

  dfp <- NULL
  #loop over all variables
  for (name in variables){
    #perform GLM with the variable fitted to the outcome variable
    formula <- paste0(outcome_var," ~ ",name)
    glmFit <- tryCatch({
        glm(formula,family=binomial,data=model)
      },
      error=function(e){
        message(e)
        return (NULL)
      }
    )
    if(is.null(glmFit)){
      next
    }

    glmFit.summary <- summary(glmFit)
    #get the p-value for this GLM univariate fit
    p = coef(glmFit.summary)[,"Pr(>|z|)"][[2]]
    #save the pvalue for this variable (name)
    temp <- data.frame(name,p)
    dfp = rbind(dfp, temp)
  }
  #return which variables have a p-value less than the specified threshold
  vars_to_use <- as.vector((dfp %>% filter(p < thres))$name)
  return (vars_to_use)
}

#' Hacky bit of code to classify product-stage (vaccine dose) combinations
#' could be done a much better way... but it works...
#' should be used on already prepared data
#' @param stage integer saying which stage of the vaccine was it when the serology measurement was taken
#' @param d1,d2,d3,d4 vaccine product used at each stage of the vaccine
#' @export
calculate_product_category <- function (stage,d1,d2,d3,d4){

  #for up to a possible 4 vaccines, given the 'stage' (defined by serology measurement date)
  # define variables such that:
  # 0 - vaccines not administered (at time of this stage)
  # 1 - AZ administered (at this stage)
  # 2 - Anything else (Pfizer, Moderna) administered (at this stage)

  d1 = ifelse(d1 == 'Covid-19 Vaccine AstraZeneca',1,2)
  d2 = ifelse(is.na(d2),0,ifelse(d2 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d3 = ifelse(is.na(d3),0,ifelse(d3 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d4 = ifelse(is.na(d4),0,ifelse(d4 == 'Covid-19 Vaccine AstraZeneca',1,2))
  d2 <- ifelse(stage<2,0,d2) # d2 has to be 0 as stage is <2
  d3 <- ifelse(stage<3,0,d3)
  d4 <- ifelse(stage<4,0,d4)

  # e.g. 2211 - two doses AZ followed by 2 doses of non-AZ
  cat <- as.integer(paste0(d4,d3,d2,d1))

  #convert to more interpretable representation
  cat <- case_when(
    cat==1 ~ 'One dose AZ',
    cat==2 ~ 'One dose Pfizer/Moderna',
    cat==11 ~ 'Two doses of AZ',
    cat==22 ~ 'Two doses of Pfizer/Moderna',
    cat==222 | cat==2222 ~ 'Mixed 3+ doses (no AZ)',
    TRUE ~ 'Mixed 3+ doses (including AZ)'
  )

  #convert to a factor
  cat <- factor(cat,
                levels = c('Two doses of Pfizer/Moderna',
                           'One dose AZ',
                           'One dose Pfizer/Moderna',
                           'Two doses of AZ',
                           'Mixed 3+ doses (including AZ)',
                           'Mixed 3+ doses (no AZ)'))
  return (cat)
}

#' Helper function to code the key variables used in the serology analysis
#' e.g. making factor variables
#' @param model
#' @return modified model
#' @export
code_vars <- function(model){

  #start coding some variable with references (first factor level)
  model <- model %>%
            mutate(Sex = as.factor(Sex), #make into a factor
                   age=ageYear, #keep a record of the original age before categorisation
                   #n_risk_gps=as.factor(as.integer(n_risk_gps)-1), #coversion from <ord> to <factor>
                   cat = calculate_product_category(stage,d1_product,d2_product,d3_product,d4_product)#,
            )
  #setup SIMD, factor version of age and factor version of BMI
  # BMI levels are defined by standard defintions of under/over weight and obese
  model <- model %>%
            mutate_at(vars(one_of('simd2020v2_sc_quintile')), ~relevel(as.factor(.x),'3')) %>% #SIMD=3 as reference
            mutate_at(vars(one_of('ageYear')), ~factor(cut(.x, breaks = c(0,20,40,60,150), right = T,
                                                       labels = c("0-19","20-39","40-59","60+")),
                                                       levels=c("40-59","0-19","20-39","60+"))) %>%
            mutate_at(vars(one_of('Q_BMI')), ~relevel(cut(ifelse((is.na(.x) | .x<5 | .x>80 ),-1,.x),
                                                      breaks = c(-2,5,18.5,25,30,80), right = F,
                                                      labels = c("Unknown","<18.5","18.5-25","25-30","30+")),
                                                      ref='18.5-25')
            )

  #create a better variable for the definition of immuno suppressed
  #we have a duplication from QCOVID and the vaccine record
  #we decided to combine them into a factor level
  #Note: we dropped the Q_COVID - Q_DIAG_IMMU variable as it double counts,
  #      but is also very rare (so was never used in any analysis)
  model <- model %>%
          mutate(immuno = case_when(
                                    (severely_immuno_supp==1)  ~ 'Severely', #| (Q_DIAG_IMMU>0)
                                    immuno_supp==1 ~ 'Yes',
                                    TRUE ~ 'No')) %>%
          mutate(immuno = factor(immuno,levels=c('No','Yes','Severely'))) %>% ##ifelse(immuno==1,0,immuno))) %>%
          select(-immuno_supp,-severely_immuno_supp)#,-Q_DIAG_IMMU)

  #make these binary variables into a Yes/No
  model <- model %>%
          mutate_at(c("insufficient_response","prior_infection","ch_resident","shielding"),
                    ~as.factor(ifelse(.==1,'Yes','No'))) %>%
          mutate_at(vars(contains('Q_DIAG')),
                    ~as.factor(ifelse(.==1,'Yes','No'))) %>%
          mutate_at(vars(contains('additional')),
                    ~as.factor(.))

  return (model)
}

#' Simple function to get the dataframe for modelA given a cohort dataframe
#' @param df input data for the cohort
#' @return cleaned and mutated version of the data for analysis
#' @export
get_modelA_df <- function(df){
  require(tidyr)
  return (df %>% code_vars())
}

#' Simple function to get the dataframe for modelB given a cohort dataframe
#' @param df input data for the cohort
#' @return cleaned and mutated version of the data for analysis
#' @export
get_modelB_df <- function(df){
  require(tidyr)
  return (df %>% filter(n_risk_gps!=0) %>% code_vars())
}




#' Perform an unadjusted gam fit for each variable in a formula given a dataframe
#' @param model input dataframe to be fitted to the GAM model
#' @param formula a full formula for a (possible) multivariate game
#' @export
calculate_unadjusted_gam <- function(model,formula){
  #split the formula (e.g. y ~ x + z) to extract all the individual terms (x,z)
  terms <- labels(terms(formula))
  #extract the 'outcome' variable, i.e. the y in y ~ x + z + ...
  outcome <- all.vars(formula)[1]
  fits <- NULL
  #loop over all terms in the formula
  for (var in terms){
    #create a new univariate formula for the variable term
    f <- as.formula(paste0(outcome," ~ ",var))
    #perform the gam fit with this formula, on the input dataframe
    fit <- perform_gam(model,f)
    #book the fit results
    fits[[var]] <- fit
  }
  #return all the fit results
  return (fits);
}

#' Perform a GAM fit of a model given a formula
#' in the serology GAM analys(es) we always use binomial logistic regression
#' .. so this is a simple convienence function
#'
#' @param model input dataframe to be used
#' @param formula formula string (e.g. y ~ x + z)
#' @export
perform_gam <- function(model,formula){
  require(mgcv)
  fit <- gam(formula,family=binomial,data=model)
  return(fit)
}

#' ..... REDUNDANT FUNCTION......
#' Extract the parametric terms from a GAM fit
#' .. these are all variables that are not splines or interactions..
#' @param fit the GAM fit object
#' @param model the input dataframe used for the model
#' @export
get_pterm_or_from_gam <- function (fit,model) {
  #extract all parametric term labels
  pterms <- labels(fit$pterms)
  or <- NULL
  #loop over each term
  for(term in pterms){
    #extract all possible levels from the varible
    #this works best when all these terms are factors (which is what we do for this analysis)
    levs <- levels(model[[term]])
    if (is.null(levs)){
      levs <- unique(sort(model[[term]]))
    }
    else{
      levs <- levels(droplevels(model[[term]]))
    }
    #get the first level as the reference
    ref <- levs[1]
    #get all other levels to calculate ORs with respect to the reference level
    values <- levs[2:length(levs)]
    #loop over the values
    for (value in values){
      #use the or_gam function to get the ORs (and CIs) of this parametic term
      tor <- or_gam(data = model, model = fit, pred = term, values=c(ref,value))
      #book this variable
      or <- or %>% rbind (tor)
    }
  }
  #do some cleaning up of the mini dataframe used to save the ORs and return
  or <- or %>% unite('names',sep='',c(predictor,value2)) %>% select(-value1)
  names(or)[2] <- 'OR'
  names(or)[3] <- 'LCL'
  names(or)[4] <- 'UCL'
  return (or)
}


#' Get a map between the coefficient name and
#'  - the original variable it came from
#'  - the level it corresponds to
#' This function is for use with factor levels in GAM Fits
#' @param gamFit GAM fit object from mgcv
#' @export
get_coeff_name_map <- function(gamFit){
  vars <- gamFit$var.summary
  retval <- list()
  for (varname in labels(vars)) {
    levs = levels(vars[[varname]])
    if(is.null(levs)){
      levs <- unique(vars[[varname]])
    }

    for (lev in levs){
      mlev = paste0(varname,lev)
      retval[[mlev]] = list(var=varname,level=lev)
    }
  }
  return (retval);
}


#' Extract ORs from a GAM fit for the parametric variables
#' @param mod the GAM fit object from mgcv
#' @export
get_or_from_gam <- function(mod, level = 0.95) {
  # a method for extracting confidence intervals and returning a tidy data frame
  # - found online via stackexchange
  require(mgcv)
  require(dplyr)

  #get the names of all variables used in the GAM
  name_map <- get_coeff_name_map(mod)
  #get unadjusted ORs that should have been already calculated and
  # added as an attribute to this model
  unadjusted <- attributes(mod)$unadjusted

  #from the name map, lookup a variable name
  get_var_name <- function(name){
    var <- name_map[name] %>% map(c('var')) %>% as.character
    return (var)
  }
  #get the level name from the name map
  get_level_name <- function(name){
    lev <- name_map[name] %>% map(c('level')) %>% as.character
    return (lev)
  }

  get_odds <- function(mod){
    mod.s <- summary(mod)
    print (mod.s)
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
      mutate(label=ifelse(is.na(label),var,label)) %>%
      #mutate(var = ifelse(is.null(var),names,var)) %>%
      as_tibble()
    return (tab)
  }

  #get the odds ratios for all parametric terms in this GAM Fit
  tab <- get_odds(mod)
  #hacky, but just get back the names of all variables we calculate ORs for
  vars <- unique(tab$var)

  utab <- NULL
  #loop over these variables
  for (var in vars){
    if (is.null(unadjusted[[var]])) next
    #get the odds for the univariate gam fits
    temp <- get_odds(unadjusted[[var]]) %>% select(names,OR,LCL,UCL) %>%
            rename(uOR=OR,uLCL=LCL,uUCL=UCL)
    utab <- rbind(temp,utab)
  }

  #hacky code for building up tables for convienence
  if(!is.null(utab)){
    tab <- tab %>% left_join(utab)
  }

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
             label=eavehelpers::get_label(var)) %>%
      mutate(label=ifelse(is.na(label),var,label))
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
  else{
    term_list[[var]] <- seq(0,300,1)
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
  #pred <- pred %>% mutate(cat=eavehelpers::get_label(paste0('cat',cat)))
  pred <- pred %>% mutate(cat='1')
  pred <- pred %>% mutate(stage=ifelse(grepl('1',cat),'Dose 1',ifelse(grepl('2',cat),'Dose 2','Dose 3')))


  return (pred);
}

#' @export
plot_gam <- function(gamFit,var, var_ref=0, cats = NULL) {


  palette <- 'main'
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
    labs(y='Odds Ratio',x=eavehelpers::get_label(var),color='',fill='') +
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
get_ratios_table <- function (tab){
  tab %>%
    mutate(uORs=ifelse(pos==1,"-",paste0(sprintf("%.2f",round(uOR,2))," (",sprintf("%.2f",round(uLCL,2))," - ", sprintf("%.2f",round(uUCL,2)), ")"))) %>%
    mutate(ORs=ifelse(pos==1,"-",paste0(sprintf("%.2f",round(Nominal,2))," (",sprintf("%.2f",round(LCL,2))," - ", sprintf("%.2f",round(UCL,2)), ")"))) %>%
    mutate(label = replace(label, duplicated(label), '')) %>%
    select(var,label,level,uORs,ORs) %>% return
}

#' @export
display_ratios_table <- function (tab){
  require(knitr)
  get_ratios_table(tab) %>%
    kable() %>%
    kable_classic(full_width = F, html_font = "Cambria")
}


#' @export
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
  toShow <- toShow %>% as_tibble() %>%
    mutate_at(c('estimate','conf.low','conf.high'),exp) %>%
    mutate(name=as.character(var)) %>%
    mutate(name=eavehelpers::get_label(name)) %>%
    mutate(name=ifelse(is.na(name),as.character(var),name))

  levels <- unique(toShow$name)
  #toShow <- toShow %>% mutate(name=factor(name,levels=levels))

  toShow <- toShow %>% mutate(level = ifelse(pos==1,paste0('<b>',as.character(level),' (ref) </b>'),as.character(level)))
  return (toShow);
}


