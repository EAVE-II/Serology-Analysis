library(ggplot2)
library(dplyr)

f_normal <- ~  scale*(exp(-0.5*(time-mu)^2/sd^2))
normal <- deriv(f_normal,
                namevec=c("scale","mu","sd"),
                function.arg=c("time","scale","mu","sd"))


f_normal2 <- ~  start + scale*(exp(-0.5*(time-mu)^2/sd^2))
normal2 <- deriv(f_normal,
                namevec=c("scale","mu","sd","start"),
                function.arg=c("time","scale","mu","sd","start"))



get_df <- function(){
  
  df <- NULL
  
  m_scale <- 150
  m_days <- 30.
  m_sd <- 3

  for (sex in c('M','F')){
    
    size <- 1000
    s_days <- 0
    if(sex == 'F'){
      size <- 1500
      s_days <- 5.
    }
    
    for (product in c('AZ','Pf')){
    
      p_sd <- 1
      p_scale <- 1.
      p_days <- 0
      if(product == 'AZ'){
        p_sd <- 1.0
        p_scale <- 0.9
        p_days <- 0.
      }
    
      x <- sample(seq(1,200,0.1),size=size,replace = T)
      
      dose <- sample(0:1,size=length(x),replace=T)
      
      df1 <- data.frame(x=x,Sex=sex,product=product,dose=dose) %>% rowwise() %>%
             mutate(
                    y= rnorm(1,m_scale*p_scale*0.5*dose,m_scale*p_scale*0.1*dose) + 
                      rnorm(1,1,0.1)*normal(log(x),
                                            m_scale*p_scale,
                                            log(m_days + s_days + p_days), 
                                            log(m_sd)#*p_sd
                                            )
             ) %>% ungroup
      
    
      df <- df %>% rbind(df1)
    }
  }
  return (df %>% mutate(id=row_number()));
}

df <- get_df()

df %>% ggplot(aes(x=x,y=y,color=Sex)) + geom_point() + facet_grid(product ~ dose)


startvec <- c(scale=150,mu=log(30),sd=log(3))
m0 <- nls(y ~ normal(log(x), scale, mu, sd),
          df,
          start=startvec
          #lower=c(start=0,scale=0,mu=-10,sd=-10),
)
startvec <- coef(m0)
startvec['start'] = 1
startvec

fit <- nlmer(y ~  start + normal(log(x), scale, mu, sd) ~
               + (scale|Sex) + (mu|Sex) + (sd|Sex) 
               + (start|dose), 
               #+ (scale|product) + (mu|product) + (sd|product),
             df,
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             verbose=1,
             start = startvec)
summary(fit)


results <- get_nlmer_results(fit,modify=T)

scaled_to_normal <- function(term,.){
  case_when(
    term=='mu' ~ exp(.),
    term=='sd' ~ exp(.),
    #term=='scale' ~ .,
    TRUE ~ .
  )
}
res <- results$results %>% mutate_at(c('estimate','conf.low','conf.high'),~ scaled_to_normal(term,.))

plot_nlmer_results(res) 


pars <- (results$results %>% filter(group=='Sex' & level=='M') %>% as.list)$estimate
scale <- pars[[1]]
mu <- pars[[2]]
sd <- pars[[3]]

data.frame(x=seq(0,200,0.1)) %>%
        mutate(y=normal(log(x),scale,mu,sd)) %>%
  ggplot(aes(x=x,y=y)) +
  geom_line(linetype='dashed')


dasd




sdas























































startvec <- c(scale=100,mu=3.,sd=1.)

m0 <- nls(y ~ normal(log(x),scale,mu,sd),
          df,
          start=startvec
          #lower=c(start=0,scale=0,mu=-10,sd=-10),
)
startvec <- coef(m0)
startvec


fit <- nlmer(y ~ normal(log(x), scale, mu, sd) ~ 
               + (Sex|scale) + (Sex|mu) ,
             df,
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)),
             verbose=1,
             start = startvec)
summary(fit)


fixed <- fixef(fit)
fixed
random <- ranef(fit)
random

results <- broom.mixed::tidy(fit, effects = "ran_vals", conf.int = TRUE)
results <- results %>% rowwise %>% mutate_at(c('estimate','conf.low','conf.high'), 
                                             ~ . + fixed[group][[1]]) %>% ungroup 


male <- results %>% filter(term=='SexF' & level=='1')
male

scale <- male$estimate[[1]]
mu <- male$estimate[[2]]
sd <- fixed['sd'][[1]]

scale.low <- male$conf.low[[1]]
mu.low <- male$conf.low[[2]]
sd.low <- fixed['sd'][[1]]

scale.high <- male$conf.high[[1]]
mu.high <- male$conf.high[[2]]
sd.high <- fixed['sd'][[1]]

x <- seq(1,100,0.1)



df2 <- data.frame(x=x) %>%
  mutate(
    y=normal(log(x),scale,mu,sd),
    ymin=normal(log(x),scale.low,mu.low,sd.low),
    ymax=normal(log(x),scale.high,mu.high,sd.high)
  )


df %>% ggplot(aes(x=x,y=y)) + 
  geom_point(aes(color=Sex)) +
  geom_ribbon(aes(x=x,y=y,ymax=ymax,ymin=ymin),data=df2) 
#geom_line(aes(x=x,y=y),data=df3)



