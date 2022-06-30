library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

install.packages("ggpubr")

colnames(df_ana_2)


get_kmeans <- function(df,nk=5){
  km <- kmeans(df,centers=nk,nstart=20)
  clusters <- km$cluster
  return (km$centers)
}

df <- df_ana_2 %>%  filter(outcome==1  & product_binary == 1 & stage != 1 & !is.na(Q_BMI) & Q_BMI>0 & Q_BMI<80 ) 
qnames <- df %>% get_qnames()
length(qnames)
qnames
df %>% group_by(Q_DIAG_CKD_LEVEL) %>% summarise(n=n())

df <- df %>% select(ageYear,Sex,Q_BMI,qnames) %>% drop_na() %>%
      mutate(Sex=ifelse(Sex=='F',0,1))#,n_risk_gps=as.integer(n_risk_gps)-1)
df %>% group_by(Q_DIAG_CKD_LEVEL) %>% summarise(n=n())


for (qname in qnames){
  centers <- get_kmeans(df %>% filter(!!as.name(qname)>0) %>% select(-qnames))
  print (qname)
  #print (nrow(df %>% filter(!!as.name(qname)>0) %>% select(-qnames)))
  print (centers)
}

df %>% group_by(Sex) %>% summarise(n=n())
mean(df$Q_BMI)


temp <- df
temp$clusters <- as.factor(clusters) 
temp%>% select(ageYear,Sex,clusters)
