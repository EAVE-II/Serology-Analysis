

pc.df_ana <- df_serology_pc %>% left_join(df_vac) %>% left_join(df_qcovid) %>% 
  mutate(ageYear = as.factor(cut(age, breaks = c(0,20,40,60,150), right = T, 
                                 labels = c("0-19","20-39","40-59","60+")))) 

qnames <- c("Q_DIAG_BLOOD_CANCER","Q_DIAG_DIABETES_1","Q_DIAG_CKD_LEVEL")

pc.model1 <- pc.df_ana  %>%
              filter(as.Date(Sampledate_iso) < d2_datetime & as.Date(Sampledate_iso) > d1_datetime) %>%
              mutate(days=as.numeric(as.Date(Sampledate_iso) - d1_datetime,units='days')) %>%
              select(days,IgG,ageYear,qnames)

pc.model2 <- pc.df_ana %>% 
  filter(as.Date(Sampledate_iso) < d3_datetime & as.Date(Sampledate_iso) > d2_datetime) %>%
  mutate(days=as.numeric(as.Date(Sampledate_iso) - d2_datetime,units='days')) %>%
  select(days,IgG,ageYear,qnames)

pc.model3 <- pc.df_ana %>% 
  filter(as.Date(Sampledate_iso) < d4_datetime & as.Date(Sampledate_iso) > d3_datetime) %>%
  mutate(days=as.numeric(as.Date(Sampledate_iso) - d3_datetime,units='days')) %>%
  select(days,IgG,ageYear,qnames)


pc.model <- pc.model1 %>% mutate(stage='1') %>% rbind(pc.model2 %>% mutate(stage='2'))%>% rbind(pc.model3 %>% mutate(stage='3'))
#model <- model2 %>% mutate(stage='2') %>% rbind(model3 %>% mutate(stage='3'))

pc.model %>% group_by(Q_DIAG_CKD_LEVEL) %>% summarise(n=n())

pc.model2 <- pc.model %>% mutate(None=Q_DIAG_DIABETES_1+Q_DIAG_BLOOD_CANCER+Q_DIAG_CKD_LEVEL) %>% 
             pivot_longer(cols = c(qnames,None), names_to = 'cond', values_to = 'Values') %>% 
             filter((Values>0 & cond!='None') | (Values==0 & cond=='None')) %>% select(-Values) 
             

pc.model2

p <- pc.model2 %>% filter(days<80 & !is.na(ageYear)) %>% 
  ggplot(aes(x=days,y=IgG)) +
  stat_summary_bin(aes(color=cond),bins=15) +
  #stat_smooth(method="gam", aes(color=cond),formula=y ~ splines::bs(x,df=6),se=T) +
  scale_colour_discrete_phs(palette='all') +
  #ylim(0,2000) + 
  #scale_y_log10() +
  labs(x='Days Since Vaccination',y='IgG',color='Risks')

p + facet_grid(stage ~ cond,scales="free")#,labeller=names)



bd.df_ana <- df_serology_bd %>% left_join(df_vac) %>%
  mutate(ageYear = as.factor(cut(age, breaks = c(0,20,40,60,150), right = T, 
                                 labels = c("0-19","20-39","40-59","60+")))) 

bd.model1 <- bd.df_ana  %>% 
  filter(as.Date(Sampledate_iso) < d2_datetime & as.Date(Sampledate_iso) > d1_datetime) %>%
  mutate(days=as.numeric(as.Date(Sampledate_iso) - d1_datetime,units='days')) %>%
  select(days,IgG)

bd.model2 <- bd.df_ana %>% 
  filter(as.Date(Sampledate_iso) < d3_datetime & as.Date(Sampledate_iso) > d2_datetime) %>%
  mutate(days=as.numeric(as.Date(Sampledate_iso) - d2_datetime,units='days')) %>%
  select(days,IgG)

bd.model3 <- bd.df_ana %>% 
  filter(as.Date(Sampledate_iso) < d4_datetime & as.Date(Sampledate_iso) > d3_datetime) %>%
  mutate(days=as.numeric(as.Date(Sampledate_iso) - d3_datetime,units='days')) %>%
  select(days,IgG)

bd.model <- bd.model1 %>% mutate(stage='1') %>% rbind(bd.model2 %>% mutate(stage='2'))%>% rbind(bd.model3 %>% mutate(stage='3'))


#model <- model %>% mutate(x=cut(days,breaks=c(0,5,10,15,20,25,30,35,45,50,55,60,65,70,75,80,85))) %>%
#  filter(!is.na(x))


pc.model <- pc.model %>% mutate(nIgG = IgG/max(IgG))
bd.model <- bd.model %>% filter(!is.na(IgG)) %>% mutate(nIgG = IgG/max(IgG))

model <- pc.model %>% mutate(dataset='Primary Care') %>% rbind(bd.model %>% mutate(dataset='Blood Donors'))



#saveRDS(bd.df_ana,"/conf/EAVE/GPanalysis/data/temp/time_analysis_serology_snbts.rds")
#saveRDS(pc.df_ana,"/conf/EAVE/GPanalysis/data/temp/time_analysis_serology_primcare.rds")




p <- model %>% filter(days<150) %>% 
              ggplot(aes(x=days,y=nIgG)) +
              stat_summary_bin(aes(color=dataset)) +
              stat_smooth(method="gam", aes(color=dataset),formula=y ~ splines::bs(x,df=6),se=T) +
              scale_colour_discrete_phs(palette='all') +
              ylim(0,1) + 
              #scale_y_log10() +
              labs(x='Days Since Vaccination',y='Normalised IgG',color='Dataset')

names <- as_labeller(
  c(`1`='Stage 1',`2`='Stage 2',`3`='Stage 3')
  )
p + facet_grid(stage ~ .,scales="free",labeller=names)



              
