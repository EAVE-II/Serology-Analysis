library(dplyr)
library(tibble)
library("tidyr")
library(ggplot2)

df_serology <- readRDS("/conf/EAVE/GPanalysis/data/serology_primcare_march22.rds") %>% as_tibble()

df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") %>% 
             as_tibble() 

df_qcovid <- df_qcovid %>% dplyr::filter(EAVE_LINKNO %in% df_serology$EAVE_LINKNO)

print (nrow(qcovid))

df_nrisks <- df_qcovid[,c('EAVE_LINKNO','n_risk_gps')]

df_ana <- df_serology %>% dplyr::left_join(df_nrisks) %>% 
          dplyr::select(c("EAVE_LINKNO","Sampledate_iso","test_result_qual","test_result_quant","n_risk_gps")) %>%
          dplyr::mutate(IgG = tidyr::extract_numeric(test_result_quant))  %>%
          dplyr::mutate(n_risk_gps = ifelse(is.na(n_risk_gps), "0", levels(n_risk_gps)[n_risk_gps]))
          

df_ana_positive <- df_ana %>% dplyr::filter(test_result_qual == "Positive")
df_ana_positive


bw <- 100
n <- nrow(df_ana_positive)
p <- ggplot(df_ana_positive, aes(x=IgG, fill=n_risk_gps)) +
     geom_histogram(binwidth = bw) + xlim(0,1000+bw)

#y=Quant_Result,position="identity",binwidth = bw) #aes(y=..density..)
#geom_histogram(position="identity",binwidth = bw) +
#scale_y_log10(
#breaks = seq(0, 0.05, 0.01) * (bw * n),
#labels = function(x) x / (bw * nrow(df_ana))
#)
p

ggplot_build(p)$data[[1]]

write.csv(ggplot_build(p)$data[[1]],'hist_data.csv')
