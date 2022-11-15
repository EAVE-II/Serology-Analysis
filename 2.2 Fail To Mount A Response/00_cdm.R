
library(dplyr)
library(tidyr)
library(jsonlite)
library(stringr)

generate_report = function(input_file,out_folder){
  #if (!dir.exists(sub_folder)){
  #  dir.create(sub_folder, recursive = T)
  #} 
  df <- readRDS(input_file) 
  if (nrow(df)>10000){
    df <- df %>% sample_n(10000)
  }
  
  data = list()
  for(col in colnames(df)){
    type <- pillar::type_sum(df[[col]])
    #print(type)
    if (type != 'chr'){
      max <- max(df[[col]])
      min <- min(df[[col]])
    }
    else {
      min <- ''
      max <- ''
    }
    
    summary <- df %>% group_by(!!sym(col)) %>% summarise(n=n()) %>% 
      filter(n>10) %>%
      mutate(frac=n/sum(n)) %>%
      arrange(desc(frac)) %>%
      head(5) %>%
      mutate(value=!!sym(col),frac=signif(frac,2)) %>%
      select(value,frac)
    
    if (nrow(summary) < 1){
      summary[1,] <- NA
    }else if(col == "EAVE_LINKNO"){
      summary <- summary %>% head(1)
      summary[1,] <- NA
    }
    summary <- tryCatch(summary %>% mutate(across(everything(), as.character)), error = function(e) NULL)
    if(is.null(summary)){
      next 
    }
    res <- summary#toJSON(summary,pretty=TRUE,na="null")
    data[[col]] = list(name=col,label='',description='',type=type,min=min,max=max,values=res)#,pretty=TRUE)
    
  }
  #col <- gsub("/", "_", col)
  fname <- paste0(out_folder,tools::file_path_sans_ext(basename(input_file)),'.RData')
  print (fname)
  saveRDS(data,fname)
}

generate_scan_report = function(input_file,out_folder){
  print (input_file)
  
  sub_folder <- paste0(out_folder,tools::file_path_sans_ext(basename(input_file)))
  if (!dir.exists(sub_folder)){
    dir.create(sub_folder, recursive = T)
  } 
  df <- readRDS(input_file) 
  if (nrow(df)>10000){
    df <- df %>% sample_n(10000)
  }
  print (colnames(df))
  for (col in colnames(df)){
    summary <- df %>% group_by(!!sym(col)) %>% summarise(n=n()) %>% 
      filter(n>10) %>%
      mutate(frac=n/sum(n)) %>%
      arrange(desc(frac)) %>%
      head(10) %>%
      mutate(value=!!sym(col)) %>%
      select(value,frac)
    
    if (nrow(summary) < 1){
      summary[1,] <- NA
    }else if(col == "EAVE_LINKNO"){
      summary <- summary %>% head(1)
      summary[1,] <- NA
    }
    summary$table <- basename(input_file)
    summary$field <- col
    summary$folder <- dirname(input_file)
    summary$date_last_modified <- file.mtime(input_file)
    
    print (summary)
    summary <- tryCatch(summary %>% mutate(across(everything(), as.character)), error = function(e) NULL)
    if(is.null(summary)){
      next 
    }
    res <- toJSON(summary,pretty=TRUE,na="null");
    col <- gsub("/", "_", col)
    fname <- paste0(sub_folder,"/",col,".json")
    write(res,fname)
    
  }
  rm(df)
}

out_folder <- "/home/calumm09/DataDictionary/"
in_folder <- "/conf/EAVE/GPanalysis/data/"
files <- list.files(path=in_folder, pattern="*.rds", full.names=T)
files


#generate_report("/conf/EAVE/GPanalysis/data/serology_snbts_july22_v3.rds",out_folder)
#data = readRDS('/home/calumm09/DataDictionary/lookups/serology_snbts_july22_v3.RData')


for (fname in files){
  generate_scan_report(fname,out_folder)
}





generate_scan_report("/conf/EAVE/GPanalysis/data/serology_primcare_march22.rds",out_folder)

rules <- '{
            "measurement_id":{"field":"LabSpecimenNo"}, 
            "person_id":"EAVE_LINKNO",
            "measurement_source_value":"test_result_qual",
            "value_as_number":"test_result_quant"
}'

mapping <- fromJSON(rules)


dataframe_to_cdm = function(df,rules){
  temp <- NULL
  for (mapping in rules){
      retval <- df
      for ( col in labels(mapping)){
        var <- mapping[[col]]
        print (var)
        if (is.list(var)){
            field <- var$field
            term_mapping <- var$term_mapping
            if(!is.null(term_mapping)){
              for (term in labels(term_mapping)){
                output <- term_mapping[[term]]
                retval <- retval %>% mutate(!!sym(col) := ifelse(!!sym(field)==term,output,NA))
              }
            }
            else{
              retval <- retval %>% mutate(!!sym(col) := !!sym(field))
            }
        }#is.list
        else {
            retval <- retval %>% mutate(!!sym(col) := var)
        }
      }#end of for loop over output columns
      retval <- retval %>% select(labels(mapping)) %>%
                 filter(if_any(contains('_concept_id'),~!is.na(.)))
      print (retval)
      temp <- rbind(temp,retval)
  }
  return (temp)
}


df_measurements <- df_serology_pc %>% dataframe_to_cdm(mapping)
df_measurements

df_qcovid <- readRDS("/conf/EAVE/GPanalysis/data/cleaned_data/QCOVID_feb22.rds") 
df_qcovid

df_temp <- df_qcovid %>% head(10)
df_temp

qnames <- df_temp %>% select(contains("Q_DIAG")) %>% colnames 
qnames


template <- '{
            "person_id":{"field":"EAVE_LINKNO"},
            "condition_concept_id":{
                                    "field":"${field}",
                                    "term_mapping":{
                                          "1":"${field}"
                                    }
            },  
            "condition_source_value":{"field":"${field}"}}'

rules <- list()
for (name in qnames){
  rules <- append(rules,str_interp(template,list(field=name)))
}
rules <- paste0('[',paste0(rules,collapse=','),']')
rules

mapping <- fromJSON(rules,simplifyVector = FALSE)


df_conditions <- df_qcovid %>% dataframe_to_cdm(mapping)
df_conditions
df_conditions %>% group_by(condition_concept_id) %>% summarise(n=n())
