pkg.env <- new.env(parent = emptyenv())
pkg.env$dd <- do.call(rbind,jsonlite::fromJSON(system.file("extdata", "data_dictionary.json", package = "eavehelpers"),flatten=T)$columns)

var.names = c(
  "insufficient_response"="Insufficient Antibodies",
  "ch_resident"="Care Home Resident",
  "immuno"="Immunosuppressed",
  "shielding"="Shielding",
  "simd2020v2_sc_quintile"="Scottish Index of Multiple Deprivation",
  "cat"="Vaccine Dose",
  "n_risk_gps"="Number of Risks",
  "Sex"="Sex",
  "Q_BMI"="Body Mass Index",
  "prior_infection"="Known Prior Infection",
  "ur6_2016_name"="Location",
  "additional"="Subsequent Vaccinations",
  "days_since_start"="Pandemic Period",
  "days_since_vac2"="Days Since 2nd Vaccination", 
  "n_risk_gpsOther"="Other Risks",
  "n_risk_gps"="Number of Risk Groups",
  "age"="Age",
  "SexM"="Sex (Male)",
  "ageYear"="Age",
  #"Q_DIAG_DIABETES_1"="Diabetes (Pre & Types I&II)", 
  #"Q_DIAG_DIABETES_2"="Diabetes (Other)",
  "Q_DIAG_DIABETES_1"="Diabetes (Type-I)", 
  "Q_DIAG_DIABETES_2"="Diabetes (Type-II)",
  "Q_DIAG_SICKLE_CELL"="Sickle Cell Disease",
  "Q_DIAG_RESP_CANCER"="Respiratory Cancer",
  "Q_DIAG_ASTHMA"="Asthma",    
  "Q_DIAG_BLOOD_CANCER"="Haematological Cancer",
  "Q_DIAG_CHD"="Coronary Heart Disease ",
  "Q_DIAG_COPD"="Chronic Obstructive Pulmonary Disease",   
  "Q_DIAG_CKD_LEVEL"="Chronic Kidney Disease",
  "Q_DIAG_CKD"="Chronic Kidney Disease",
  "Q_DIAG_CKD3"="Chronic Kidney Disease (Level 3)",
  "Q_DIAG_CKD4"="Chronic Kidney Disease (Level 4)",
  "Q_DIAG_CKD5"="Chronic Kidney Disease (Level 5)",
  "Q_DIAG_AF"="Atrial Fibrillation",
  "Q_DIAG_CCF"="Heart Failure",
  "Q_DIAG_EPILEPSY"="Epilepsy",     
  "Q_DIAG_FRACTURE"="A prior fracture of hip, wrist, spine or humerus",     
  "Q_DIAG_IMMU"="Immune Deficiency",         
  "Q_DIAG_NEURO"="Rare Neurone Disease",       
  "Q_DIAG_PARKINSONS"="Parkinsons",  
  "Q_DIAG_PULM_RARE"="Cystic Fibrosis or Bronchiectasis or Alveolitis",   
  "Q_DIAG_PVD"="Peripheral Vascular Disease",         
  "Q_DIAG_RA_SLE"="Rheumatoid Arthritis",      
  "Q_DIAG_SEV_MENT_ILL"="Severe Mental Health Illness",
  "Q_DIAG_STROKE"="Stroke",      
  "Q_DIAG_SICKLE_CELL"="Sickle Cell Disease",
  "Q_DIAG_VTE"="Thrombosis or Pulmonary Embolus",
  "Q_DIAG_CEREBRALPALSY"="Cerebral Palsy",
  "Q_DIAG_CIRRHOSIS"="Cirrhosis",
  "Q_DIAG_CONGEN_HD"="Congenital Heart Disease",
  "Q_DIAG_HIV_AIDS"="HIV/Aids",
  "Q_DIAG_PULM_HYPER"="Pulmonary Hypertension",
  "Q_DIAG_DEMENTIA"="Dementia"
  )

pkg.env$var.names = var.names


#' Get the Label
#'
#' Pass a variable name (as a string) and get the label(s) found in
#' the data dictionary returned
#' 
#' @param the variable name
#' @return a label (or multiple, if multiple found) of the variable
#' @export
get_label <- function(var){
    return (as.character(pkg.env$var.names[var]));
    #return (pkg.env$dd[pkg.env$dd$name == var,]$label);
}

#' Get the Description
#' 
#'
#' Pass a variable name (as a string) and get the description(s) found in
#' the data dictionary returned
#' 
#' @param the variable name
#' @return a description (or multiple, if multiple found) of the variable
#' @export
get_description <- function(var){
    return (pkg.env$dd[pkg.env$dd$name == var,]$description);
}

#' Look up a variable in the data dictionary
#'
#' Return all information found about this variable by name
#' 
#' @param the variable name
#' @return a data.frame of information found about the variable input
#' @export
get_var <- function(var){
    return (pkg.env$dd[pkg.env$dd$name == var,]);
}

#' Get the Data Dictionary
#'
#' Return the full data dictionary
#'
#' @return the full EAVE data dictionary as a data.frame
#' @export
get_data_dictionary <- function(){
    return (pkg.env$dd);
}