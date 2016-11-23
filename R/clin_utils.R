
source('~/paths.R')
source(paths.utils('survival_utils.R'))
#source(paths.utils('firehose_utils.R'))
library(survival)

clinu.build <- function(clin.data){
  clin.pats <- as.character(clin.data['patient.bcr_patient_barcode',])
  dups.pats <- duplicated(clin.pats, fromLast=T)
  clin.data <- clin.data[,!dups.pats]
  
  ## prepare clin matrix
  days_to_last_followup = apply(clin.data[grep('patient.days_to_last_followup', rownames(clin.data)),],2,function(x){
    if(any(!is.na(x))){
      as.numeric(max(x[!is.na(x)]))
    }else{NA}
  })
  
  days_to_death=apply(clin.data[grep('days_to_death', rownames(clin.data)),],2,function(x){
    if(any(!is.na(x))){
      as.numeric(unique(x[!is.na(x)]))
    }else{NA}
  })
  days_to_death = sapply(days_to_death, max)
#  days_to_death=unlist(days_to_death)
  dupNames = names(days_to_last_followup[!is.na(days_to_last_followup)])[names(days_to_last_followup[!is.na(days_to_last_followup)]) %in% names(days_to_death[!is.na(days_to_death)])]
  days_to_last_followup[dupNames] = NA

  clin = list(patient = list(bcr_patient_barcode = toupper(as.character(clin.data['patient.bcr_patient_barcode',])),
                             days_to_last_followup = days_to_last_followup,
                             days_to_death=days_to_death,
                             gender=as.character(clin.data['patient.gender',]), 
                             anatomic_neoplasm_subdivision = as.character(clin.data['patient.anatomic_neoplasm_subdivision',]),
                             days_to_birth=as.numeric(clin.data['patient.days_to_birth',]),
                             number_of_first_degree_relatives_with_cancer_diagnosis = as.numeric(clin.data['patient.number_of_first_degree_relatives_with_cancer_diagnosis',]),
                             stage=as.character(clin.data['patient.stage_event.pathologic_stage',])),
              radiations = list(radiations=as.character(clin.data['patient.radiations',]),
                                radiation.bcr_radiation_barcode=as.character(clin.data['patient.radiations.radiation.bcr_radiation_barcode',])),
              drugs = list(drugs=as.character(clin.data['patient.drugs',]),
                           drug.bcr_drug_barcode=as.character(clin.data['patient.drugs.drug.bcr_drug_barcode',])))          
  return(clin)
}



cu.table1 = function(df, study = 'study') {
# for example:
#  dat = data.frame(study = c(rep('MSK', 10), rep('TM', 8)), gender = c(rep('M', 8), rep('F', 7), rep(NA, 3)), stage = factor(sample(1:4, 18, replace=T)), kras = factor(runif(18)>0.5))

  r = list()
  p = list()
  for( s in unique(df[, study])) {
    r[[s]] = apply( subset(df, df[, study] == s, ), 2, table, useNA='always') 
    p[[s]] = lapply (r[[s]], function (x) return (100*x/sum(x)))
  }

  clin.tab = data.frame(count = unlist(r), percent = unlist(p))
  return (clin.tab)
}

