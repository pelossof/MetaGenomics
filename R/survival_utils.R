library(survival)

survu.plotkm = function(time, death, groups=c(), col=c(), xlab = 'Time (Months)', ylab = 'Survival Probability', noise = runif(length(time)), add=F, ...) {
## time: named numeric vector of time
## death: named boolean vector of death event
  
  fit = list()

  #remove patients without survival or followup data
#  is.no.data = time == 0 & death == 0
#  time = time[!is.no.data]
#  death = death[!is.no.data]
  time = time + noise
  groups = lapply(groups, intersect, names(time))

  if(length(groups) == 0){
    groups = list('all'=names(time))
  }
  
  if(is.null(names(groups))){
    names = as.character(seq(1:length(groups)))
  }else{
    names=names(groups)
  }
  if(length(col) == 0){
    col = as.character(seq(1:length(groups)))
  }
  
  if(add){
    for (i in 1:length(groups)){
      fit[[names(groups)[i]]] = survfit(Surv(time[unlist(groups[[i]])], death[unlist(groups[[i]])])~1, conf.type="none")
      lines(fit[[names(groups)[i]]], xlab = xlab, ylab = ylab, col=col[i])
    }
  }else{
    for (i in 1:length(groups)){
      fit[[names(groups)[i]]] = survfit(Surv(time[unlist(groups[[i]])], death[unlist(groups[[i]])])~1, conf.type="none")
      if(i == 1){
        if ('xlim' %in% names(list(...)))
          plot(fit[[names(groups)[i]]], xlab = xlab, ylab = ylab, col=col[i], ...)
        else
          plot(fit[[names(groups)[i]]], xlab = xlab, ylab = ylab, col=col[i], xlim=c(0,max(time[unlist(groups)] + 1, na.rm=T)), ...)
      }else{
        lines(fit[[names(groups)[i]]], xlab = xlab, ylab = ylab, col=col[i])
      }
    }
    
    if(length(groups) > 1){
      labels = unlist(lapply(1:length(groups), function(x){rep(x, length(groups[[x]]))}))
      diff = survdiff(Surv(time[unlist(groups)], death[unlist(groups)]) ~ labels)
      p.val = 1 - pchisq(diff$chisq, length(groups)-1)
      legend("bottomleft", c(paste(names, " N=", unlist(lapply(groups,length)),sep=""), 
                             paste("p-val", sprintf("%.2e", p.val))), col=c(col, rgb(1,1,1,1)), pch=16, bty='n', cex=0.8)
    }
  }

  return (fit)
}



survu.plotkm.old = function(data, groups=c(), col=c(), add=F, ...) {
  time = data$patient$days_to_last_followup
  time[is.na(time)] = data$patient$days_to_death[is.na(time)]
  names(time) = data$patient$bcr_patient_barcode
  event = as.numeric(is.na(data$patient$days_to_last_followup))
  names(event) = data$patient$bcr_patient_barcode
  
  #remove patients without survival or followup data
  is.no.data = time == 0 & event == 0
  time = time[!is.no.data]
  event = event[!is.no.data]
  time = time + runif(length(time))
  groups = lapply(groups, intersect, names(time))

  if(length(groups) == 0){
    groups = list('a'=names(time))
  }
  
  if(is.null(names(groups))){
    names = as.character(seq(1:length(groups)))
  }else{
    names=names(groups)
  }
  if(length(col) == 0){
    col = as.character(seq(1:length(groups)))
  }
  
  if(add){
    for (i in 1:length(groups)){
      test = survfit(Surv(time[unlist(groups[[i]])]/30.43, event[unlist(groups[[i]])])~1, conf.type="none")
      lines(test, xlab="Time (months)", ylab="Survival Probability", col=col[i])
    }
  }else{
    for (i in 1:length(groups)){
      test = survfit(Surv(time[unlist(groups[[i]])]/30.43, event[unlist(groups[[i]])])~1, conf.type="none")
      if(i == 1){
        plot(test, xlab="Time (months)", ylab="Survival Probability", col=col[i], xlim=c(0,max(time[unlist(groups)]/30.43 + 1, na.rm=T)), ...)
      }else{
        lines(test, xlab="Time (months)", ylab="Survival Probability", col=col[i])
      }
    }
    
    if(length(groups) > 1){
      labels = unlist(lapply(1:length(groups), function(x){rep(x, length(groups[[x]]))}))
      diff = survdiff(Surv(time[unlist(groups)], event[unlist(groups)]) ~ labels)
      p.val = 1 - pchisq(diff$chisq, length(groups)-1)
      legend("bottomleft", c(paste(names, " N=", unlist(lapply(groups,length)),sep=""), 
                             paste("p-val", sprintf("%.2e", p.val))), col=c(col, rgb(1,1,1,1)), pch=16, bty='n', cex=0.8)
    }
  }
}







