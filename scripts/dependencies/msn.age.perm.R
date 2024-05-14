msn_age_perm = function(MRI_table, str){
  str.l.sl.t = str.l.sl = vector(length=nroi)
  roi.ranef = matrix(nrow = nroi, ncol = length(unique(MRI_table$id)))
  edf = vector(length = nroi)
  p=list()
  for (r in 1:nroi) {
    df.roi =  data.frame(id = MRI_table$id, center= MRI_table$site, age = MRI_table$age, sex = MRI_table$sex, str =str[r,])
    
    # Linear Models
    lm.roi = tryCatch({ 
      lm.roi <- lme(str~age+sex+center, random = ~1|id,data = df.roi,na.action = na.omit) }
      , error = function(e) {an.error.occured <<- TRUE})
    
    if(!(class(lm.roi)=='lme')){next}else{
      str.l.sl.t[r] = summary(lm.roi)$tTable[2,4]
      str.l.sl[r] = lm.roi$coefficients$fixed[2]
    }
  }
  return(list(t=str.l.sl.t, sl= str.l.sl))
}
