plotglobal=function(MRI_table,newdat,sm.long,sm.str){
  newdat <- expand.grid(age=c(min(MRI_table$age),
                              max(MRI_table$age)),
                        sex = unique(MRI_table$sex))
  
  df=MRI_table[,c(c('id','age','sex',sm.str))]
  colnames(df)[grep(sm.str,colnames(df))] = 'cur'
  lm.glob = lm(cur~age+sex, data = df)
  lm.pred = predict.lm(lm.glob, newdata=newdat, se.fit = TRUE)
  plotob=ggplot(df,aes(x=age, y=cur, color=sex))+
    geom_point()+
    ggtitle(sm.long)+
    geom_line(data=newdat, aes(y=lm.pred$fit),size=1)+
    geom_ribbon(data=newdat[newdat$sex=='Female',],aes(x=age,
                                                     min=lm.pred$fit[newdat$sex=='Female']-2*lm.pred$se.fit[newdat$sex=='Female'],
                                                     ymax=lm.pred$fit[newdat$sex=='Female']+2*lm.pred$se.fit[newdat$sex=='Female']),
                alpha=0.5, inherit.aes=F) +
    geom_ribbon(data=newdat[newdat$sex=='Male',],aes(x=age,
                                                     ymin=lm.pred$fit[newdat$sex=='Male']-2*lm.pred$se.fit[newdat$sex=='Male'],
                                                     ymax=lm.pred$fit[newdat$sex=='Male']+2*lm.pred$se.fit[newdat$sex=='Male']),
                alpha=0.5, inherit.aes=F) +
    ylab('')+theme_minimal(base_family = "Gill Sans")+theme(plot.title = element_text(hjust = 0.5))+labs(color='')+
    scale_color_manual(values=c(colors.morph[7],colors.morph[1]))
  
  return(plotob)
}

