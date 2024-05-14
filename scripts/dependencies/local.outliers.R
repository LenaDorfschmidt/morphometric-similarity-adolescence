local.outliers = function(MM,MRI_table, sm, sm.long, outstr, nm.ggseg){
  nsub=nrow(MRI_table)
  nmorph = length(sm)
  nroi = dim(MM)[1]
  
  dir.create(outstr, recursive = T)
  # Local Outliers
  print('Estimating local outliers for each phenotype.')
  # calculate global number of median absolute deviations (MAD) per MM
  loc.mad = loc.mad.sign = loc.z.sign =coor =array(NA, dim=c(nmorph, nroi, nsub))
  med= MAD= SD = mu = array(NA, dim=c(nmorph, nroi))
  for (m in 1:nmorph) {
    for (r in 1:nroi){
      MAD[m,r]=mad(MM[r,m,], na.rm=T) # 29 Vals = 0 
      coor[m,r,]=(MM[r,m,]-median(MM[r,m,])) # 90000 Vals = 0 
      loc.z.sign[m,r,] = (MM[r,m,]-mean(MM[r,m,]))/sd(MM[r,m,])
      loc.mad.sign[m,r,] = (MM[r,m,]-median(MM[r,m,]))/mad(MM[r,m,], na.rm=T)
    }
  }
  loc.outl = abs(loc.mad.sign)
  
  coor.sum = apply(coor==0, c(1,2),sum) # as.vector(t(coor.sum))[(nroi+1):(nroi*2)] == coor.sum[2,] | as.vector(t(MAD))[(nroi+1):(nroi*2)]==MAD[2,]
  df.mad.issues = data.frame(med=as.vector(t(coor.sum)), mad= as.numeric(as.vector(t(MAD))==0),mm=rep(sm,each=nroi), label = rep(nm.ggseg, nmorph))
  df.mad.issues$med[df.mad.issues$med<=1]=NA; df.mad.issues$mad[df.mad.issues$mad==0]=NA
  
  sub.outl = apply(loc.outl, c(1,3), function(x) sum(x>=5, na.rm=T))
  df.outl.subj =data.frame(N = colSums(sub.outl, na.rm=T), MAD = abs(colSums(sub.outl, na.rm=T)-median(colSums(sub.outl, na.rm=T)))/mad(colSums(sub.outl, na.rm=T)))
  
  del.roi = which(colSums((MAD[-which(rowSums(MAD==0)>2),]==0))>0)
  if(any(del.roi)){
    nm.ggseg=nm.ggseg[-del.roi]
    loc.outl = loc.outl[,-del.roi,]
    MM = MM[-del.roi,,]
    nroi=nroi-length(del.roi)
  }
 
  print('Plotting local outliers for each phenotype.')
  p.local.outl = list()
  for (m in 1:nmorph) {
    df.local.outl = data.frame(value = rowSums(loc.outl[m,,]>=3)/nsub, label = nm.ggseg)
    p.local.outl[[m]] = ggseg(.data=df.local.outl, atlas=glasser,mapping=aes(fill=value), position='stacked',hemisphere =c('left','right'))+
      scale_fill_gradientn(colours = paletteer_c("viridis::plasma",n=100), limits=c(0,0.05), na.value='transparent', oob=squish,breaks =c(0,0.025,0.05))+
      ggtitle(sm.long[m])+
      ylab('')+xlab('')+
      labs(fill='#Subjects MAD>=3')+theme_minimal(base_family = "Gill Sans")+theme(legend.position = 'bottom')+
      guides(fill = guide_colorbar(title.position = "top", hjust=0.5))
  }
  
  g.loc<-ggarrange(plotlist = p.local.outl,
                   nrow=3, ncol=4, common.legend = TRUE, legend = 'bottom')
  ggsave(paste0(outstr,'/local.MAD.png'),g.loc, width=9, height=5, bg='white')
  
  # # PLot heatmap of Subject x ROI MAD values 
  # p.local.outl.heatmap = p.local.outl.roi.hist = p.local.outl.sub.hist = list()
  # lab.vals = c('0',paste0(seq(1,25,5),'-',seq(5,25,5)),'>25')
  # print('Plotting SUBJECT X REGIONS outliers as heatmap.')
  # for (m in 1:nmorph) {
  #   roi.by.sub.z = matrix(.bincode(loc.outl[m,,],c(3,4,5,ifelse(max(loc.outl[m,,])>5,max(loc.outl[m,,]) ,5))), nrow = nroi, ncol=nsub)
  #   df.loc.heatmap = data.frame(roi = rep(seq(1,nroi), nsub), sub = rep(seq(1:nsub),each=nroi),
  #                               outl = .bincode(loc.outl[m,,],c(3,4,5,ifelse(max(loc.outl[m,,])>5,max(loc.outl[m,,]) ,5))), age = rep(MRI_table$age,each=nroi))
  #   
  #   p.local.outl.heatmap[[m]] = ggplot(df.loc.heatmap, aes(x=roi, y = sub, fill=as.factor(outl)))+geom_tile()+xlab('Regions')+ylab('Subjects')+
  #     labs(fill='Z')+theme_minimal()+
  #     scale_fill_manual(values = c('yellow','orange','red'), na.value = 'transparent',labels=c(3,4,5))+ggtitle(sm.long[m])
  #   
  #   outl.counts.roi = df.loc.heatmap %>% group_by(roi) %>% dplyr::summarize(outl3=ceiling(sum(outl==1, na.rm=T)/5),outl4=ceiling(sum(outl==2, na.rm=T)/5),outl5=ceiling(sum(outl==3, na.rm=T)/5))
  #   outl.counts.sub = df.loc.heatmap %>% group_by(sub) %>% dplyr::summarize(outl3=ceiling(sum(outl==1, na.rm=T)/5),outl4=ceiling(sum(outl==2, na.rm=T)/5),outl5=ceiling(sum(outl==3, na.rm=T)/5))
  #   outl.counts.roi = gather(outl.counts.roi, outl, count, outl3:outl5, factor_key = T)
  #   outl.counts.sub = gather(outl.counts.sub, outl, count, outl3:outl5, factor_key = T)
  #   outl.counts.roi$count[outl.counts.roi$count>5] = 6
  #   outl.counts.sub$count[outl.counts.sub$count>5] = 6
  #   
  #   p.local.outl.roi.hist[[m]] <- ggplot(outl.counts.roi, aes(alpha=ifelse(count==0,0.2,1),fill = as.factor(outl), x = factor(count, labels=lab.vals[sort(unique(count)+1)])))+  
  #     geom_bar(stat='count', color='grey20')+ggtitle(sm.long[m])+
  #     scale_fill_manual(values = c('yellow','orange','red'), na.value = 'transparent',labels=c(3,4,5))+
  #     theme_minimal(base_family = "Gill Sans")+xlab('# Regions') + ylab('# Subjects')+
  #     scale_alpha(guide = 'none')+labs(fill='# outl')+
  #     theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))#+scale_x_continuous(breaks = seq(0,mod(max(outl.counts.roi$count),5)*5,5))
  #   
  #   p.local.outl.sub.hist[[m]] <- ggplot(outl.counts.sub, aes(alpha=ifelse(count==0,0.2,1),fill = as.factor(outl), x = factor(count, labels=lab.vals[sort(unique(count)+1)])))+  
  #     geom_bar(stat='count', color='grey20')+ggtitle(sm.long[m])+
  #     scale_fill_manual(values = c('yellow','orange','red'), na.value = 'transparent',labels=c(3,4,5))+
  #     theme_minimal(base_family = "Gill Sans")+ylab('# Regions') + xlab('# Subjects')+
  #     scale_alpha(guide = 'none')+labs(fill='# Z')+
  #     theme(legend.position = 'bottom',axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))
  # }
  # g.loc.heatmap<-ggarrange(plotlist = p.local.outl.heatmap,
  #                          nrow=3, ncol=4, common.legend = T, legend = 'bottom')
  # ggsave(paste0(outstr,'/local.heatmap.Z.png'),g.loc.heatmap, width=9, height=6)
  # 
  # g.loc.hist.roi.outl<-ggarrange(plotlist = p.local.outl.roi.hist,
  #                                nrow=3, ncol=4,common.legend = T, legend = 'bottom')
  # ggsave(paste0(outstr,'/local.hist.roi.Z.png'),g.loc.hist.roi.outl, width=9, height=6)
  # 
  # g.loc.hist.sub.outl<-ggarrange(plotlist = p.local.outl.sub.hist, 
  #                                nrow=3, ncol=4, common.legend = T, legend = 'bottom')
  # ggsave(paste0(outstr,'local.hist.sub.Z.png'),g.loc.hist.sub.outl, width=9, height=6)
  # 
  # 
  return(list(del.roi=del.roi, loc.outl=loc.outl))
}