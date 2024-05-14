global.outliers = function(MM,MRI_table,sm,sm.group, outstr){
  dir.create(outstr)
  nsub = nrow(MRI_table)
  nmorph = dim(MM)[2]
    
  # Correlation between morphometric measures
  print('Estimating correlation between phenotypes')
  MM.cor.mat = rcorr(t(apply(MM, c(2,3), mean)))
  mat=MM.cor.mat$r
  rownames(mat) = sm; colnames(mat) = sm 

  print('Estimating global morphometric effects')
  # Global Morphometric Effects
  glob.morph = matrix(ncol = nmorph, nrow = nsub)
  for (m in 1:nmorph) {
    print(m)
    if (sm[m] == 'SA' | sm[m] == 'GM') {
      glob.morph[,m] = colSums(MM[,m,])
    }else{
      glob.morph[,m] = colMeans(MM[,m,])
    }
  }
  colnames(glob.morph) = paste0('glob.',sm)

  # Outlier detection
  # loop over morphological variables
  # calculate global median absolute deviations (MAD)/ Z Score
  glob.z.sign = glob.mad.sign=  array(NA, dim=c(nmorph, nsub))
  for (m in 1:nmorph) {
    glob.mad.sign[m,] = (glob.morph[,m]-median(glob.morph[,m],na.rm=T))/mad(glob.morph[,m],na.rm=T)
    glob.z.sign[m,] = (glob.morph[,m]-mean(glob.morph[,m],na.rm=T))/sd(glob.morph[,m],na.rm=T)
  }
  glob.z = abs(glob.z.sign)
  glob.mad = abs(glob.mad.sign)
  
  # Plot global morphometric development by age; color in outliers as MAD>c(3,4,5)
  p.glob.outl = list(); i=1
  for (outlm in list(glob.z,glob.mad)){
    for (m in 1:nmorph) {
      df.glob.outl = data.frame(outl.glob = glob.morph[,m], mm = rep(sm[m],nsub), age = MRI_table$age, 
                                col = as.factor(.bincode(outlm[m,],c(min(outlm[m,])-1,3,4,5,ifelse(max(outlm[m,])>5,max(outlm[m,]) ,5)))))
      p.glob.outl[[m]] <- ggplot(df.glob.outl, aes(x=age, y = outl.glob, col=col))+
        geom_point(alpha=0.3)+
        scale_color_manual(values = c('grey20','yellow','orange','red'), na.value = 'transparent',labels=c(NA,3,4,5))+
        labs(col= ifelse(i==1,'Z-Score','MAD'))+
        ggtitle(sm.long[m])+
        ylab(paste0('Global ',sm[m]))+ xlab('Age')+
        theme_minimal(base_family = "Gill Sans")+theme(legend.position = 'bottom')
    }
    g<-ggarrange(plotlist=p.glob.outl, nrow=3, ncol=4, common.legend = TRUE, legend = 'bottom')
    ggsave(ifelse(i==1,paste0(outstr,"global.Z.png"),paste0(outstr,"global.MAD.png")),g, width=10, height=5,bg = 'white')
    i=i+1
  }

  # Exclude subjects that are outliers:
  excl.sub.mad = unique(unlist(apply(glob.mad,1,function(x) which(x>5)))) # This would allow to excl based on individual mm
  print(paste0('Excluding ',length(excl.sub.mad),' scans Because of Global Outliers.'))

  return(list(glob.morph=glob.morph, 
              excl.sub.mad=excl.sub.mad))
}