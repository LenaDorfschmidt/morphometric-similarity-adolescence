#### Load required packages
library(tidyr)        # data wrangling

library(ggplot2)      # plotting
library(ggseg)        # brain surface plotting
library(ggsegGlasser) # brain surface plotting
library(paletteer)    # color scales
library(ggpubr)       # arrange plots
library(Hmisc)        # compute correlation on matrix rows 
library(scales)       # squish function
library(xtable)       # saving tables for latex doc

#### Load dependencies
source('./scripts/dependencies/glob.outliers.R')
source('./scripts/dependencies/local.outliers.R')
source('./scripts/external/rotate_parcellation-master/R/perm.sphere.p.R')

# Load Data
load(paste0('./data/hcpd.RData')) 
load('data/perm.id.360.RData')

##################################
#         ANALYSIS SETUP
##################################
# Set parameters here
outpath='results/SI/sensitivity-analyses/out-of-sample-replication/'
outpath.SI = paste0(outpath,'/SI/')
outpath.main = paste0(outpath,'/main/')

# Create output directory structure
dir.create(outpath, recursive = T)
dir.create(outpath.main)
dir.create(outpath.SI)

# Generate color scales
colors.morph = paletteer_d('rcartocolor::Tropic') #paletteer_d('Redmonder::dPBIYlBu')
colors.msn.age = paletteer_d('Redmonder::dPBIPuOr')

printout.stats=list() # To make our live easier we will print loads of info into this

# Subset for subjects that are at least 14 years old to match NSPN
choose.subs = which(MRI_table$age>=14)
SM = SM[,,choose.subs]
MRI_table = MRI_table[choose.subs,]

# Save region names for matching with NSPN
write.table(nm,file=paste0('results/SI/replication/hcp-d/nm.hcp.txt'),row.names = F)

##################################
#    GLOBAL OUTLIER DETECTION
##################################
glob.outl.output = global.outliers(SM,MRI_table,sm,sm.group, paste0(outpath.SI,'outliers/'))
excl.sub.mad = glob.outl.output$excl.sub.mad
MRI_table=cbind(MRI_table, glob.outl.output$glob.morph)

printout.stats[[1]] <- paste0('Number of subjects excluded due to global outliers: ', length(excl.sub.mad))

# Delete subjects based on local outliers
if(any(excl.sub.mad)){
  SM = SM[,,-excl.sub.mad]
  MRI_table = MRI_table[-excl.sub.mad,]
}

nsub = dim(MRI_table)[1]

##################################
#      LOCAL OUTLIER DETECTION
##################################
loc.outl.output = local.outliers(SM,MRI_table, sm, sm.long,paste0(outpath.SI,'outliers/'),nm.ggseg)
loc.outl = loc.outl.output$loc.outl; del.roi = loc.outl.output$del.roi

printout.stats[[2]] <- paste0('Number of regions excluded due to local outliers: ', length(del.roi))

# Delete ROIs based on local outliers
if(any(del.roi)){
  SM = SM[-del.roi,,]
  #loc.mad.sign = loc.mad.sign[,-del.roi,]
  nroi = dim(SM)[1]
  nm.ggseg = nm.ggseg[-del.roi]
  nmorph = length(sm)
  loc.outl = loc.outl[,-del.roi,]
  Yeo=Yeo[-del.roi]
  mesulam=mesulam[-del.roi]
  vonEconomo=vonEconomo[-del.roi]
}

# Set all outl>=4 to NA
loc.outl = aperm(loc.outl, c(2,1,3))
SM = replace(SM, loc.outl >= 5, NA)

printout.stats[[3]] <- paste0('Outliers percentage (number of ROIs set to NA across all subjects, regions, phenotypes): ',sum(loc.outl >= 5, na.rm=T)/length(loc.outl)*100,'%')

##################################
#    MORPHOMETRIC DEVELOPMENT
##################################
#morph.dev.output = morphometric.development(SM, MRI_table,sm.long, nm.ggseg,colors.morph, paste0(outpath.SI,'/morphometric-development/'))
#glob.morph=morph.dev.output$glob.morph; df.roi.morph.effects=morph.dev.output$df.roi.morph.effects; df.roi.morph.effects.glob.corr=morph.dev.output$df.roi.morph.effects.glob.corr
SM.param=aperm(apply(SM,c(3,2),scale),c(1,3,2)) # z-score and reorder dimensions

# Data Overview
# Acquisition overview
acq.overview=MRI_table %>% group_by(id) %>% mutate(firstscan=order(age), firstage=min(age))
relevel=acq.overview[acq.overview$firstscan==1,]$id[order(acq.overview[acq.overview$firstscan==1,]$age)]
acq.overview$id = factor(acq.overview$id, levels = relevel)
acq.overview$y = as.numeric(acq.overview$id)

# SI Figure 16, Panel A
plot.acq.overview = ggplot(acq.overview, aes(y=y,x=age, color=firstage))+
  geom_point()+
  labs(color='Baseline Age')+
  theme_bw()+ylab('Participant')+
  scale_color_gradientn(colours = paletteer_c('viridis::magma',100))
ggsave(plot=plot.acq.overview, filename=paste0(outpath.SI, 'acquisition.overview.png'),height=4,width=4)


##################################
#             FIGURE 1
################################## 
# Global plots
nmorph=length(sm)
newdat <- expand.grid(age=c(min(MRI_table$age),
                            max(MRI_table$age)),
                      sex = unique(MRI_table$sex))
glob.res=matrix(NA, nrow=nmorph, ncol=4); rownames(glob.res)=sm.long; colnames(glob.res)=c('t.age','p.age','t.sex','p.sex')
glob.res = data.frame(glob.res)
for (m in 1:nmorph) {
  sm.str=paste0('glob.',sm[m])
  df=MRI_table[,c(c('id','age','sex',sm.str))]
  colnames(df)[grep(sm.str,colnames(df))] = 'cur'
  lm.glob = lm(cur~age+sex+age*sex, data = df)
  glob.res[m,1:2]=summary(lm.glob)$coefficients[2,3:4];
  glob.res[m,3:4]=summary(lm.glob)$coefficients[3,3:4]; 
  #lm.pred = predictSE.lme(lm.glob, level=0, newdata=newdat)
}
glob.res$p.age.fdr = p.adjust(glob.res$p.age, method='fdr')
glob.res$p.sex.fdr = p.adjust(glob.res$p.sex, method='fdr')
print(xtable(glob.res, type = "latex"), file=paste0(outpath.SI, 'global.results.tex'))

# Regional morphometric development
str.morph.l.sl.p = str.morph.l.fu25 =str.morph.l.sl.t = str.morph.l.bl=str.morph.l.sl.sex.p = str.morph.l.sl.sex.t = err.idx = matrix(nrow = nmorph, ncol = nroi)
for (m in 1:nmorph) {
  print(paste0('Local development for phenotype: ',sm[m]))
  for (r in 1:nroi) {
    df = data.frame(age = MRI_table$age, sex = MRI_table$sex, id = MRI_table$id, morph = SM[r,m,], 
                    glob = MRI_table[,grep(paste0('glob.',sm[m]),colnames(MRI_table))])
    possibleError <- tryCatch({
      lm.morph.loc = lm(morph~age+sex,data=df, na.action = na.omit)
      lm.morph.loc.glob = lm(morph~age+sex+glob,data=df, na.action = na.omit)
    },
    error=function(e) {
      print(paste("Oops! --> Error in Feature ",m," Region ",r,sep = ""))
      err.idx[m,r] = TRUE
    }
    )
    if(inherits(possibleError, "error")) next
    str.morph.l.sl.p[m,r] = summary(lm.morph.loc)$coefficients[2,4]
    str.morph.l.sl.t[m,r] = summary(lm.morph.loc)$coefficients[2,3]
    str.morph.l.bl[m,r]=summary(lm.morph.loc)$coefficients[1,1]+14*summary(lm.morph.loc)$coefficients[1,2]+0.5*summary(lm.morph.loc)$coefficients[1,3]
  }
}
str.morph.l.sl.p.fdr=apply(str.morph.l.sl.p,1, function(i) p.adjust(i, method='fdr'))

write.csv(str.morph.l.sl.p, file=paste0(outpath.main,'str.morph.l.sl.p.csv'), row.names=F)
write.csv(str.morph.l.sl.t, file=paste0(outpath.main,'str.morph.l.sl.t.csv'), row.names=F)
write.csv(str.morph.l.bl, file=paste0(outpath.main,'str.morph.l.bl.csv'), row.names=F)

# SI Fig 16, Panel B
globeff = gather(MRI_table[,c('id','age','sex',paste0('glob.',sm))], measure, value, paste0('glob.',sm),factor_key = T)
globeff$measure=factor(globeff$measure,labels=sm.long)
globeff$measure=factor(globeff$measure,levels=sm.long[unlist(lapply(1:3,function(x) c(x,x+3)))])
df = globeff
df$measure = factor(df$measure, levels=sm.long)
newdat = expand.grid(age=range(MRI_table$age),
                     sex= unique(MRI_table$sex))
collect.glob.plots = list()
for (i in 1:nmorph) {
  m = unique(df$measure)[i]
  df.m = subset(df, measure ==m)
  lm.m = lm(value~age+sex,data=df.m)
  pred= predict.lm(lm.m, newdata=newdat, se.fit=T, level=0) 
  p = ggplot(df.m, aes(x = age, y = value, color = sex)) +
    geom_point(aes(alpha=0.5)) + 
    geom_ribbon(data=newdat[newdat$sex == 'Female',], aes(x=age, ymin=pred$fit[which(newdat$sex == 'Female')]-2*pred$se.fit[which(newdat$sex == 'Female')], ymax=pred$fit[which(newdat$sex == 'Female')]+2*pred$se.fit[which(newdat$sex == 'Female')]), alpha=0.5, inherit.aes=F) +
    geom_ribbon(data=newdat[newdat$sex == 'Male',], aes(x=age, ymin=pred$fit[which(newdat$sex == 'Male')]-2*pred$se.fit[which(newdat$sex == 'Male')], ymax=pred$fit[which(newdat$sex == 'Male')]+2*pred$se.fit[which(newdat$sex == 'Male')]), alpha=0.5, inherit.aes=F) +  
    geom_line(data=newdat, aes(y=pred$fit),size=1)+
    theme_minimal() +
    ylab('') +
    xlab('Age') +
    theme(text=element_text(size=14,  family="GillSans"),
          plot.title = element_text(hjust = 0.5), legend.position = 'bottom')+
    ggtitle(m)+ 
    scale_color_manual(values=c(colors.morph[7],colors.morph[1]))
  ggsave(plot = p, filename = paste0(outpath.main,gsub(' ','_',m),'.png'),width = 3,height=3, bg='white')
}

printout.stats[[4]] <- paste(paste0(paste0(' Significant ROIs in ',sm,': '),colSums(str.morph.l.sl.p.fdr<0.05),'; '), collapse='')

# Plot regional morphometric development
df.roi.morph.effects = data.frame(age = as.vector(t(str.morph.l.sl.t)), 
                                  age.p = as.vector(t(str.morph.l.sl.p)),
                                  age.fdr = as.vector(str.morph.l.sl.p.fdr),
                                  age.thresh = ifelse(as.vector(str.morph.l.sl.p.fdr)<0.05,as.vector(t(str.morph.l.sl.t)),NA),
                                  bl = as.vector(t(str.morph.l.bl)),
                                  mm = rep(factor(sm.long, levels = sm.long[unlist(lapply(1:3,function(x) c(x,x+3)))]), each=nroi), 
                                  label = rep(nm.ggseg,nmorph)) 
limvar=8
# SI Fig. 16, Panel C
panel.morphroi <- df.roi.morph.effects %>%
  group_by(mm)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=age))+
  scale_fill_gradientn(colors=colors.morph,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~mm, ncol=2)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

##################################
#    MORPHOMETRIC SIMILARITY
##################################
#transform features non-parametrically using MAD
SM.mad = array(NA,dim=dim(SM))
for (m in 1:nmorph) {
  print(m)
  for (s in 1:nsub) {
    SM.mad[,m,s] = (SM[,m,s]-median(SM[,m,s], na.rm = T))/mad(SM[,m,s], na.rm=T)
  }
}

# Calculate MSNs
MSN.mad = array(NA, dim=c(nroi,nroi,nsub)) # new "non-parametric" MSN = MAD + Pearson's r
for (s in 1:nsub) {
  MSN.mad[,,s] = cor(t(SM.mad[,,s]))
  MSN.mad[,,s][seq(1,nroi^2,by=nroi+1)] = NA
}
MSN = MSN.mad
MSN.group = apply(MSN,c(1,2),function(i) mean(i, na.rm=TRUE))
str = apply(MSN, 3, function(i) rowMeans(i,na.rm=T))

# ROI MSN Development
newdat <- expand.grid(age=c(min(MRI_table$age),
                            max(MRI_table$age)), 
                      sex = unique(MRI_table$sex))

str.l.sl.p = str.l.sl.t = str.l.sl = str.l.14 = str.l.sl.sex.p = str.l.sl.sex.t  = vector(length=nroi)
roi.ranef = matrix(nrow = nroi, ncol = length(unique(MRI_table$id)))
edf = vector(length = nroi)
p=list()
for (r in 1:nroi) {
  
  print(paste0('ROI: ',r))
  df.roi =  data.frame(id = MRI_table$id, age = MRI_table$age, sex = MRI_table$sex, str =str[r,])
  
  # Linear Models
  lm.roi = tryCatch({ 
    lm.roi <- lm(str~age+sex, data = df.roi,na.action = na.omit) }
    , error = function(e) {an.error.occured <<- TRUE})
  
  if(!(class(lm.roi)=='lm')){next}else{

    str.l.sl.p[r] = summary(lm.roi)$coefficients[2,4]     # p-value for effect of age
    str.l.sl.t[r] = summary(lm.roi)$coefficients[2,3]
    str.l.sl.sex.p[r] = summary(lm.roi)$coefficients[3,4] # p-value for effect of sex
    str.l.sl.sex.t[r] = summary(lm.roi)$coefficients[3,3]
    str.l.sl[r] = lm.roi$coefficients[2]
    str.l.14[r] = lm.roi$coefficients[1] + lm.roi$coefficients[2]*14 + 1/2 * lm.roi$coefficients[3]
    
  }
}
# Multiple comparisons correction
str.l.sl.p.fdr = p.adjust(str.l.sl.p, method='fdr')
str.l.sl.sex.p.fdr = p.adjust(str.l.sl.sex.p, method='fdr')

# Save results 
write.table(str.l.sl.t, (paste0(outpath.main,'/str.l.sl.t.csv')), row.names = F, col.names =F)
write.table(str.l.sl.p, (paste0(outpath.main,'/str.l.sl.p.csv')), row.names = F, col.names =F)

# Number of regions with significant rate of change in MSN degree
printout.stats[[6]] <- paste0('Number of regions with significant age effects on MSN: ', sum(str.l.sl.p.fdr<0.05))
printout.stats[[7]] <- paste0('Number of regions with significant sex effects on MSN: ', sum(str.l.sl.sex.p.fdr<0.05))


## Age effects on MSN
## SI Fig. 17, Panel A
limvar=5 
panel.msnage <- data.frame(label=rep(gsub('R_','rh_R_',gsub('L_','lh_L_',nm)),2),
                           value=c(str.l.sl.t,ifelse(str.l.sl.p.fdr<0.05,str.l.sl.t,NA)), 
                           group = factor(rep(c('non-corrected','FDR-corrected'),each = nroi)))  %>%
  group_by(group) %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=value))+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-5,5), oob=squish)+
  facet_wrap(~group, ncol=1)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  #labs(fill='Age (t-value)')+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# Rate of change by mesulam zone
df.mesulam.age = data.frame(age.t = str.l.sl.t,class = as.factor(mesulam$label2))
classmu.age = sapply(unique(mesulam$label2),function(class) mean(str.l.sl.t[mesulam$label2==class],na.rm=TRUE))
for (class in unique(mesulam$label2)) {
  df.mesulam.age$classmu[df.mesulam.age$class == class] = classmu.age[class]
}

## SI Fig. 17, Panel B
panel.roimesulam = ggplot(subset(df.mesulam.age, class != 'cortical_wall' & class != 'no label'), 
                          aes(x=class, y=age.t, fill = classmu)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white") +
   scale_fill_gradientn(colours=colors.msn.age, limits = c(-5, 5), oob=squish) +
  coord_flip()+
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.y = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3")),
        axis.title=element_text(size=12), legend.position = 'bottom')+
  ylab(expression(Delta~MS[14-26]))+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

write.csv(df.mesulam.age, file=(paste0(outpath.main,'mesulam-age-effect.csv')))
mean.comp.roimesulam = compare_means(age.t~class,subset(df.mesulam.age, class != 'cortical_wall' & class != 'no label'))
mean.comp.roimesulam=mean.comp.roimesulam[,-1]
print(xtable(mean.comp.roimesulam, type = "latex"), file = paste0(outpath.SI,'sign.mesulam.morph.tex'))

# Are these means significantly different from 0?
p.adjust(unlist(lapply(unique(df.mesulam.age$class)[1:4],function(x) t.test(subset(df.mesulam.age, class==x)$age.t,mu=0)$p.value)), method='fdr')<0.05

# Divergence
df=data.frame(#cor=cor(t(morph.dev.output$str.morph.l.sl.t[choose.SM[[w]],]),str.l.sl.t),
  cor=apply(t(str.morph.l.sl.t), 2, function(x) cor.test(x,str.l.sl.t)$estimate),
  p=apply(t(str.morph.l.sl.t), 2, function(x) cor.test(x,str.l.sl.t)$p.value),
  p.spin=apply(t(str.morph.l.sl.t), 2, 
               function(x) perm.sphere.p(x,str.l.sl.t,perm.id,corr.type = 'spearman')),
  sm=sm,group=sm.group)
df$p.fdr = p.adjust(df$p, method='fdr')
df$p.spin.fdr = p.adjust(df$p.spin, method='fdr')
df$sign = NA; df$sign[df$p.fdr<0.05] = 1; df$sign[df$p.spin.fdr<0.05] = 2;
df$sign= as.factor(df$sign)
df$sm=factor(df$sm, levels=df$sm[order(df$cor)])
df$group=factor(df$group, levels=c('Myelination','Grey Matter'),labels = c('microstructure','morphology'))

# SI Fig. 17, Panel C
divergence.panel=ggplot(df,aes(x=cor,y=sm, fill=group))+
  geom_col()+
  scale_fill_manual(values = c(colors.msn.age[1],colors.msn.age[7]))+
  scale_shape_manual(values = c(1, 16), labels=c(expression(P[FDR]*"<0.05"),expression(P[spin]*"<0.05")),name=) + 
  theme_minimal(base_family = "Gill Sans")+
  ylab('')+
  xlab('Correlation with MSN age effect')+
  labs(fill='')+
  theme(legend.position = 'bottom')




