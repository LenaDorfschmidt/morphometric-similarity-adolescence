#### Load required packages
library(dplyr)        # data wrangling
library(tidyr)        # data wrangling
library(kableExtra)   # data wrangling

library(nlme)         # linar mixed effects models
library(AICcmodavg)   # need predictSE.lme for Lme confidence interval

library(ggplot2)      # plotting
library(ggseg)        # brain surface plotting
library(ggsegGlasser) # brain surface plotting
library(paletteer)    # color scales
library(ggpubr)       # arrange plots
library(fields)       # image.plot function 
library(wordcloud)    # for Neurosynth result wordcloud
library(extrafont)    # to use Gill Sans font, first time use: library(extrafont); font_import(); loadfonts()

library(matrixStats)  # matrix manipulation
library(Hmisc)        # compute correlation on matrix rows 
library(scales)       # squish function
library(igraph)       # network abalysis

library(xtable)       # saving tables for latex doc

# Set working directory
setwd("~/Documents/Cambridge/NSPN_MSN/github-code/morphometric-similarity-adolescence/")

#### Load data
load((paste0('data/NSPN.RData')))

#### Load dependencies
source('scripts/dependencies/glob.outliers.R')                            # Global outlier detection
source('scripts/dependencies/local.outliers.R')                           # Local outlier detection
source('scripts/external/rotate_parcellation-master/R/perm.sphere.p.R')   # Spin-test p-value estimation

##################################
#         ANALYSIS SETUP
##################################
# Create colorscales for later plotting use
colors.morph = paletteer_d('rcartocolor::Tropic') 
colors.msn.age = paletteer_d('Redmonder::dPBIPuOr')
colors.msn.str = paletteer_d('Redmonder::dPBIRdGy')
colors.cor = paletteer_d('Redmonder::dPBIRdBu')

# Setup output directories 
outpath='results/'
outpath.SI = paste0(outpath,'/SI/')
outpath.main = paste0(outpath,'/main/')
latex.tables = paste0(outpath,'/latex-tables/')
dir.create(outpath, recursive = T)
dir.create(outpath.main)
dir.create(outpath.SI)
dir.create(latex.tables)

# To make our live easier we will print lots of info into this
printout.stats=list() 

#################################################
#            GLOBAL OUTLIER DETECTION
#################################################
# Estimate global outliers in morphometric features
# Creates SI Fig 2
glob.outl.output = global.outliers(SM,MRI_table,sm,sm.group, paste0(outpath.SI,'outliers/'))

# Finds outliers subjects for exlcusion
excl.sub.mad = glob.outl.output$excl.sub.mad
MRI_table=cbind(MRI_table, glob.outl.output$glob.morph)

printout.stats[[1]] <- paste0('Number of subjects excluded due to global outliers: ', length(excl.sub.mad))
  
# Delete subjects based on global outliers
SM = SM[,,-excl.sub.mad]              # morphometric features
FC = FC[,,-excl.sub.mad]              # functional connectivity matrices
MRI_table = MRI_table[-excl.sub.mad,] # demographics
nsub = dim(MRI_table)[1]

#################################################
#           LOCAL OUTLIER DETECTION
#################################################
# Estimate global outliers in morphometric features
# Creates SI Fig 3
loc.outl.output = local.outliers(SM,MRI_table, sm,sm.long, paste0(outpath.SI,'outliers/'), nm.ggseg)
loc.outl = loc.outl.output$loc.outl; del.roi = loc.outl.output$del.roi

printout.stats[[2]] <- paste0('Number of regions excluded due to local outliers: ', length(del.roi))
# We will save region names for replication in HCP-D (since we are excluding regions here, that might not be exlcuded in HCP-D)
write.csv(nm[-del.roi], file=paste0(outpath.main,'NSPN.nm.txt'))

# Delete ROIs based on local outliers
if(any(del.roi)){
  SM = SM[-del.roi,,] 
  nroi = dim(SM)[1]
  nm.ggseg = nm.ggseg[-del.roi]
  nmorph = length(sm)
  loc.outl = loc.outl[,-del.roi,]
  Yeo=Yeo[-del.roi]
  mesulam=mesulam[-del.roi]
  vonEconomo=vonEconomo[-del.roi]
}

# Set all outl>=5 to NA
loc.outl = aperm(loc.outl, c(2,1,3))
SM = replace(SM, loc.outl >= 5, NA)

printout.stats[[3]] <- paste0('Outliers percentage (number of ROIs set to NA across all subjects, regions, phenotypes): ',sum(loc.outl >= 5, na.rm=T)/length(loc.outl)*100,'%')

#################################################
#            FINAL DATASET
#################################################
# Based on outliers above, we exclude curvature features
choose.features = which(is.na(match(sm,c('MC','CI','FI','IC'))))

SM=SM[,choose.features,]
nmorph = dim(SM)[2]
sm=sm[choose.features] # c("FA" "MD" "MT" "SA" "GM" "CT")
sm.long=sm.long[choose.features]
sm.group = sm.group[choose.features]
SM.param=aperm(apply(SM,c(3,2),scale),c(1,3,2)) # z-score and reorder dimensions

# Data Overview - which structural scans have FC data, too
MRI_table$session = factor(MRI_table$session)
MRI_table$fc.exists=!(apply(FC,3, function(x) sum(is.na(x)))==129600)
demographics.table = MRI_table %>% group_by(session,sex) %>% dplyr::summarize(N=n(), hasFC= sum((fc.exists)))
print(xtable(demographics.table, type = "latex"), file = paste0(latex.tables,'has.fc.tex'))

# I had to post-edit the following table because I just couldn't figure out how to output it prettyly, ugh
table.scans<- ftable(xtabs( ~ sex + session + agebin, MRI_table ),col.vars=c(2,1))
hasfc<- ftable(xtabs(fc.exists ~ sex + agebin, MRI_table ), col.vars=c(1))
hasfc<- ftable(xtabs(fc.exists ~ sex + session + agebin, MRI_table ),col.vars=c(2,1))
print(xtableFtable(table.scans, type = "latex",method = 'compact',), file = paste0(latex.tables,'scans.per.session.tex'))

# Acquisition overview
# SI Fig 1
acq.overview=MRI_table %>% group_by(id) %>% mutate(firstscan=order(age), firstage=min(age))
relevel=acq.overview[acq.overview$firstscan==1,]$id[order(acq.overview[acq.overview$firstscan==1,]$age)]
acq.overview$id = factor(acq.overview$id, levels = relevel)
acq.overview$y = as.numeric(acq.overview$id)
plot.acq.overview = ggplot(acq.overview, aes(y=y,x=age))+
  geom_line(aes(group=id,color=firstage))+
  geom_point()+
  labs(color='Baseline Age')+
  theme_bw()+ylab('Participant')+
  scale_color_gradientn(colours = paletteer_c('viridis::magma',100))
ggsave(plot=plot.acq.overview, filename=paste0(outpath.SI, 'acquisition.overview.png'),height=4,width=4)


#################################################
#              PERMUTATION TESTING
################################################# 
# For permutation analyses, either uncomment and run the following lines (takes a while, not recommended)
# set.seed(1234)
# source('scripts/external/rotate_parcellation-master/R/rotate.parcellation.R')
# # read-in spherical coordinates
# hcp.centroids = as.matrix(read.table('scripts/external/rotate_parcellation-master/R/coords_MMP_sphere_geodesic.txt')) # available from https://github.com/frantisekvasa/rotate_parcellation
# # obtain coordinates of regions on L and R hemispheres
# coord.l = hcp.centroids[1:180,];              # coordinates of L hemisphere
# coord.r = hcp.centroids[181:360,];            # coordinates of R hemisphere
# # delete dropout regions
# coord.l = coord.l[-del.roi[del.roi<=180],];
# coord.r = coord.r[-(del.roi[del.roi>180]-180),];
# # run spherical permutation script
# perm.id = rotate.parcellation(coord.l,coord.r,nrot=10000)
# save(perm.id, file = 'scripts/replication/data/perm.id.RData')
# 
# # delete dropout regions
# coord.l = hcp.centroids[fc.keep.id[fc.keep.id<=180],];
# coord.r = hcp.centroids[fc.keep.id[fc.keep.id>180],];
# # run spherical permutation script
# perm.id.fc = rotate.parcellation(coord.l,coord.r,nrot=10000)
# save(perm.id.fc, file = 'scripts/replication/data/perm.id.fc.RData')

# Alternatively, import the results
load('data/perm.id.RData')

#################################################
#       GLOBAL MORPHOMETRIC DEVELOPMENT
################################################# 
# Global plots
newdat <- expand.grid(age=c(min(MRI_table$age),
                            max(MRI_table$age)),
                      sex = unique(MRI_table$sex),
                      site = factor(which.max(table(MRI_table$site))))
glob.res=matrix(NA, nrow=nmorph, ncol=4); rownames(glob.res)=sm.long; colnames(glob.res)=c('t.age','p.age','t.sex','p.sex')
glob.res = data.frame(glob.res)
for (m in 1:nmorph) {
  sm.str=paste0('glob.',sm[m])
  df=MRI_table[,c(c('id','age','sex','site',sm.str))]
  colnames(df)[grep(sm.str,colnames(df))] = 'cur'
  lm.glob = lme(cur~age+sex+site+age*sex, random = ~1|id,data = df)
  glob.res[m,1:2]=summary(lm.glob)$tTable[2,4:5];
  glob.res[m,3:4]=summary(lm.glob)$tTable[3,4:5]; 
  #lm.pred = predictSE.lme(lm.glob, level=0, newdata=newdat)
}
glob.res$p.age.fdr = p.adjust(glob.res$p.age, method='fdr')
glob.res$p.sex.fdr = p.adjust(glob.res$p.sex, method='fdr')
print(xtable(glob.res, type = "latex"), file=paste0(latex.tables, 'global.results.tex'))

# Global morphometric development
# Fig 2, Panel A
globeff = gather(MRI_table[,c('id','age','sex','site',paste0('glob.',sm))], measure, value, paste0('glob.',sm),factor_key = T)
globeff$measure=factor(globeff$measure,labels=sm.long)
globeff$measure=factor(globeff$measure,levels=sm.long[unlist(lapply(1:3,function(x) c(x,x+3)))])
df = globeff
df$site = MRI_table$site
df$measure = factor(df$measure, levels=sm.long)
newdat = expand.grid(age=range(MRI_table$age),
                     sex= unique(MRI_table$sex),
                     site=as.factor(3))
for (i in 1:nmorph) {
  m = unique(df$measure)[i]
  df.m = subset(df, measure ==m)
  lm.m = lme(value~age+sex+site, random=~1|id,data=df.m)
  pred= predictSE.lme(lm.m, newdata=newdat, se.fit=T, level=0) 
  p = ggplot(df.m, aes(x = age, y = value, color = sex)) +
    geom_point(aes(alpha=0.5)) + 
    geom_line(aes(group=id,alpha=0.5)) +
    geom_ribbon(data=newdat[newdat$sex == 0,], aes(x=age, ymin=pred$fit[which(newdat$sex == 0)]-2*pred$se.fit[which(newdat$sex == 0)], ymax=pred$fit[which(newdat$sex == 0)]+2*pred$se.fit[which(newdat$sex == 0)]), alpha=0.5, inherit.aes=F) +
    geom_ribbon(data=newdat[newdat$sex == 1,], aes(x=age, ymin=pred$fit[which(newdat$sex == 1)]-2*pred$se.fit[which(newdat$sex == 1)], ymax=pred$fit[which(newdat$sex == 1)]+2*pred$se.fit[which(newdat$sex == 1)]), alpha=0.5, inherit.aes=F) +  
    geom_line(data=newdat, aes(y=pred$fit),size=1)+
    theme_minimal() +
    ylab('') +
    xlab('Age') +
    theme(text=element_text(size=14,  family="GillSans"),
          plot.title = element_text(hjust = 0.5), legend.position = 'bottom')+
    ggtitle(m)+ 
    scale_color_manual(values=c(colors.morph[7],colors.morph[1]))
  ggsave(plot = p, filename = paste0(outpath.main,'glob.morph.development.',sm[i],'.png'),width = 3,height=3, bg='white')
}

global.effect.sizes=lapply(sm.long, function(x) summary(lme(value~age+sex+site, random= ~1|id,data = subset(globeff,measure==x)))$tTable[2,4:5])
printout.stats[[5]] = paste0(paste(rep(sm, each=2),c('t-value: ','p-value: ')),unlist(global.effect.sizes),';', collapse='')

#################################################
#       REGIONAL MORPHOMETRIC DEVELOPMENT
################################################# 
# Estimate regional morphometric feature development
# Creates SI Fig. 5 and  SI Fig. 6
# Regional morphometric development
str.morph.l.sl.p = str.morph.l.fu25 =str.morph.l.sl.t = str.morph.l.bl=str.morph.l.sl.sex.p = str.morph.l.sl.sex.t = err.idx = matrix(nrow = nmorph, ncol = nroi)
str.morph.l.sl.p.glob.corr = str.morph.l.sl.t.glob.corr = err.idx.glob.corr =  matrix(nrow = nmorph, ncol = nroi)
for (m in 1:nmorph) {
  print(paste0('Local development for phenotype: ',sm[m]))
  for (r in 1:nroi) {
    # Regional development, not correct for global feature value
    df = data.frame(age = MRI_table$age, sex = MRI_table$sex, center = MRI_table$site, 
                    id = MRI_table$id, morph = SM[r,m,], 
                    glob = MRI_table[,grep(paste0('glob.',sm[m]),colnames(MRI_table))])
    
    possibleError <- tryCatch({
      lm.morph.loc = lme(morph~age+sex+center,random=~1|id,data=df, na.action = na.omit)
    },
    error=function(e) {
      print(paste("Oops! --> Error in Feature ",m," Region ",r,sep = ""))
      err.idx[m,r] = TRUE
    }
    )
    if(inherits(possibleError, "error")) next
    str.morph.l.sl.p[m,r] = summary(lm.morph.loc)$tTable[2,5]
    str.morph.l.sl.t[m,r] = summary(lm.morph.loc)$tTable[2,4]
    str.morph.l.bl[m,r]=summary(lm.morph.loc)$coefficients$fixed[1]+14*summary(lm.morph.loc)$coefficients$fixed[1]+0.5*summary(lm.morph.loc)$coefficients$fixed[3]

    # Correct for global feature value
    df = data.frame(age = MRI_table$age, sex = MRI_table$sex, center = MRI_table$site, id = MRI_table$id, 
                    morph = SM[r,m,], morph.glob =  MRI_table[,grep(paste0('glob.',sm[m]),colnames(MRI_table))])
    
    possibleError <- tryCatch({
      lm.morph.loc = lme(morph~age+sex+center+morph.glob,random=~1|id,data=df, na.action = na.omit)
    },
    error=function(e) {
      print(paste("Oops! --> Error in Feature ",m," Region ",r,sep = ""))
      err.idx.glob.corr[m,r] = TRUE
    }
    )
    if(inherits(possibleError, "error")) next
    str.morph.l.sl.p.glob.corr[m,r] = summary(lm.morph.loc)$tTable[2,5]
    str.morph.l.sl.t.glob.corr[m,r] = summary(lm.morph.loc)$tTable[2,4]
    }
}
str.morph.l.sl.p.fdr=apply(str.morph.l.sl.p,1, function(i) p.adjust(i, method='fdr'))
str.morph.l.sl.p.fdr.glob.corr=apply(str.morph.l.sl.p.glob.corr,1, function(i) p.adjust(i, method='fdr'))
printout.stats[[4]] <- paste(paste0(paste0(' Significant ROIs in ',sm,': '),colSums(str.morph.l.sl.p.fdr<0.05),'; '), collapse='')

write.csv(str.morph.l.sl.p, file=paste0(outpath.main,'str.morph.l.sl.p.csv'), row.names=F)
write.csv(str.morph.l.sl.t, file=paste0(outpath.main,'str.morph.l.sl.t.csv'), row.names=F)
write.csv(str.morph.l.bl, file=paste0(outpath.main,'str.morph.l.bl.csv'), row.names=F)

# Plot regional morphometric development
df.roi.morph.effects = data.frame(age = as.vector(t(str.morph.l.sl.t)), 
                                  age.p = as.vector(t(str.morph.l.sl.p)),
                                  age.fdr = as.vector(str.morph.l.sl.p.fdr),
                                  age.thresh = ifelse(as.vector(str.morph.l.sl.p.fdr)<0.05,as.vector(t(str.morph.l.sl.t)),NA),
                                  mm = rep(factor(sm.long, levels = sm.long[unlist(lapply(1:3,function(x) c(x,x+3)))]), each=nroi), 
                                  label = rep(nm.ggseg,nmorph)) 
limvar=8
panel.morphroi <- df.roi.morph.effects %>%
  group_by(mm)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=age),hemisphere = ifelse(average.hemi,'left',c('left','right')))+
  scale_fill_gradientn(colors=colors.morph,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~mm, ncol=2)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

# Fig 2, Panel B
ggsave(plot=panel.morphroi, filename=paste0(outpath.main, 'regional.morphometric.development.png'),height=5,width=5,bg='white')

# Regional morphometric development, thresholded for significance
# SI Fig. 5
SI.age.thresh <- df.roi.morph.effects %>%
  group_by(mm)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=age.thresh),hemisphere = ifelse(average.hemi,'left',c('left','right')))+
  scale_fill_gradientn(colors=colors.morph,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~mm, ncol=2)+
  #theme_void()+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

ggsave(filename=paste0(outpath.SI,'msn.age.regional.png'), plot=SI.age.thresh,width = 8,height=4, bg='white')

# Regional morphometric development, corrected for global trend
df.roi.morph.effects.corr = data.frame(age = as.vector(t(str.morph.l.sl.t.glob.corr)), 
                                  age.p = as.vector(t(str.morph.l.sl.p.glob.corr)),
                                  age.fdr = as.vector(str.morph.l.sl.p.fdr.glob.corr),
                                  age.thresh = ifelse(as.vector(str.morph.l.sl.p.fdr.glob.corr)<0.05,as.vector(t(str.morph.l.sl.t.glob.corr)),NA),
                                  mm = rep(factor(sm.long, levels = sm.long[unlist(lapply(1:3,function(x) c(x,x+3)))]), each=nroi), 
                                  label = rep(nm.ggseg,nmorph)) 
# SI Fig. 6, Panel A
SI.age.corr = df.roi.morph.effects.corr %>%
  group_by(mm)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=age),hemisphere = ifelse(average.hemi,'left',c('left','right')))+
  scale_fill_gradientn(colors=colors.morph,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~mm, ncol=2)+
  #theme_void()+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom', text = element_text(size = 12))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

# SI Fig. 6, Panel B
SI.age.thresh.corr = df.roi.morph.effects.corr %>%
  group_by(mm)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=age.thresh),hemisphere = ifelse(average.hemi,'left',c('left','right')))+
  scale_fill_gradientn(colors=colors.morph,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~mm, ncol=2)+
  #theme_void()+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom', text = element_text(size = 12))+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

# SI Fig. 6
ggsave(plot = ggarrange(SI.age.corr,SI.age.thresh.corr, 
                        nrow=1, vjust=0.7,hjust = 0,
                        labels = c('A | Regional morphometric development',
                                   'B | Regional development thresholded'),
                        font.label = list(size = 13, family="Gill Sans"))+
         theme(plot.margin = margin(0.5,0,0,0, "cm")),
       filename = paste0(outpath.SI,'morph.age.global.corrtion.png'), width=8,height=5, bg='white')

#################################################
#       WITHIN-SUBJECT REGIONAL CHANGE
################################################# 
# Within-subject individual feature effects
BLFUDIFF=MRI_table %>%
  filter(session %in% c('ses-baseline', 'ses-1stFollowUp')) %>%
  group_by(id) %>%
  mutate(age_baseline = ifelse(session == "ses-baseline", age, NA)) %>%
  mutate(site_baseline = ifelse(session == "ses-baseline", site, NA)) %>%
  summarise(across(starts_with("glob."), ~ diff(.)),
            sex = first(sex),
            age_baseline = first(age_baseline),
            site_baseline = first(site_baseline))

BLFUDIFF.LONG = BLFUDIFF %>%
  pivot_longer(cols = starts_with("glob."),
               names_to = "variable",
               values_to = "value")

BLFUDIFF.LONG = BLFUDIFF.LONG[which(!is.na(match(BLFUDIFF.LONG$variable,
                                                 c('glob.FA','glob.MD','glob.MT','glob.GM','glob.CT','glob.SA')))),] 
BLFUDIFF.LONG$variable = factor(BLFUDIFF.LONG$variable, 
                                levels=c('glob.CT','glob.GM','glob.SA',
                                         'glob.FA','glob.MT','glob.MD'),
                                labels=c('CT','GM','SA','FA','MT','MD'))

PLOT.DIFF.MORPH = list() 
for (m in levels(BLFUDIFF.LONG$variable)){
  SUBSET=subset(BLFUDIFF.LONG,variable==m)
  lim.max = max(abs(SUBSET$value))
  mean_value <- mean(SUBSET$value)
  se_value <- sd(SUBSET$value) / sqrt(nrow(SUBSET))
  CI.low = mean_value + 1.96*se_value
  CI.high = mean_value - 1.96*se_value
  
  PLOT.DIFF.MORPH[[m]] = SUBSET %>%
    ggplot(aes(y = value, x = variable)) +
    geom_jitter(aes(color=value))+
    geom_boxplot(aes(fill=mean(value))) +
    theme_minimal(base_family = 'Gill Sans') +
    ylab('') +
    scale_y_continuous(limits = function(x) c(-max(abs(x)), max(abs(x))))+
    xlab('')+
    scale_color_gradientn(colours=colors.morph,limits = c(-lim.max,lim.max))+
    scale_fill_gradientn(colours=colors.morph,limits = c(-lim.max,lim.max))
}

# SI Fig. 4
plot.within.sub = ggarrange(plotlist = PLOT.DIFF.MORPH, legend = F,nrow=1)
ggsave(plot=plot.within.sub, filename=paste0(outpath.SI, 'within.subject.change.png'),height=2.5,width=5)

# percentage of subjects with specific within-subject trend
BLFUDIFF.LONG %>% group_by(variable) %>% summarise(Nneg=sum(value<0)/n()*100) 

#################################################
#        REGIONAL DEVLOPMENT BY MESULAM ZONES
################################################# 

mm.effect = as.data.frame(t(str.morph.l.sl.t)); colnames(mm.effect) = sm
mm.effect$ve = vonEconomo
mm.ve.mean = mm.effect%>% group_by(ve) %>%  summarise(across(everything(), mean))
mm.ve.mean$name = ve.name 

df.regional.mesulam = data.frame(t=as.vector(t(str.morph.l.sl.t)),feature=rep(sm,each=nroi),class=rep(mesulam,6))

mesulam.morph = subset(df.regional.mesulam, class != 'cortical_wall' & class != 'no label') %>% 
  group_by(feature, class) %>% mutate(mean_t=mean(t)) %>% ggplot(aes(x=class, y=t,fill=mean_t)) + 
  facet_grid(~factor(feature, levels=sm[order(sm.group)]))+
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white")+scale_fill_gradientn(colours=colors.morph, limits = c(-5, 5), oob=squish) +
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.x = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3"),angle = 90),
        axis.title=element_text(size=12),text = element_text(size=12), legend.position = 'bottom')+
  ylab(expression(Delta~MS[14-26]))+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig. 7
ggsave(plot=mesulam.morph, filename=paste0(outpath.SI, 'regional.morphometric.development.mesulam.png'),height=5,width=5,bg='white')

# Between-mesulam zone morphometric development
comp.mean.morph.mesulam = lapply(unique(subset(df.regional.mesulam, class != 'cortical_wall' & class != 'no label')$feature),
                                 function(i) 
                                   compare_means(t~class,
                                                 subset(subset(df.regional.mesulam, 
                                                               class != 'cortical_wall' & 
                                                                 class != 'no label' & 
                                                                 feature==i))))

comp.mean.morph.mesulam = lapply(comp.mean.morph.mesulam,function(i) i[,c(2,3,4,5,7)])
comp.mean.morph.mesulam = data.table::rbindlist(comp.mean.morph.mesulam)
comp.mean.morph.mesulam$feature = rep(sm,each=6)
comp.mean.morph.mesulam= comp.mean.morph.mesulam[,c(6,1:5)]

out.latex.tab.morph.mesulam = comp.mean.morph.mesulam %>% kbl(
  format = "latex",
  row.names = F,
  booktabs = T,
  escape = F,
  linesep = c("", "", "","","", "\\midrule")
) %>% collapse_rows(columns = 1, latex_hline = "none")

print(out.latex.tab.morph.mesulam, file = paste0(latex.tables,'tab.morph.mesulam.tex'))

#################################################
#            MORPHOMETRIC SIMILARITY
################################################# 
#transform features non-parametrically using MAD
SM.mad = array(NA,dim=dim(SM))
for (m in 1:nmorph) {
  print(m)
  for (s in 1:nsub) {
    SM.mad[,m,s] = (SM[,m,s]-median(SM[,m,s], na.rm = T))/mad(SM[,m,s], na.rm=T)
  }
}

# Calculate MSNs
MSN.param = array(NA, dim=c(nroi,nroi,nsub)) # standard "parametric" MSN = Z-score + Pearson's r
MSN.mad = array(NA, dim=c(nroi,nroi,nsub)) # new "non-parametric" MSN = MAD + Pearson's r
MSN.mad.mad = array(NA, dim=c(nroi,nroi,nsub)) # new "non-parametric" MSN = MAD + Spearman's rho
for (s in 1:nsub) {
  MSN.param[,,s] = cor(t(SM.param[,,s]))
  MSN.param[,,s][seq(1,nroi^2,by=nroi+1)] = NA
  MSN.mad[,,s] = cor(t(SM.mad[,,s]))
  MSN.mad[,,s][seq(1,nroi^2,by=nroi+1)] = NA
  MSN.mad.mad[,,s] = cor(t(SM.mad[,,s]), method='spearman')
  MSN.mad.mad[,,s][seq(1,nroi^2,by=nroi+1)] = NA
}
MSN = MSN.mad
MSN.group = apply(MSN,c(1,2),function(i) mean(i, na.rm=TRUE))
str = apply(MSN, 3, function(i) rowMeans(i,na.rm=T))

# Regional MSN Development
newdat <- expand.grid(age=c(min(MRI_table$age),
                            max(MRI_table$age)), 
                      sex = unique(MRI_table$sex),
                      center = factor(names(which.max(table(MRI_table$mri_center)))))

str.l.sl.beta = str.l.sl.p = str.l.sl.t = str.l.sl = str.l.14 = str.l.site.p= str.l.site.F= str.l.sl.sex.p = str.l.sl.sex.t = str.l.sl.interaction.p= str.l.sl.interaction.t = p.int = t.int = aic = vector(length=nroi)
roi.ranef = matrix(nrow = nroi, ncol = length(unique(MRI_table$id)))
edf = vector(length = nroi)
p=list()
for (r in 1:nroi) {
  print(paste0('ROI: ',r))
  df.roi =  data.frame(id = MRI_table$id, center= MRI_table$site, age = MRI_table$age, sex = MRI_table$sex, str =str[r,])
  # Linear Models
  lm.roi = tryCatch({ 
    lm.roi <- lme(str~age+sex+center, random = ~1|id,data = df.roi,na.action = na.omit) }
    , error = function(e) {an.error.occured <<- TRUE})
  
  if(!(class(lm.roi)=='lme')){next}else{
    roi.ranef[r,match(rownames(lm.roi$coefficients$random$id),unique(MRI_table$id))] = as.vector(lm.roi$coefficients$random$id)
    str.l.sl.p[r] = summary(lm.roi)$tTable[2,5]     # p-value for effect of age
    str.l.sl.t[r] = summary(lm.roi)$tTable[2,4]
    str.l.sl.sex.p[r] = summary(lm.roi)$tTable[3,5] # p-value for effect of sex
    str.l.sl.sex.t[r] = summary(lm.roi)$tTable[3,4]
    str.l.sl[r] = lm.roi$coefficients$fixed[2]
    str.l.14[r] = lm.roi$coefficients$fixed[1] + lm.roi$coefficients$fixed[2]*14 + 1/2 * lm.roi$coefficients$fixed[3] + lm.roi$coefficients$fixed[4]*(1/3) + lm.roi$coefficients$fixed[5]*(1/3)
    str.l.site.p[r]=anova(lm.roi)[4,4]
    str.l.site.F[r]=anova(lm.roi)[4,3]
    str.l.sl.beta[r] = lm.roi$coefficients$fixed[2]
  }
}
str.l.sl.p.fdr = p.adjust(str.l.sl.p, method='fdr')         # multiple comparisons correction
str.l.sl.sex.p.fdr = p.adjust(str.l.sl.sex.p, method='fdr') # multiple comparisons correction
write.table(str.l.sl.t, paste0(outpath.main,'/str.l.sl.t.csv'), row.names = F, col.names =F)
write.table(str.l.sl.p, paste0(outpath.main,'/str.l.sl.p.csv'), row.names = F, col.names =F)

# Output for further analyses
save(MRI_table,MSN,str.l.sl.t,str.l.sl.beta,str.l.sl.p.fdr,str.l.14,mesulam,mesulam.col,nm,nm.ggseg,nsub,perm.id,
     colors.msn.age,outpath.main,outpath.SI,FC,fc.keep.id,del.roi,Yeo,nroi,file = paste0('data/for.further.analyses.RData'))

# Table of signifiant regions, with Mesulam zone assignment
df.sign.msn.age = data.frame(regionrty=gsub('_',' ',nm)[str.l.sl.p.fdr<0.05],tvalue=str.l.sl.t[str.l.sl.p.fdr<0.05],pval.fdr=str.l.sl.p.fdr[str.l.sl.p.fdr<0.05],network=mesulam[str.l.sl.p.fdr<0.05])
df.sign.msn.age=df.sign.msn.age[order(df.sign.msn.age$network),]
print(xtable(df.sign.msn.age, type = "latex"), file = paste0(outpath.SI,'sign.msn.age.mesulam.tex'),include.rownames=FALSE)

printout.stats[[6]] <- paste0('Number of regions with significant age effects on MSN: ', sum(str.l.sl.p.fdr<0.05))
printout.stats[[7]] <- paste0('Number of regions with significant sex effects on MSN: ', sum(str.l.sl.sex.p.fdr<0.05))

## Plot age effects on MSN
limvar=5
panel.msnage <- data.frame(label=rep(nm.ggseg,2),
                           value=c(str.l.sl.t,ifelse(str.l.sl.p.fdr<0.05,str.l.sl.t,NA)), 
                           group = factor(rep(c('non-corrected','FDR-corrected'),each = nroi)))  %>%
  group_by(group) %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=value))+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~group, ncol=1)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# Fig. 3, Panel A
ggsave(plot=panel.msnage, filename=paste0(outpath.main, 'msn.age.effect.png'),height=3,width=5, bg='white')

# MSN age effect by mesulam zone
df.mesulam.age = data.frame(age.t = str.l.sl.t,class = as.factor(mesulam))
classmu.age = sapply(unique(mesulam),function(class) mean(str.l.sl.t[mesulam==class],na.rm=TRUE))
for (class in unique(mesulam)) {
  df.mesulam.age$classmu[df.mesulam.age$class == class] = classmu.age[class]
}

write.csv(df.mesulam.age, file=paste0(outpath.main,'mesulam-age-effect.csv'))

mean.comp.roimesulam = compare_means(age.t~class,subset(df.mesulam.age, class != 'cortical_wall' & class != 'no label'))
mean.comp.roimesulam=mean.comp.roimesulam[,-1]
print(xtable(mean.comp.roimesulam, type = "latex"), file = paste0(latex.tables,'sign.mesulam.morph.tex'))

panel.roimesulam = ggplot(subset(df.mesulam.age, class != 'cortical_wall' & class != 'no label'), aes(x=class, y=age.t, fill = classmu)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white") +
  scale_fill_gradientn(colours=colors.msn.age, limits = c(-3, 3), oob=squish) +
  coord_flip()+
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.y = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3")),
        axis.title=element_text(size=12), legend.position = 'bottom')+
  ylab(expression(Delta~MS[14-26]))+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# Fig. 3, Panel B
ggsave(plot=panel.msnage, filename=paste0(outpath.main, 'msn.age.effect.png'),height=3,width=5, bg='white')

# Edgewise MSN Development
mat.str.l.sl.p = mat.str.l.sl.t = mat.str.l.sl = mat.str.l.14 = mat.str.l.sl.sex.p = mat.str.l.sl.sex.t = matrix(NA,nrow=nroi,ncol=nroi)
p=list()
for (i in 1:(nroi-1)) {
  print(paste0('ROI: ',i))
  for (j in (i+1):nroi) {
    df.edge =  data.frame(id = MRI_table$id, center= MRI_table$site, age = MRI_table$age, sex = MRI_table$sex, str =MSN[i,j,])
    # Linear Models
    lm.roi = tryCatch({
      lm.roi <- lme(str~age+sex+center, random = ~1|id,data = df.edge,na.action = na.omit) }
      , error = function(e) {an.error.occured <<- TRUE})
    if(!(class(lm.roi)=='lme')){next}else{
      mat.str.l.sl.p[i,j] = summary(lm.roi)$tTable[2,5]     # p-value for effect of age
      mat.str.l.sl.t[i,j] = summary(lm.roi)$tTable[2,4]
      mat.str.l.sl.sex.p[i,j] = summary(lm.roi)$tTable[3,5] # p-value for effect of sex
      mat.str.l.sl.sex.t[i,j] = summary(lm.roi)$tTable[3,4]
      mat.str.l.sl[i,j] = lm.roi$coefficients$fixed[2]
      mat.str.l.14[i,j] = lm.roi$coefficients$fixed[1] + lm.roi$coefficients$fixed[2]*14 + 1/2 * lm.roi$coefficients$fixed[3] + lm.roi$coefficients$fixed[4]*(1/3) + lm.roi$coefficients$fixed[5]*(1/3)
    }
  }
}

# Only estimated upper triangle above, create full matrix
T.MAT = mat.str.l.sl.t
T.MAT[which(is.na(T.MAT))] = 0
T.MAT = T.MAT+t(T.MAT)
P.MAT = mat.str.l.sl.p
P.MAT[which(is.na(P.MAT))] = 0
P.MAT = P.MAT+t(P.MAT)

assignment = mesulam
assignment[which(is.na(assignment))]='z'

# Fig. 3, Panel C
png(paste0(outpath.main,'msn.age.effect.matrix.png'))
plot.smooth.mat(T.MAT,P.MAT,as.numeric(as.factor(assignment)),
                as.numeric(as.factor(assignment)),
                colors.msn.age,mesulam.col,mesulam.col,col.nm=NULL)
dev.off()


# Divergence in morphometric similarity development
df=data.frame(cor=apply(t(str.morph.l.sl.t), 2, function(x) cor.test(x,str.l.sl.t)$estimate),
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
divergence.panel=ggplot(df,aes(x=cor,y=sm, fill=group))+
  geom_col()+
  geom_point(aes(shape=sign),size=4)+
  scale_fill_manual(values = c(colors.msn.age[1],colors.msn.age[7]))+
  scale_shape_manual(values = c(1, 16), labels=c(expression(P[FDR]*"<0.05"),expression(P[spin]*"<0.05")),name=) + 
  theme_minimal(base_family = "Gill Sans")+
  ylab('')+
  xlab('Correlation with MSN age effect')+
  labs(fill='')+
  theme(legend.position = 'bottom')

# Fig. 3, Panel D
ggsave(plot=panel.msnage, filename=paste0(outpath.main, 'msn.age.effect.png'),height=3,width=5, bg='white')


plot.mm.msn = data.frame(mm=as.vector(t(str.morph.l.sl.t)), 
                         msn=rep(str.l.sl.t,nmorph),
                         group = factor(rep(sm,each=nroi)),
                         pfdr = rep(str.l.sl.p.fdr,nmorph))
plot.mm.msn$group = factor(plot.mm.msn$group,levels=sm[order(sm.group)])

corvals.mm.msn= plot.mm.msn %>% 
  ggplot(aes(x=msn, y=mm)) +
  facet_wrap(~group, ncol=3) +
  geom_point(data = subset(plot.mm.msn, pfdr >= 0.05), 
             aes(fill = msn, color = msn), 
             alpha=0.8, size = 2, shape = 21, stroke = 1) +  # Fill points with pfdr < 0.05 according to msn, outline in black
  geom_point(data = subset(plot.mm.msn, pfdr < 0.05), 
             aes(fill=msn), 
             alpha=0.8, shape = 21, stroke = 1,color = 'grey30',) +  # Plot other points normally
  geom_hline(yintercept = 0, linetype='dotted') +
  scale_color_gradientn(colors=colors.msn.age, limits = c(-5,5), oob=squish) +
  scale_fill_gradientn(colors=colors.msn.age, limits = c(-5,5), oob=squish) +  # Adjust fill scale if necessary
  geom_smooth(method='lm', color='black') +
  theme_minimal(base_family = "Gill Sans") +
  ylab('Feature age effect (t-value)') +
  xlab('MSN age effect (t-value)') +
  labs(color='', fill='') +
  theme(legend.position = 'none')

# Fig. 3, Panel D; SI. Fig. 18
ggsave(corvals.mm.msn, filename = paste0(outpath.SI,'corvals.mm.msn.png'),width = 6, height=2.8, bg='white')


## Co-location metabolism
sydnor = read.csv('data/external/statistical.maps.HCP.csv')
sydnor = sydnor[,-match(c('x','y','z','curvature','thickness','T1T2','CBF','PC1_AHBA'),colnames(sydnor))]
if(any(del.roi)){sydnor = sydnor[-del.roi,]}
mean.age.effect = rowMeans(cbind(str.l.sl.t[1:179],str.l.sl.t[180:358]))

perm.id.hemi = perm.id[1:179,]
df.sydnor = data.frame(cor = cor(sydnor, mean.age.effect), 
                       map = c('Arterial spin labelling','Cerebral metabolic rate 02','Cerebral metabolic rate glucose','Cerebral blood volume',
                               'Cerebral blood flow','Evolutionary expansion','Developmental expension','Aerobic glycolysis',
                               'Externopyramidisation','PC1 neurosynth','G1 fMRI','Allometric scaling','Geodesic distance'),
                       sign = cut(unlist(lapply(1:dim(sydnor)[2],function(r) cor.test(sydnor[,r],mean.age.effect)$p.value)), 
                                  breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")),
                       sign.spin= cut(unlist(lapply(1:dim(sydnor)[2], function(r) perm.sphere.p(sydnor[,r],mean.age.effect,perm.id.hemi,corr.type = 'spearman'))),
                                      breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))  # Create column of significance labels
df.sydnor$map = factor(df.sydnor$map, levels =  df.sydnor$map[order(df.sydnor$cor)])

prep.for.stats = df.sydnor[which(df.sydnor$sign.spin=='***' | df.sydnor$sign.spin=='*'),]
printout.stats[[8]] = paste(paste0(prep.for.stats[,2],' - Correlation: ',prep.for.stats[,1],';  - p-value: ',prep.for.stats[,4],' | '),collapse='')

ex1=data.frame(y=sydnor$glasser_CMRGlu,x=mean.age.effect) %>% 
  ggplot(aes(x=x,y=y,color=x))+geom_point()+geom_smooth(method='lm')+
  ylab('MR glucose')+
  xlab(expression(Delta~MS[14-26]))+
  labs(color=expression(Delta~MS[14-26]))+
  scale_color_gradientn(colours = colors.msn.age,limits=c(-4,4))+  
  theme_minimal(base_family = "Gill Sans")+ theme(text = element_text(size=12), legend.position = 'bottom')+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))
  
ex2=data.frame(y=sydnor$glasser_CBV,x=mean.age.effect) %>% 
  ggplot(aes(x=x,y=y,color=x))+geom_point()+geom_smooth(method='lm')+
  ylab('Cerebral BV')+
  xlab(expression(Delta~MS[14-26]))+
  labs(color=expression(Delta~MS[14-26]))+
  scale_color_gradientn(colours = colors.msn.age,limits=c(-4,4))+  
  theme_minimal(base_family = "Gill Sans")+ theme(text = element_text(size=12), legend.position = 'bottom')+
  guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Plot everything
sydnor.panel=ggplot(aes(y=map,x=1, fill=cor), data=df.sydnor)+geom_tile() + 
  scale_fill_gradientn(colours=colors.msn.age, limits = c(-0.5, 0.5), oob=squish) +
  geom_text(aes(label=sign.spin), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Correlation") + 
  theme_minimal(base_family = "Gill Sans")+ theme(axis.text.x=element_blank(),
                                                  axis.ticks.x=element_blank())#, 
# Fig. 3, Panel E
ggsave(ggarrange(sydnor.panel,ggarrange(ex1,ex2, ncol=1,common.legend = T, legend='bottom'),
          ncol=2,widths=c(1,1),
          vjust=0.5,hjust = 0,
          font.label = list(size = 12, family="Gill Sans"))+
  theme(plot.margin = margin(0.8,0,0,0, "cm")), 
  filename = paste0(outpath.main,'sydnor.png'),width = 6.5, height=3.5,bg='white')


# Sex effect on morphometric similarity
msn.sex.effect <- data.frame(label=rep(nm.ggseg,2),value=c(str.l.sl.sex.t,ifelse(str.l.sl.sex.p.fdr<0.05,str.l.sl.sex.t,NA)), 
                             group = factor(rep(c('non-corrected','FDR-corrected'),each = nroi))) %>%
  group_by(group) %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=value))+
  scale_fill_gradientn(colors=rev(paletteer_d('RColorBrewer::RdBu')),limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~group, ncol=1)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Sex (t-value)')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,))

# SI Fig. 8
ggsave(paste0(outpath.SI,'msn-sex-effect.png'),plot = msn.sex.effect,width = 10, height=4,bg='white')

# Site effect on morphometric similarity
str.l.site.p.fdr=p.adjust(str.l.site.p,method='fdr')
paste0('Regions that have significant MSN and site effects: ',nm[which(rowSums(cbind(str.l.site.p.fdr, str.l.sl.p.fdr)<0.05)==2)])

msn.site.effect <- data.frame(label=rep(nm.ggseg,2),value=c(str.l.site.F,ifelse(str.l.site.p.fdr<0.05,str.l.site.F,NA)), 
                              group = factor(rep(c('non-corrected','FDR-corrected'),each = nroi))) %>%
  group_by(group) %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=value))+
  scale_fill_gradientn(colors=rev(paletteer_d('RColorBrewer::RdBu')),limits=c(-15,15), oob=squish)+
  facet_wrap(~group, ncol=1)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Site (F-value)')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig. 9
ggsave(paste0(outpath.SI,'msn-site-effect.png'),plot = msn.site.effect,width = 10, height=4,bg='white')


# Baseline morphometric similarity
baseline.msn <- ggseg(.data=data.frame(label=nm.ggseg,value=str.l.14),mapping=aes(fill=value), atlas=glasser)+  
  scale_fill_gradientn(colors=colors.msn.str,limits=c(-0.1,0.1), oob=squish)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Age (t-value)')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig. 21
ggsave(paste0(outpath.SI,'baseline.MSN.png'),plot = baseline.msn,width = 10, height=4)


# Yeo network decoding
df.Yeo.age = data.frame(age.t = str.l.sl.t,class = as.factor(Yeo))
classmu.age = sapply(1:length(unique(Yeo)),function(class) mean(str.l.sl.t[Yeo==class],na.rm=TRUE))
yeo.name = tolower(yeo.name)
for (class in 1:length(unique(Yeo))) {
  df.Yeo.age$classmu[df.Yeo.age$class == class] = classmu.age[class]
}
# SI Fig 13, Panel A
panel.roiYeo <- ggplot(na.omit(df.Yeo.age), aes(x=class, y=age.t, fill = classmu)) +
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white") +
  scale_fill_gradientn(colours=colors.msn.age, limits = c(-1.5, 1.5), oob=squish) +
  ylab('MSN age-effect') + xlab('') +
  scale_x_discrete(labels=yeo.name) +
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.y = element_text(color =  yeo.col, size = 12,hjust=1),
        axis.title=element_text(size=12), legend.position = 'bottom')+
  labs(fill='Age (t-value)')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))+
  coord_flip()

# Von Economo cell type decoding
df.vonEconomo.age = data.frame(age.t = str.l.sl.t,class = as.factor(vonEconomo))
classmu.age = sapply(1:length(unique(vonEconomo)),function(class) mean(str.l.sl.t[vonEconomo==class],na.rm=TRUE))
for (class in 1:length(unique(vonEconomo))) {
  df.vonEconomo.age$classmu[df.vonEconomo.age$class == class] = classmu.age[class]
}
# SI Fig 13, Panel B
panel.roivonEconomo = ggplot(df.vonEconomo.age, aes(x=class, y=age.t, fill = classmu)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white") +
  scale_fill_gradientn(colours=colors.msn.age, limits = c(-3, 3), oob=squish) +
  coord_flip()+
  xlab('')+
  scale_x_discrete(labels=ve.name) +
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.y = element_text(color =  ve.col, size = 12,hjust=1),
        axis.title=element_text(size=12), legend.position = 'bottom')+
  ylab(expression(Delta~MS[14-26]))+
  labs(fill=expression(Delta~MS[14-26]))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig 13
ggsave(ggarrange(panel.roivonEconomo,panel.roiYeo,
                 ncol=2,vjust=-0.2,hjust = 0,common.legend=T, legend='bottom',
                 labels = c('A | Cytoarchitectonic classes',
                            'B | rsfMRI modules'),
                 font.label = list(size = 12, family="Gill Sans"))+
         theme(plot.margin = margin(1,0,0,0, "cm")),
       filename = paste0(outpath.SI,'yeo-VE.png'),width = 6, height=4,bg='white')


# Mesulam class specific changes
msn.net.str = array(NA,dim=c(length(unique(mesulam)),nroi,nsub))
for (n in 1:nsub) {
  for (v in 1:length(unique(mesulam))) {
    msn.net.str[v,,n] = rowMeans(MSN[,,n][,which(mesulam==unique(mesulam)[v])],na.rm=T) # cortical
  }
}

msn.net.str.p.fdr =msn.net.str.sl = msn.net.str.p = array(NA,dim=c(length(unique(mesulam)),nroi))
for (v in 1:length(unique(mesulam))) {
  print(v)
  for (r in 1:nroi) {
    try({
      df=data.frame(msn=msn.net.str[v,r,], age=MRI_table$age, sex= MRI_table$sex, site=MRI_table$site, id=MRI_table$id )
      l = lme(msn ~ age+sex+site, random=~1|id, data=na.omit(df))
      msn.net.str.sl[v,r] = summary(l)$tTable[2,4]
      msn.net.str.p[v,r] = summary(l)$tTable[2,5]
    })
  }
  msn.net.str.p.fdr[v,] = p.adjust(msn.net.str.p[v,], method='fdr')
}

df.net.msn.age = data.frame(msn=as.vector(t(msn.net.str.sl)), 
                            p=as.vector(t(msn.net.str.p)),
                            group=factor(rep(unique(mesulam),each=nroi)),
                            label=rep(nm.ggseg,length(unique(mesulam))))
df.net.msn.age$p.fdr=as.vector(apply(msn.net.str.p,1,function(x) p.adjust(x,method='fdr')))

# Network specific MSN change
# SI Fig. 15, Panel B, top
msn.net.plot.top = subset(data.frame(value=as.factor(mesulam), 
                                     group=factor(mesulam),
                                     label=nm.ggseg),group != 'cortical_wall' & group != 'no label') %>%
  group_by(group) %>%
  ggseg(atlas='glasser',mapping=aes(fill=value),hemisphere = 'left')+
  facet_grid(~group)+
  scale_fill_manual(values=mesulam.col)+
  theme_void(base_family = "Gill Sans")+
  theme(legend.position = 'none')

# SI Fig. 15, Panel B, bottom
msn.net.plot.bottom = subset(df.net.msn.age, group != 'cortical_wall' & group != 'no label') %>% 
  group_by(group) %>% 
  ggseg(atlas='glasser',mapping=aes(fill=msn),hemisphere = 'left')+
  facet_grid(~group)+
  scale_fill_gradientn(colours=colors.msn.age, limits=c(-5,5))+
  theme_void(base_family = "Gill Sans")+
  theme(legend.position = 'bottom',strip.text.x = element_blank())+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,title="Age (t-value)"))

# SI Fig. 15
ggsave(plot = ggarrange(msn.net.plot.top,
                        msn.net.plot.bottom,
                        nrow=2, vjust=0.5,hjust = 0,widths=c(1,1),
                        #labels = c('A | Metabolic enrichment',
                        #          ''),
                        font.label = list(size = 13, family="Gill Sans"))+
         theme(plot.margin = margin(0.5,0,0,0, "cm")),
       filename = paste0(outpath.SI,'msn.network.change.mesulam.png'), width=12,height=4, bg='white')


######################################################
#           SITE-SPECIFIC TRAJECTORIES
######################################################

# SI Fig 10, Panel A
SITE.ACQ = ggplot(acq.overview, aes(y=y,x=age, color=firstage))+
  facet_wrap(~site)+
  geom_point()+
  geom_line(aes(group=id,color=firstage))+
  labs(color='Baseline Age')+
  theme_minimal(base_family = 'Gill Sans')+
  ylab('Participant')+
  scale_color_viridis(option='magma')

# Estimate regional maps of rate of change in MSN degree for each site separately
dir.create(paste0(outpath.SI,'/sensitivity-analyses/'),recursive = T)
site.str.l.sl.p = site.str.l.sl.t = matrix(NA, nrow=3, ncol=nroi)
# Loops through sites
for (s in 1:3) { 
  site.idx = which(MRI_table$site ==s)
  # Loops over regions
  for (r in 1:nroi) {
    print(paste0('ROI: ',r))
    df.roi =  data.frame(id = MRI_table$id[site.idx], age = MRI_table$age[site.idx], 
                         sex = MRI_table$sex[site.idx], str =str[r,site.idx])
    
    # Linear Models
    lm.roi = tryCatch({ 
      lm.roi <- lme(str~age+sex, random = ~1|id,data = df.roi,na.action = na.omit) }
      , error = function(e) {an.error.occured <<- TRUE})
    
    if(!(class(lm.roi)=='lme')){next}else{
      site.str.l.sl.p[s,r] = summary(lm.roi)$tTable[2,5]     # p-value for effect of age
      site.str.l.sl.t[s,r] = summary(lm.roi)$tTable[2,4]     # t-value for effect of age
    }
  }
}

df.repl.site = data.frame(msn.t = c(site.str.l.sl.t[1,],site.str.l.sl.t[2,],site.str.l.sl.t[3,],str.l.sl.t), 
                          label = rep(nm.ggseg, 4),
                          site=rep(c('CBSU','UCL','WBIC','Original'),each=nroi)) 
df.repl.site$site = factor(df.repl.site$site, levels = c('Original','WBIC','UCL','CBSU'))

# SI Fig 10, Panel C
SITE.MAPS = df.repl.site %>%
  group_by(site)  %>%
  ggseg(atlas=glasser,mapping=aes(fill=msn.t))+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-limvar,limvar), oob=squish)+
  facet_wrap(~site, ncol=1)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  labs(fill='t-value')

combinations = combn(levels(df.repl.site$site),2)
SITE.CORR = data.frame(site1 = rep(NA, ncol(combinations)),
                       site2 = rep(NA, ncol(combinations)),
                       pval = rep(NA, ncol(combinations)),
                       rho = rep(NA, ncol(combinations)))
# Estimate correlation of site-specific with original maps
for (i in 1:ncol(combinations)) {
  print(i)
  map1 = subset(df.repl.site, site==combinations[1,i])$msn.t
  map2 = subset(df.repl.site, site==combinations[2,i])$msn.t
  rho = cor.test(map1,map2,method = 'spearman')$estimate
  pspin = perm.sphere.p(map1,map2, perm.id = perm.id)
  SITE.CORR$site1[i] = combinations[1,i]; SITE.CORR$site2[i] = combinations[2,i]
  SITE.CORR$pval[i] = pspin; SITE.CORR$rho[i] = rho; 
}

SITE.CORR$site1 = factor(SITE.CORR$site1,levels = c('Original','WBIC','UCL','CBSU'))
SITE.CORR$site2 = factor(SITE.CORR$site2,levels = c('Original','WBIC','UCL','CBSU'))
SITE.CORR$pvalprint = ifelse(SITE.CORR$pval<0.05 & SITE.CORR$pval>0.001,'P<0.05',paste0('P=',round(SITE.CORR$pval,2)))
SITE.CORR$pvalprint[which(SITE.CORR$pval<0.001)] = 'P<0.001'

# SI Fig 10, Panel B
SITE.CORRPLOT = SITE.CORR %>% 
  ggplot(aes(x=site1,y=site2, fill=rho, label=pvalprint))+
  geom_tile()+
  scale_fill_gradientn(colours = brewer.pal(9,'Reds'))+
  geom_text(aes(family='Gill Sans'))+
  theme_minimal(base_family='Gill Sans')+
  ylab('')+xlab('')+labs(fill='Spearman Corr.')

# Site-specific acquisition overview 
acq.overview=MRI_table %>% group_by(id) %>% mutate(firstscan=order(age), firstage=min(age))
relevel=acq.overview[acq.overview$firstscan==1,]$id[order(acq.overview[acq.overview$firstscan==1,]$age)]
acq.overview$id = factor(acq.overview$id, levels = relevel)
acq.overview$y = as.numeric(acq.overview$id)
acq.overview$site = factor(acq.overview$site, levels=c(1,2,3),labels=c('CBSU','UCL','WBIC'))


# SI Fig 10
ggsave(plot = ggarrange(ggarrange(plotlist = list(SITE.ACQ,SITE.CORRPLOT),ncol = 1),SITE.MAPS),
       filename = paste0(outpath.SI,'/sensitivity-analyses/site.specific.maps.png'),bg='white')


######################################################
#       WITHIN-SUBJECT CHANGES IN MSN DEGREE
######################################################
# Within-subject regional MSN change
TMPTABLE = MRI_table
STRTMP=data.frame(t(str))
colnames(STRTMP) = paste0('regional.',seq(1,nroi))
STRTMP$id = TMPTABLE$id
STRTMP$age_days = TMPTABLE$age_days
TMPTABLE = TMPTABLE %>% left_join(STRTMP)

REGBLFUDIFF=TMPTABLE %>%
  filter(session %in% c('ses-baseline', 'ses-1stFollowUp')) %>%
  group_by(id) %>%
  mutate(age_baseline = ifelse(session == "ses-baseline", age, NA)) %>%
  mutate(site_baseline = ifelse(session == "ses-baseline", site, NA)) %>%
  summarise(across(starts_with("regional."), ~ diff(.)),
            sex = first(sex),
            age_baseline = first(age_baseline),
            site_baseline = first(site_baseline))

REGBLFUDIFF.LONG = REGBLFUDIFF %>%
  pivot_longer(cols = starts_with("regional."),
               names_to = "variable",
               values_to = "value")

REGBLFUDIFF.LONG$label = factor(REGBLFUDIFF.LONG$variable, 
                                levels=paste0('regional.',seq(1,nroi)),
                                labels=nm.ggseg)

MUSUBJCHANGE = REGBLFUDIFF.LONG %>% group_by(label) %>% summarise(mu=mean(value, na.rm=T))

cor.test(MUSUBJCHANGE$mu, str.l.sl.t)
perm.sphere.p(MUSUBJCHANGE$mu,str.l.sl.t,perm.id=perm.id)

P1 = MUSUBJCHANGE %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=mu),position = 'stacked')+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-0.01,0.01),breaks=c(-0.01,0,0.01), oob=squish)+
  theme_void(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Mean within-person change in regional MSN degree')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

P2 = data.frame(within=MUSUBJCHANGE$mu,between=str.l.sl.t) %>% 
  ggplot(aes(x=between,y=within, color=between))+
  scale_color_gradientn(colors=colors.msn.age,limits=c(-5,5), oob=squish)+
  geom_point()+
  geom_smooth(method='lm')+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  ylab('Mean within-subject change')+xlab('Between-subject change')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

ggarrange(plotlist=list(P1,P2))


######################################################
#       RANK-STABILITY OF REGIONAL MSN DEGREE
######################################################
MRI_table$idx=1:nrow(MRI_table)
SUBSET = MRI_table[,c('id','sex','session','idx')]
SUBSET.WIDE = spread(SUBSET,session, idx)
colnames(SUBSET.WIDE)[c(3,5)] = c('ses1','ses2')
IDX = subset(SUBSET.WIDE, !is.na(ses1) & !is.na(ses2))

BLFUSET = MRI_table[c(IDX$ses1,IDX$ses2),]
regional.ranks.bl = apply(str[,IDX$ses1],1,rank)
regional.ranks.fu = apply(str[,IDX$ses2],1,rank)
LOC.RANK =LOC.RANK.P= vector(length=nroi)
for (r in 1:nroi) {
  TMP=BLFUSET
  TMP$rankROI = (c(regional.ranks.bl[,r],regional.ranks.fu[,r]))
  LOC.MSN.WIDE = spread(TMP[,c('id','session','rankROI')],session,rankROI)
  LOC.RANK[r] = cor.test(LOC.MSN.WIDE$`ses-baseline`,LOC.MSN.WIDE$`ses-1stFollowUp`,method='spearman')$estimate
  LOC.RANK.P[r] = cor.test(LOC.MSN.WIDE$`ses-baseline`,LOC.MSN.WIDE$`ses-1stFollowUp`,method='spearman')$p.value
}

# Estimate MSN global rank for each scan
MRI_table$glob.MSN = colMeans(str,na.rm=T)
# Reshape from long to wide 
BLFUSET.GLOB = MRI_table %>% group_by(session) %>% mutate(rankMSN = rank(glob.MSN),
                                                          rankGM = rank(glob.GM))
# MSN global rank dataframe
GLOBMSN.WIDE = pivot_wider(data = BLFUSET.GLOB,
                           id_cols = id,
                           names_from = session,
                           values_from = c(rankMSN, age,sex))
cor.test(GLOBMSN.WIDE$`rankMSN_ses-baseline`,GLOBMSN.WIDE$`rankMSN_ses-1stFollowUp`,method='spearman')

GLOBMSN.WIDE$agediff = GLOBMSN.WIDE$`age_ses-baseline`-GLOBMSN.WIDE$`age_ses-1stFollowUp`
lm(`rankMSN_ses-1stFollowUp`~ `rankMSN_ses-baseline`+agediff+`sex_ses-baseline`,data=GLOBMSN.WIDE)

# Correlation between baseline and follow up MSN rank
# SI Fig. 24, Panel A
plot.glob.rank = GLOBMSN.WIDE %>% 
  ggplot(aes(x=`rankMSN_ses-baseline`,y=`rankMSN_ses-1stFollowUp`))+
  geom_point()+
  theme_minimal(base_family='Gill Sans')+
  ylab('Global MSN rank at follow up')+xlab('Global MSN rank at baseline')+
  geom_smooth(method = 'lm')

# Regional correlation between baseline and follow up MSN degree rank
# SI Fig. 24, Panel B
plot.loc.rank = data.frame(label = rep(nm.ggseg,2),
           rankval = c(LOC.RANK,
                       ifelse(p.adjust(LOC.RANK.P, method='fdr')<0.05,LOC.RANK,NA)),
           group = rep(c('Uncorrected','FDR-corrected'),each=nroi)) %>%
  group_by(group) %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=rankval))+
  facet_wrap(~group, ncol=1)+
  scale_fill_gradientn(colors=pals::ocean.balance(20),limits=c(-0.5,0.5), oob=squish)+
  theme(legend.position = 'bottom')+
  labs(fill=expression(rho))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))+
  theme_minimal()

# SI Fig 24
ggsave(plot = ggarrange(plot.glob.rank,plot.loc.rank,widths = c(1,2)),height=3,
       filename = paste0(outpath.SI,'/sensitivity-analyses/rank-stability.png'),bg='white')


# Correlation between regional change in MSN degree and rank stability
df=data.frame(rankval = LOC.RANK,
           ranksign = p.adjust(LOC.RANK.P, method='fdr')<0.05,
           msnval = str.l.sl.t,
           msnsign = str.l.sl.p.fdr<0.05,
           label=nm.ggseg)
df$size=ifelse(df$ranksign, 2, 1)
df$shape = 20
df$shape[df$ranksign] = 17
df$shape[df$msnsign] = 21
df$shape[which(df$msnsign & df$ranksign)] = 24

# Correlation and spin test p-value
cor.test(df$rankval,df$msnval, method = 'spearman')
perm.sphere.p(df$rankval,df$msnval,perm.id = perm.id,corr.type = 'spearman')

# SI Fig 25, Panel A
plot.cor.rank.msn.age = df%>%
  ggplot(aes(x=msnval,y=rankval,size=size,fill=msnval,shape=as.factor(shape)))+
  geom_point()+
  scale_size_identity()+
  scale_shape_manual(values=c(17,20,21,24))+
  scale_fill_gradientn(colors = colors.msn.age, limits=c(-5,5))+
  theme_minimal(base_family='Gill Sans')+
  labs(fill='', shape='')+
  xlab('Change in MSN degree')+
  ylab('Rank correlation')

# SI Fig 25, Panel B, top
plot.cor.rank.msn.age.sign = df %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=ifelse(shape==24,msnval,NA)))+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-5,5), oob=squish)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Change in MSN degree')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig 25, Panel B, bottom
plot.cor.rank.msn.age.not.sign =df %>%
  ggseg(atlas=glasser,
        mapping=aes(fill=ifelse(shape==21,msnval,NA)))+
  scale_fill_gradientn(colors=colors.msn.age,limits=c(-5,5), oob=squish)+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'bottom')+
  labs(fill='Change in MSN degree')+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))


# SI Fig 25
ggsave(ggarrange(plot.cor.rank.msn.age,
                 ggarrange(plot.cor.rank.msn.age.sign,plot.cor.rank.msn.age.not.sign, 
                           ncol = 1,common.legend = T),widths = c(1.2,1)), 
       file=paste0(outpath.SI,'sensitivity-analyses/rank.vs.msn.age.png'), bg='white',height = 3)

#########################################################################################
###                          NEUROSYNTH fMRI DECODING                                ####
#########################################################################################
# Change this path to where ever you have downloaded the Neurosynth terms to
neurosynth=read.csv('data/neurosynth/neurosynth.thresholded.csv')

# SI Fig. 14, Panel A
words.pos = subset(neurosynth, corr.>0) # choosing only positive weights
words.pos = words.pos[order(words.pos$corr.,decreasing = T),]
words.pos=words.pos[1:50,]
freq.pos = words.pos$corr./sum(words.pos$corr.)*100 # transform into frequency
pdf(paste0(outpath.SI,'wordcloud_pos.thresh.pdf')) # save image as pdf file
wordcloud(words=words.pos$term, 
          freq=freq.pos, 
          min.freq = min(freq.pos), 
          random.order = FALSE,
          scale=c(2,.1), 
          colors = colors.msn.age[10])
dev.off()

# SI Fig. 14, Panel B
words.neg = subset(neurosynth, corr.<0) # choosing only positive weights
words.neg = words.neg[order(words.neg$corr.),]
words.neg=words.neg[1:50,]
freq.neg = words.neg$corr./sum(words.neg$corr.)*100 # transform into frequency
pdf(paste0(outpath.SI,'wordcloud_neg_thresh.pdf')) # save image as pdf file
wordcloud(words=words.neg$term, 
          freq=freq.neg, 
          #  max.words = nWords, 
          min.freq = min(freq.neg), 
          random.order = FALSE,
          scale=c(2,.1), 
          colors = colors.msn.age[2])
dev.off()

df.neurosynth = data.frame(t=ifelse(str.l.sl.p.fdr<0.05, str.l.sl.t, NA),
           label=nm.ggseg) 

subset(df.neurosynth, t<0) %>%
  ggseg(atlas=glasser,mapping=aes(fill=t), position='stacked')+
  scale_fill_gradientn(colors=colors.msn.age, limits=c(-5,5))+
  theme_void(base_family = "Gill Sans")+
  labs(fill='t-value')+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

subset(df.neurosynth, t>0) %>%
  ggseg(atlas=glasser,mapping=aes(fill=t), position='stacked')+
  scale_fill_gradientn(colors=colors.msn.age, limits=c(-5,5))+
  theme_void(base_family = "Gill Sans")+
  labs(fill='t-value')+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


sink(paste0(outpath.main,'all.stats.txt'))
print(printout.stats)
sink()



