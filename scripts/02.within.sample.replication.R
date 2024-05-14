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
library(extrafont)    # to use Gill Sans font, first time use: library(extrafont); font_import(); loadfonts()
library(matrixStats)  # matrix manipulation
library(Hmisc)        # compute correlation on matrix rows 
library(scales)       # squish function
library(xtable)       # saving tables for latex doc

#### Load dependencies
source('scripts/dependencies/msn.age.perm.R')                             # Within-sample-replication permutation

#### Load data
load(paste0('data/for.further.analyses.RData'))

#################################################################
#                  WITHIN-SAMPLE REPLICATION
#################################################################
# Permutation testing for sensitivity analysis
## Setup permutation matrix
nperm <- 1000
MAT <- matrix( NA_integer_, nrow=NROW(MRI_table), ncol=nperm )
for(IDX in 1:NCOL(MAT) ) {
  for( sex in c('female','male')) {
    for( grp in 1:5 ) {
      WHICH <- which( MRI_table$sex==sex & MRI_table$agebin==grp)
      N= round(9/10*length(WHICH))
      BIN = sample(WHICH,N, replace=F)
      MAT[ BIN, IDX ] = 1
    }
  }
}

# Careful, this takes a *long* time to run
T.PERM = SL.PERM = matrix(NA,nrow=nroi,ncol=nperm)
for (i in 1:NCOL(MAT)) {
  print(paste0('Permation: ',i))
  idx=which(!is.na(MAT[,i]))
  perm.out=msn_age_perm(MRI_table[idx,],str[,idx])
  T.PERM[,i] = perm.out$t
  SL.PERM[,i] = perm.out$sl
}

# Estimate standard deviation across permuations as confidence intervals
SD.SL.PERM = rowSds(SL.PERM)
CI.SL = cbind(str.l.sl.beta+(-1.96*SD.SL.PERM),
              str.l.sl.beta+1.96*SD.SL.PERM)

# Plot high/low confidence interval maps
df.perm.maps = data.frame(values = c(CI.SL[,1], CI.SL[,2],str.l.sl),
                          group = factor(rep(c('CI Low','CI High','Real'),each=nroi), levels=c('CI Low','Real','CI High')),
                          label=rep(nm.ggseg,3))

# SI Fig. 11, Panel C
lower.upper.map.plot = df.perm.maps%>%  
  group_by(group) %>%
  ggseg(atlas=glasser,mapping=aes(fill=values))+
  facet_wrap(.~group, nrow=3)+
  labs(fill='Rate of Change')+
  theme_minimal(base_family = "Gill Sans")+
  scale_fill_gradientn(colors=colors.msn.age, limits= c(-0.0025,0.0025),oob=squish)

# Correlation between original and lower/upper confidence interval maps
df.perm.maps.scatter = data.frame(repl = c(CI.SL[,1], CI.SL[,2]),
                                  orig = rep(str.l.sl,2),
                                  group = factor(rep(c('CI Low','CI High'),each=nroi), levels=c('CI Low','CI High')),
                                  label=rep(nm.ggseg,2),
                                  thresh=rep(str.l.sl.p.fdr<0.05,2))
# SI Fig. 11, Panel D
lower.upper.plot = df.perm.maps.scatter %>% ggplot(aes(x=orig, y=repl, group=group, color=repl))+
  geom_point(aes(size=thresh))+geom_smooth(method='lm', color='black')+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = colors.msn.age, limits=c(-0.0025,0.0025))+
  theme_minimal(base_family = "Gill Sans")+
  ylab('Lower and upper bound')+
  labs(size='significant')

# Ablation analysis. Leave out subsequently larger amounts of data
THRESH = c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5)
nthresh <- length(THRESH)
MATTHRESH <- matrix( NA_integer_, nrow=NROW(MRI_table), ncol=nthresh )
for(IDX in 1:NCOL(MATTHRESH) ) {
  for( sex in c('female','male')) {
    for( grp in 1:5 ) {
      WHICH <- which( MRI_table$sex==sex & MRI_table$agebin==grp)
      N= round(THRESH[IDX]*length(WHICH))
      BIN = sample(WHICH,N, replace=F)
      MATTHRESH[ BIN, IDX ] = 1
    }
  }
}

write.table(MATTHRESH,paste0(outpath.SI,'permutation.matrix.txt'),row.names = FALSE,col.names = FALSE)

# Re-estimate MSN degree change based on smaller samples
T.THRESH = SL.THRESH = matrix(NA,nrow=nroi,ncol=nthresh)
for (i in 1:NCOL(MATTHRESH)) {
  print(paste0('Threshold: ',i))
  idx=which(!is.na(MATTHRESH[,i]))
  perm.out=msn_age_perm(MRI_table[idx,],str[,idx])
  T.THRESH[,i] = perm.out$t
  SL.THRESH[,i] = perm.out$sl
}

df.thresh.flat = data.frame(vals = as.vector(T.THRESH), 
                            map = rep(paste0(100*THRESH,'%'),each=nroi),
                            labels = rep(nm.ggseg,NCOL(MATTHRESH)),
                            orig = rep(str.l.sl,NCOL(MATTHRESH)),
                            sig = rep((str.l.sl.p.fdr<0.05),NCOL(MATTHRESH)))

df.thresh.cor.flat = data.frame(vals = as.vector(cor(SL.THRESH)),
                                X = rep(paste0(100*(THRESH),'%'),NCOL(MATTHRESH)),
                                Y = rep(paste0(100*(THRESH),'%'),each=NCOL(MATTHRESH)))

# SI Fig 11, Panel A
# Correlation between original and ablation analysis maps
rand.heatmap.plot = ggplot(df.thresh.cor.flat,aes(x=X, y=Y, fill=vals))+geom_tile()+
  scale_fill_gradientn(colours = paletteer_c("viridis::magma",100)) +
  scale_x_discrete(limits=rev,position = "top")+
  theme_minimal(base_family = "Gill Sans")+
  theme(legend.position = 'right')+
  labs(fill='Correlation')+
  xlab('')+ylab('')

# Correlation between original and smaller sample maps
# SI Fig 11, Panel B
rand.plot = ggplot(df.thresh.flat, aes(x=orig, y=vals, group=map))+
  geom_point(alpha=0.5, aes(color=map, size=as.factor(sig)))+
  geom_smooth(method='lm',aes(color=map))+
  scale_color_manual(values = paletteer_c("viridis::viridis",10))+
  theme_minimal(base_family = "Gill Sans")+
  geom_hline(yintercept=0)+
  labs(color='threshold', size='significant')+
  xlab('Original')+ylab('Random at various thresholds')


# SI Fig 11
ggsave(ggarrange(ggarrange(rand.heatmap.plot,rand.plot,ncol=2,nrow=1,widths = c(1,1), heights=c(0.6,0.9))+
                   theme(plot.margin = margin(0.5,0,0,0, "cm")),
                 ggarrange(lower.upper.map.plot, lower.upper.plot, widths = c(1.8,1.2), heights = c(1,0.6))+
                   theme(plot.margin = margin(0.5,0,0,0, "cm")),
                 vjust=0,hjust = 0,nrow=2, ncol=1,
                 font.label = list(size = 12, family="Gill Sans")),
       filename = paste0(outpath.SI,'sensitivity-analyses/within.sample.replication.png'),width=9, height=7, bg='white')

