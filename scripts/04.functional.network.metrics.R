#### Load required packages
library(igraph)       # network abalysis
library(nlme)         # linar mixed effects models
library(tidyr)        # data wrangling
library(ggplot2)      # plotting
library(ggseg)        # brain surface plotting
library(ggsegGlasser) # brain surface plotting
library(dplyr)        # data wrangling
library(forcats)      # fct_reorder
library(ggpubr)       # arrange plots


#### Load dependencies
source('scripts/dependencies/network.metrics.R')                          # Graph property estimation
source('scripts/external/rotate_parcellation-master/R/perm.sphere.p.R')   # Spin tests


#### Load dependencies
load('data/for.further.analyses.RData')

#########################################################################################
###                        FUNCTIONAL NETWORK METRICS                                ####
#########################################################################################
set.seed(993)

# Create MSN group level matrix to derive louvain communities
MSN.group = apply(MSN,c(1,2),function(i) mean(i, na.rm=TRUE))
MSN.group.pos = graph.adjacency(MSN.group, mode = "undirected", weighted = TRUE, diag = FALSE)
nedge=nroi*(nroi-1)/2
MSN.group.pos =  delete.edges(MSN.group.pos, which(rank(E(MSN.group.pos)$weight)<(0.85*nedge)))
louvain.pos = cluster_louvain(MSN.group.pos)

all.fc.graphs = list()
fc.matched = FC[-del.roi,-del.roi,]
for(s in 1:nsub){
  fc.pos = fc.matched[,,s]; fc.pos[fc.pos<0] = NA
  all.fc.graphs[[s]] <- graph.adjacency(fc.pos, mode = "undirected", weighted = TRUE, diag = FALSE)
  all.fc.graphs[[s]] = delete.edges(all.fc.graphs[[s]], which(is.na(E(all.fc.graphs[[s]])$weight)))
  all.fc.graphs[[s]]$id <- MRI_table$id[s]
}
all.fc.graphs <- assign_communities(all.fc.graphs, louvain.pos, "community")

fc.network.metrics = list()
for (s in 1:length(all.fc.graphs)) {
  print(s)
  fc.network.metrics[[s]] = calc_nodal_metric(all.fc.graphs[[s]],'community')
  #fc.network.metrics[[s]] = calc_nodal_metric(all.fc.graphs[[s]],'mesulam')
}

weighted.degree.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$strength))
part.coef.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$part.coeff))
local.clustering.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$local.clustering))
efficiency.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$efficiency))
betweenness.centrality.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$betweenness.centrality))
between.module.deg.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$between.module.deg))
within.module.deg.all=do.call(rbind,lapply(fc.network.metrics, function(l) l$within.module.deg))

WEIGHTED.DEGREE.AGE.T = WEIGHTED.DEGREE.AGE.P = WEIGHTED.DEGREE.BL = array(NA, dim = c(nroi,1))
PARTCOEF.AGE.T = PARTCOEF.AGE.P = PARTCOEF.BL = array(NA, dim = c(nroi,1))
CLUST.AGE.T = CLUST.AGE.P = CLUST.BL = array(NA, dim = c(nroi,1))
EFFICIENCY.AGE.T = EFFICIENCY.AGE.P = EFFICIENCY.BL = array(NA, dim = c(nroi,1))
BETWEENNESS.AGE.T = BETWEENNESS.AGE.P = BETWEENNESS.BL = array(NA, dim = c(nroi,1))
WITHIN.AGE.T = WITHIN.AGE.P = WITHIN.BL = array(NA, dim = c(nroi,1))
BETWEEN.AGE.T = BETWEEN.AGE.P = BETWEEN.BL = array(NA, dim = c(nroi,1))
for (r in 1:nroi) {
  try({
    df.part = data.frame(y=part.coef.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.part = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.part))
    PARTCOEF.AGE.T[r]=summary(lm.part)$tTable[2,4]
    PARTCOEF.AGE.P[r]=summary(lm.part)$tTable[2,5]
    PARTCOEF.BL[r]=summary(lm.part)$tTable[1,1]+14*summary(lm.part)$tTable[2,1]+1/2*summary(lm.part)$tTable[3,1]+1/3*summary(lm.part)$tTable[4,1]+1/3*summary(lm.part)$tTable[5,1]
  }, silent = T)
  try({
    df.clust = data.frame(y=local.clustering.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.clust = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.clust))
    CLUST.AGE.T[r]=summary(lm.clust)$tTable[2,4]
    CLUST.AGE.P[r]=summary(lm.clust)$tTable[2,5]
    CLUST.BL[r]=summary(lm.clust)$tTable[1,1]+14*summary(lm.clust)$tTable[2,1]+1/2*summary(lm.clust)$tTable[3,1]+1/3*summary(lm.clust)$tTable[4,1]+1/3*summary(lm.clust)$tTable[5,1]
  }, silent = T)
  try({
    df.eff = data.frame(y=efficiency.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.eff = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.eff))
    EFFICIENCY.AGE.T[r]=summary(lm.eff)$tTable[2,4]
    EFFICIENCY.AGE.P[r]=summary(lm.eff)$tTable[2,5]
    EFFICIENCY.BL[r]=summary(lm.eff)$tTable[1,1]+14*summary(lm.eff)$tTable[2,1]+1/2*summary(lm.eff)$tTable[3,1]+1/3*summary(lm.eff)$tTable[4,1]+1/3*summary(lm.eff)$tTable[5,1]
  }, silent = T)
  try({
    df.betw = data.frame(y=betweenness.centrality.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.betw = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.betw))
    BETWEENNESS.AGE.T[r]=summary(lm.betw)$tTable[2,4]
    BETWEENNESS.AGE.P[r]=summary(lm.betw)$tTable[2,5]
    BETWEENNESS.BL[r]=summary(lm.betw)$tTable[1,1]+14*summary(lm.betw)$tTable[2,1]+1/2*summary(lm.betw)$tTable[3,1]+1/3*summary(lm.betw)$tTable[4,1]+1/3*summary(lm.betw)$tTable[5,1]
  }, silent = T)
  try({
    df.within = data.frame(y=within.module.deg.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.within = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.within))
    WITHIN.AGE.T[r]=summary(lm.within)$tTable[2,4]
    WITHIN.AGE.P[r]=summary(lm.within)$tTable[2,5]
    WITHIN.BL[r]=summary(lm.within)$tTable[1,1]+14*summary(lm.within)$tTable[2,1]+1/2*summary(lm.within)$tTable[3,1]+1/3*summary(lm.within)$tTable[4,1]+1/3*summary(lm.within)$tTable[5,1]
  }, silent = T)
  try({
    df.between = data.frame(y=between.module.deg.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.between = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.between))
    BETWEEN.AGE.T[r]=summary(lm.between)$tTable[2,4]
    BETWEEN.AGE.P[r]=summary(lm.between)$tTable[2,5]
    BETWEEN.BL[r]=summary(lm.between)$tTable[1,1]+14*summary(lm.between)$tTable[2,1]+1/2*summary(lm.between)$tTable[3,1]+1/3*summary(lm.between)$tTable[4,1]+1/3*summary(lm.between)$tTable[5,1]
  }, silent = T)
  try({
    df.weighted.degree = data.frame(y=weighted.degree.all[,r],age=MRI_table$age,sex=MRI_table$sex, id = MRI_table$id, center=MRI_table$site)
    lm.weighted.degree = lme(y ~ age + sex + center, random=~1|id,data=na.omit(df.weighted.degree))
    WEIGHTED.DEGREE.AGE.T[r]=summary(lm.weighted.degree)$tTable[2,4]
    WEIGHTED.DEGREE.AGE.P[r]=summary(lm.weighted.degree)$tTable[2,5]
    WEIGHTED.DEGREE.BL[r]=summary(lm.weighted.degree)$tTable[1,1]+14*summary(lm.weighted.degree)$tTable[2,1]+1/2*summary(lm.weighted.degree)$tTable[3,1]+1/3*summary(lm.weighted.degree)$tTable[4,1]+1/3*summary(lm.weighted.degree)$tTable[5,1]
  }, silent = T)
}

cor.test(BETWEEN.AGE.T,str.l.sl.t, na.rm=T) #0.08
cor.test(WITHIN.AGE.T,str.l.sl.t, na.rm=T) # 0.12
cor.test(BETWEENNESS.AGE.T,str.l.sl.t, na.rm=T)# 0.02
cor.test(EFFICIENCY.AGE.T,str.l.sl.t, na.rm=T)# -0.12
cor.test(CLUST.AGE.T,str.l.sl.t, na.rm=T) #0.02
cor.test(PARTCOEF.AGE.T,str.l.sl.t, na.rm=T)# -0.23
cor.test(WEIGHTED.DEGREE.AGE.T,str.l.sl.t) # 0.049

ALL.T = list(BETWEEN.AGE.T,WITHIN.AGE.T,EFFICIENCY.AGE.T,BETWEENNESS.AGE.T,CLUST.AGE.T,PARTCOEF.AGE.T,WEIGHTED.DEGREE.AGE.T)

all.corr.vals = data.frame(p=unlist(lapply(ALL.T, function(i) cor.test(i,str.l.sl.t, na.rm=T)$p.value)),
                           t = unlist(lapply(ALL.T, function(i) cor.test(i,str.l.sl.t, na.rm=T)$estimate)),
                           p.spin = unlist(lapply(ALL.T,function(i) perm.sphere.p(as.vector(i),str.l.sl.t,perm.id,corr.type = 'spearman'))))


# Age effect on MSN degree vs functional network metrics
all.corr.vals$p.spin.fdr=p.adjust(all.corr.vals$p.spin,method='fdr')
all.corr.vals$mm = c('between network','within network','efficiency',
                     'betweenness centrality','clustering coefficient',
                     'participation coefficient','weighted degree')
# Fig 4, Panel A
net.cor = all.corr.vals %>% 
  mutate(mm = fct_reorder(mm,t)) %>%
  ggplot(aes(y=mm,x=t, fill=t))+geom_col()+
  xlim(-0.3,0.3)+
  scale_fill_gradientn(colours=colors.msn.age,lim=c(-0.3,0.3))+
  ylab('')+
  labs(fill='r')+
  xlab('Pearson correlation')+
  geom_point(shape=8,size=5,aes(x=ifelse(p.spin<0.05,t,NA),y=mm))+
  theme_minimal(base_family = "Gill Sans")


#########################################################################################
###                         PARTICIPATION COEFFICIENT                                ####
#########################################################################################

p.spin.partcoef.msn = perm.sphere.p(PARTCOEF.AGE.T,str.l.sl.t,perm.id,corr.type = 'spearman')
df.partcoef = data.frame(t=PARTCOEF.AGE.T[,1], 
                         p=p.adjust(PARTCOEF.AGE.P[,1], method='fdr'), 
                         label=nm.ggseg,
                         mesulam=mesulam)

# Surface plot of adolescent changes in participation coefficient
# Fig 4, Panel 
partcoeffc.pl = df.partcoef %>%
  ggseg(atlas=glasser,mapping=aes(fill=t), position='stacked')+
  scale_fill_gradientn(colors=c('darkblue','white','darkred'), limits=c(-3,3))+
  theme_void(base_family = "Gill Sans")+
  labs(fill='t-value')+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Adolescent changes in participation coefficient vs in MSN degree
# Fig 4, Panel D
scatter.partcoef.msn = data.frame(partcoef=PARTCOEF.AGE.T[,1],
                                  msn=str.l.sl.t) %>%
  ggplot(aes(x=msn,y=partcoef, color=msn))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_gradientn(colours=colors.msn.age)+
  theme_minimal(base_family = "Gill Sans")+
  labs(color='t-value')+
  xlab(expression(t[age]~MSN))+
  ylab(expression(t[age]~FC~participation~coefficient))

# Fig 4
ggsave(ggarrange(net.cor,partcoeffc.pl,scatter.partcoef.msn,
                 nrow=1, ncol=3,widths = c(1.7,1,1.2),
                 vjust=0,hjust = 0,
                 labels = c('A | FC vs MSN network changes',
                            'C | Age changes in PC',
                            'D | Age effect in PC vs MSN'),
                 font.label = list(size = 12, family="Gill Sans"))+
         theme(plot.margin = margin(0.5,0,0,0, "cm")),
       filename = paste0(outpath.main,'partcoef.png'),width=13, height=3, bg='white')

# Participation by mesulam zone
# SI Fig 23
partcoef.mesulam = subset(df.partcoef, mesulam != 'cortical_wall' & mesulam != 'no label') %>% 
  group_by(mesulam) %>% 
  mutate(mean_t=mean(t,na.rm=T)) %>% 
  ggplot(aes(x=mesulam, y=t,fill=mean_t)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white")+
  scale_fill_gradientn(colours=c('darkblue','white','darkred'), limits = c(-2, 2), oob=squish) +
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.x = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3"),angle = 90),
        axis.title=element_text(size=12),text = element_text(size=12), legend.position = 'bottom')+
  ylab("Rate of Change in Participation Coefficient")+
  labs(fill="Rate of Change in Participation Coefficient")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

ggsave(partcoef.mesulam,file = paste0(outpath.SI,'partcoef.mesulam.png'),width=3,height=5)

out.latex.partcoef.mesulam = compare_means(t~mesulam,subset(df.partcoef, mesulam != 'cortical_wall' & mesulam != 'no label') )
out.latex.partcoef.mesulam = out.latex.partcoef.mesulam[,c(2:5,7)]
print(xtable(out.latex.partcoef.mesulam, type = "latex",), file = paste0(outpath.SI,'partcoef.mesulam.tex'),include.rownames=FALSE)


# Rate of change in weighted degree
# SI Fig 22, Panel A
weighted.degree.pl = data.frame(t=WEIGHTED.DEGREE.AGE.T[,1], 
                                p=p.adjust(WEIGHTED.DEGREE.AGE.P[,1], method='fdr'), 
                                label=nm.ggseg) %>%
  ggseg(atlas=glasser,mapping=aes(fill=t), position='stacked')+
  scale_fill_gradientn(colors=c('darkblue','white','darkred'), limits=c(-5,5))+
  theme_void(base_family = "Gill Sans")+
  labs(fill='t-value')+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Rate of change in weighted degree vs rate of change in MSN degree
# SI Fig 22, Panel B
scatter.weighted.degree.msn = data.frame(partcoef=WEIGHTED.DEGREE.AGE.T[,1],
                                         msn=str.l.sl.t) %>%
  ggplot(aes(x=msn,y=partcoef, color=msn))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_gradientn(colours=colors.msn.age)+
  theme_minimal(base_family = "Gill Sans")+
  labs(color='t-value')+
  xlab(expression(t[age]~MSN))+
  ylab(expression(t[age]~FC~weighted~degree))

# SI Fig 22
ggsave(ggarrange(weighted.degree.pl,scatter.weighted.degree.msn,
                 nrow=1, ncol=2,
                 vjust=0,hjust = 0,
                 labels = c('A | Age effect on FC weighted degree',
                            'B | MSN age vs FC weighted degree age'),
                 font.label = list(size = 12, family="Gill Sans"))+
         theme(plot.margin = margin(0.5,0,0,0, "cm")),
       filename = paste0(outpath.main,'FCdegree.png'),width=8, height=3, bg='white')




