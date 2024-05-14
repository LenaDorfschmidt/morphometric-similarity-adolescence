#### Load required packages

#### Load dependencies


#### Load data
load(paste0('data/for.further.analyses.RData'))

#### Other setup
colors.coupling.bl = paletteer_c("scico::cork",n=100)
colors.coupling.age = paletteer_d('Redmonder::dPBIRdBu')


######################################################
#           STRUCTURE-FUNCTION-COUPLING
######################################################

# Global structure=function coupling
triup = upper.tri(matrix(nrow=nroi,ncol=nroi))
fc=FC[-del.roi,-del.roi,]
fc.edge = err = phenomat = array(NA, dim=nsub) 
for (s in 1:nsub) {
  tryCatch({
    fc.edge[s] = cor.test(MSN[,,s][triup],fc[,,s][triup],na.action=na.rm(),method='spearman')$estimate
  },
  error=function(e) {
    err[s] = TRUE
  })
}

lm.edge=lme(fcedge~age+sex+center, random=~1|id, data=df.str.fc, na.action=na.exclude)
summary(lme(fcstr~age+sex+center, random=~1|id, data=df.str.fc, na.action=na.exclude))$tTable
printout.stats[[9]] = paste0('Effect of age on global structure-function coupling - t = ',summary(lm.edge)$tTable[2,4],'; t = ',summary(lm.edge)$tTable[2,4])

# Fig. 4, Panel C
struc.fc.age=ggplot(df.str.fc,aes(x=age, y=fcstr, color=sex))+
  geom_point(alpha=0.7)+
  geom_line(aes(group=id),alpha=0.3) +
  scale_color_manual(values=c(colors.morph[7],colors.morph[1]))+
  geom_smooth(method='lm', color='black')+
  theme_minimal(base_family = "Gill Sans")+xlab('Age')+ylab('Structure-Function Correlation')+
  labs(color='')+
  theme(legend.position = 'bottom')

# Edgewise FC
str.fc.sim = str.fc.sim.r2=  matrix(NA,nrow=dim(fc)[1], ncol=nsub)
for (s in 1:nsub) {
  for (r in 1:dim(fc)[1]) {
    tmp=try({cor.test(MSN[r,-r,s],fc[r,-r,s], method='spearman', na.action=na.rm())$estimate})
    if(class(tmp)!="try-error"){
      str.fc.sim[r,s] = tmp
    } else{
      str.fc.sim[r,s] = NA
    }
  }
}

t.fc.str=p.fc.str=bl.fc.str=sl.fc.str=array(NA, dim=c(dim(fc)[1]))
for (r in 1:dim(fc)[1]) {
  try({
    df = data.frame(edge = str.fc.sim[r,],age = MRI_table$age, sex = MRI_table$sex, site= MRI_table$site, id = MRI_table$id)
    lm = lme(edge~age+sex+site, random = ~1|id,data = df, na.action = na.exclude)
    t.fc.str[r]=summary(lm)$tTable[2,4]
    p.fc.str[r]=summary(lm)$tTable[2,5]
    bl.fc.str[r]=lm$coefficients$fixed[1] + lm$coefficients$fixed[2]*14 + 1/2 * lm$coefficients$fixed[3] + lm$coefficients$fixed[4]*(1/3) + lm$coefficients$fixed[5]*(1/3)
    sl.fc.str[r]=lm$coefficients$fixed[2]
  })
}

df.fc.struct.bl.delta = data.frame(bl=bl.fc.str,
                                   sl=sl.fc.str,
                                   t.coupl=t.fc.str,
                                   mean.coupl.strength = rowMeans(str.fc.sim, na.rm=T),
                                   label=nm.ggseg,
                                   msn=str.l.sl.t,
                                   msn.bl=str.l.14,
                                   p.fdr=p.adjust(p.fc.str, method='fdr'),
                                   class = as.factor(mesulam))

printout.stats[[10]] = paste0('N regions with age on local structure-function coupling p_FDR < 0.05 = ',sum(df.fc.struct.bl.delta$p.fdr<0.05, na.rm=TRUE))
printout.stats[[11]] = paste0('Signifianct regions: ',paste(df.fc.struct.bl.delta$label[which(df.fc.struct.bl.delta$p.fdr<0.05)],collapse = '; '))

# Baseline structure-function coupling
# Fig. 4 Panel D, left
bl.coupling = df.fc.struct.bl.delta %>% 
  ggseg(atlas=glasser,mapping=aes(fill=bl), position = 'stacked')+
  scale_fill_gradientn(colors=magma(20), oob=squish) +
  theme_void(base_family = "Gill Sans")+
  ylab('')+xlab('')+
  labs(fill='Baseline Coupling')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  theme(legend.position = 'bottom')

# Rate of change in structure-function coupling
# Fig. 4 Panel D, right
sl.coupling = df.fc.struct.bl.delta %>% 
  ggseg(atlas=glasser,mapping=aes(fill=t.coupl), position = 'stacked')+
  scale_fill_gradientn(colors=rev(colors.coupling.age), limits=c(-5,5), breaks=c(-5,0,5), oob=squish) +
  theme_void(base_family = "Gill Sans")+
  ylab('')+xlab('')+
  labs(fill='Coupling Rate of Change')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
  theme(legend.position = 'bottom')

# Rate of change in structure-function coupling vs rate of change in MSN degree
# Fig. 4 Panel E
p.spin.couplingDelta.vs.msnDelta = perm.sphere.p(sl.fc.str[!is.na(bl.fc.str)],str.l.sl.t[!is.na(bl.fc.str)],perm.id.fc,corr.type = 'spearman')
r.couplingDelta.vs.msnDelta = cor.test(str.l.sl.t,sl.fc.str)$estimate
printout.stats[[12]] = paste0('Correlation coupling rate of change with MSN rate of change: ',r.couplingDelta.vs.msnDelta,'; p-value = ',p.spin.couplingDelta.vs.msnDelta)
msnDelta.vs.couplingDelta = df.fc.struct.bl.delta %>% ggplot(aes(x=msn,y=sl,color=msn))+
  geom_point()+geom_smooth(method='lm', color='black') +
  theme_minimal(base_family = "Gill Sans")+
  xlab('Change in MSN degree')+
  ylab('Coupling Rate of Change')+
  geom_hline(yintercept = 0,linetype='dashed')+
  geom_vline(xintercept = 0,linetype='dashed')+
  scale_color_gradientn(colors=colors.msn.age, limits=c(-5,5), oob=squish) +
  labs(color='Change in MSN degree')+
  theme(legend.position = 'bottom')+
  guides(color = guide_colourbar(title.position="top", title.hjust = 0.5))


# Looks funny, but all the ggarranges make sure the legend is in the same position
# Fig. 4 
ggsave(
  ggarrange(ggarrange(struc.fc.age, common.legend = T,legend = 'bottom'),
            NULL,
            ggarrange(bl.coupling,sl.coupling, nrow=1,common.legend = F, legend = 'bottom'),
            NULL,
            ggarrange(msnDelta.vs.couplingDelta, nrow=1, common.legend = T, legend = 'bottom'),
            nrow=1, widths = c(1,0.1,2,0.1,1.05),
            vjust=0,hjust = 0,
            labels = c('A | Global coupling',
                       '',
                       'B | Parameters of adolescent change in coupling',
                       '',
                       'C | Coupling rate of change'),
            #  'D | Relative Change '),
            font.label = list(size = 12, family="Gill Sans"))+
    theme(plot.margin = margin(0.5,0,0,0, "cm")),
  filename = paste0(outpath.main,'Fig4.png'),width=8, height=3, bg='white')


# Structure-function coupling at baseline by Mesulam zone 
# SI Fig. 19 Panel A
bl.coupl.mesulam = subset(df.fc.struct.bl.delta, class != 'cortical_wall' & class != 'no label') %>% 
  group_by(class) %>% 
  mutate(mean_t=mean(bl,na.rm=T)) %>% 
  ggplot(aes(x=class, y=bl,fill=mean_t)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white")+
  scale_fill_gradientn(colours=magma(20), limits = c(-0.1, 0.3), oob=squish) +
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.x = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3"),angle = 90),
        axis.title=element_text(size=12),text = element_text(size=12), legend.position = 'bottom')+
  ylab("Baseline Coupling")+
  labs(fill="Baseline Coupling")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# Rate of change in structure-function coupling by Mesulam zone 
# SI Fig. 19 Panel B
sl.coupl.mesulam = subset(df.fc.struct.bl.delta, class != 'cortical_wall' & class != 'no label') %>% 
  group_by(class) %>% 
  mutate(mean_t=mean(t.coupl,na.rm=T)) %>% 
  ggplot(aes(x=class, y=t.coupl,fill=mean_t)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.2,fill="white")+scale_fill_gradientn(colours=colors.coupling.age, limits = c(-5, 5), oob=squish) +
  xlab('')+
  geom_hline(yintercept = 0,linetype='dashed',size=0.5)+
  theme_minimal(base_family = "Gill Sans")+
  theme(axis.text.x = element_text(size = 12,hjust=1, color=c("sienna1","royalblue1","red3","seagreen3"),angle = 90),
        axis.title=element_text(size=12),text = element_text(size=12), legend.position = 'bottom')+
  ylab("Rate of Change in Coupling")+
  labs(fill="Rate of Change in Coupling")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# SI Fig. 19
ggsave(
  ggarrange(bl.coupl.mesulam, sl.coupl.mesulam,
            nrow=1, 
            vjust=0,hjust = 0,
            labels = c('A | Baseline coupling',
                       'B | Rate of change in coupling'),
            font.label = list(size = 12, family="Gill Sans"))+
    theme(plot.margin = margin(0.5,0,0,0, "cm")),
  filename = paste0(outpath.SI,'coupl.mesulam.png'),width=5, height=5, bg='white')

# Between-zone differences in structure-function coupling
comp.mean.bl.coupl = compare_means(bl~class,subset(df.fc.struct.bl.delta, class != 'cortical_wall' & class != 'no label'))[,c(1:5,7)]
comp.mean.bl.coupl[,1]='Baseline Coupling'
comp.mean.sl.coupl = compare_means(sl~class,subset(df.fc.struct.bl.delta, class != 'cortical_wall' & class != 'no label'))[,c(1:5,7)]
comp.mean.sl.coupl[,1]='Rate of Change in Coupling'

latex.out.comp.mean.coupl = rbind(comp.mean.bl.coupl,comp.mean.sl.coupl) %>% kbl(
  format = "latex",
  row.names = F,
  booktabs = T,
  escape = F,
  linesep = c("", "", "","","", "\\midrule")
) %>% collapse_rows(columns = 1, latex_hline = "none")

print(latex.out.comp.mean.coupl, file = paste0(latex.tables,'tab.coupl.mesulam.tex'))

# Spin tests permuations have been pre-calculated
load('data/perm.id.fc.RData')


# Thresholded change in structure-function coupling
# SI Fig. 20 Panel A
sl.thresh.coupling = df.fc.struct.bl.delta %>% 
  ggseg(atlas=glasser,mapping=aes(fill=ifelse(p.fdr<0.05,t.coupl,NA)), position = 'stacked')+
  scale_fill_gradientn(colors=rev(colors.coupling.age), limits=c(-5,5), breaks=c(-5,0,5), oob=squish) +
  theme_minimal(base_family = "Gill Sans")+
  ylab('')+xlab('')+
  labs(fill='Coupling Rate of Change')+
  theme(legend.position = 'bottom')+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Correlation between baseline coupling and baseline MSN degree
# SI Fig. 20 Panel B
p.spin.coupling14.vs.msn14 = perm.sphere.p(bl.fc.str[!is.na(bl.fc.str)],str.l.14[!is.na(bl.fc.str)],perm.id.fc,corr.type = 'spearman')
r.coupling14.vs.msn14 = cor.test(str.l.14,bl.fc.str)$estimate
printout.stats[[14]] = paste0('Correlation baseline coupling with MSN rate of change: ',r.coupling14.vs.msn14,'; p-value = ',p.spin.coupling14.vs.msn14)
coupling14.vs.msn14 = df.fc.struct.bl.delta %>% ggplot(aes(x=msn.bl,y=bl,color=msn.bl))+
  geom_point()+geom_smooth(method='lm', color='black') +
  theme_minimal(base_family = "Gill Sans")+
  xlab('MSN at Baseline')+
  ylab('Coupling at Baseline')+
  scale_color_gradientn(colors=colors.msn.str, limits=c(-0.07,0.07), oob=squish) +
  labs(color='')+
  theme(legend.position = 'bottom')

# Correlation between baseline coupling and rate of change in MSN degree
# SI Fig. 20 Panel C
p.spin.coupling14.vs.msnDelta = perm.sphere.p(bl.fc.str[!is.na(bl.fc.str)],str.l.sl.t[!is.na(bl.fc.str)],perm.id.fc,corr.type = 'spearman')
r.coupling14.vs.msnDelta = cor.test(str.l.sl.t,bl.fc.str)$estimate
printout.stats[[13]] = paste0('Correlation baseline coupling with MSN rate of change: ',r.coupling14.vs.msnDelta,'; p-value = ',p.spin.coupling14.vs.msnDelta)
coupling14.vs.msnDelta = df.fc.struct.bl.delta %>% ggplot(aes(x=msn,y=bl,color=msn))+
  geom_point()+geom_smooth(method='lm', color='black') +
  theme_minimal(base_family = "Gill Sans")+
  xlab('Age effect on MSN')+
  ylab('Coupling at Baseline')+
  scale_color_gradientn(colors=colors.msn.age, limits=c(-5,5), oob=squish) +
  labs(color='')+
  theme(legend.position = 'bottom')

# SI Fig. 20
ggsave(
  ggarrange(sl.thresh.coupling,NULL,coupling14.vs.msn14,NULL,coupling14.vs.msnDelta,
            nrow=1,widths = c(1,0.3,1,0.3,1),
            vjust=0,hjust = 0,
            labels = c('A | Baseline coupling vs baseline MSN',
                       '',
                       'B | Baseline coupling vs MSN age effect',
                       '',
                       'C | Coupling age effect vs MSN age effect'),
            font.label = list(size = 12, family="Gill Sans"))+
    theme(plot.margin = margin(0.5,0,0,0, "cm")),
  filename = paste0(outpath.SI,'coupling.png'),width=10, height=4, bg='white')
