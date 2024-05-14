#### Load required packages
library(ggplot2)

#### Load dependencies
source('./scripts/external/rotate_parcellation-master/R/perm.sphere.p.R')

#### Load data
nm.NSPN = read.csv('results/main/NSPN.nm.txt', row.names = 1)[,1]                                      # NSPN regions names for matching
t.NSPN = read.csv('results/main/str.l.sl.t.csv',header=F)[,1]                             # NSPN rate of change in MSN degree
t.morph.NSPN = read.csv('results/main/str.morph.l.sl.t.csv',header=T)                     # NSPN rate of change mophometric feaures

nm.HCPD = read.table('results/SI/sensitivity-analyses/out-of-sample-replication/main/HCPD.nm.txt', header=T)[,1]             # HCP-D regions names for matching
t.HCPD = read.csv('results/SI/sensitivity-analyses/out-of-sample-replication/main/str.l.sl.t.csv',header = F)[,1]      # HCP-D rate of change in MSN degree
t.morph.HCPD = read.csv('results/SI/sensitivity-analyses/out-of-sample-replication/main/str.morph.l.sl.t.csv',header=T)    # HCP-D rate of change mophometric feaures

reorder.idx = match(nm.NSPN,nm.HCPD) # check this match(nm.HCPD[match(nm.NSPN,nm.HCPD)],nm.NSPN)
load('data/perm.id.RData')

# Estimate correlation between original replication rates of change in individual morphometric features
sm = c("FA", "MD", "MT", "SA", "GM", "CT")
nmorph=6
CORR=SPIN=vector(length=nmorph)
P=list()
for (i in 1:nmorph) {
  CORR[i] = cor.test(t(t.morph.NSPN[i,]),t(t.morph.HCPD[i,reorder.idx]))$estimate
  SPIN[i] = perm.sphere.p(as.vector(t(t.morph.NSPN[i,])),as.vector(t(t.morph.HCPD[i,reorder.idx])),perm.id = perm.id)
  P[[i]] = ggplot(data.frame(NSPN = as.vector(t(t.morph.NSPN[i,])), 
                             HCPD= as.vector(t(t.morph.HCPD[i,reorder.idx]))),
                  aes(x=NSPN,y=HCPD))+
    geom_point()+
    geom_smooth(method='lm')+
    theme_bw()
}

# Estimate correlation between original replication rate of change in MSN degree
df = data.frame(NSPN = t.NSPN,HCPD=t.HCPD[reorder.idx])
cor.test(t.NSPN,t.HCPD[reorder.idx])$estimate
p.spin.msn = perm.sphere.p(t.NSPN,t.HCPD[reorder.idx], perm.id)

panel.correlation = df %>% ggplot(aes(x=NSPN,y=HCPD))+
  geom_point()+
  theme_minimal(base_family = 'Gill Sans')+
  geom_smooth(method='lm')

# SI Fig 17, Panel D
ggsave(plot=panel.correlation, filename=paste0('results/SI/sensitivity-analyses/out-of-sample-replication/correlation.hcpd.nspn.png'),height=4,width=4)


