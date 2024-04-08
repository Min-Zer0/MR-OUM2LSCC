library(ggplot2)
setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer/")

#p.dat <- read.csv("./UVMR Prescription Opioid Use (Only)  Tobacco use.csv",header = T)
#p.dat <- read.csv("./UVMR Tobacco use  Prescription Opioid Use.csv",header = T)
p.dat <- read.csv("./MVMR Tobacco use  Prescription Opioid Use.csv",header = T)

p.dat <- na.omit(p.dat)
p.dat$exp <- factor(p.dat$exp,levels = rev(p.dat$exp))
OR <- p.dat$or
CI_LL <- p.dat$or_lci95
CI_HL <- p.dat$or_uci95
P <- 
  ggplot(data=p.dat, aes(x=or,y=exp, fill=-Q_pval))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL), color=c("#C75146","#AD2E24","#81171B","#540804"),height=0.2,size=1.2)+ 
  geom_point(aes(x=OR,y=exp),size=8,shape=23)+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+
  theme(axis.text = element_text(size = 12,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
# P <- 
  ggplot(data=p.dat, aes(x=or,y=exp))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL),height=0.2,size=0.8)+ 
  geom_point(aes(x=OR,y=exp),size=6,shape=22,fill = "#6a8fce")+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+
  theme(axis.text.x = element_text(size = 20,face = "bold"), 
        axis.text.y = element_text(size = 14,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
P  
ggsave(P,filename = "MVMR.Smoking.OUM-only.png", width = 8, height = 3.5,dpi = 300)


p.dat.all <- read.csv("./updata.csv",header = T)
p.dat.all$exp <- paste0(p.dat.all$exp,p.dat.all$out)
p.dat.all$exp <- factor(p.dat.all$exp,levels = rev(p.dat.all$exp))
OR <- p.dat.all$or
CI_LL <- p.dat.all$or_lci95
CI_HL <- p.dat.all$or_uci95
# P <- 
ggplot(data=p.dat.all, aes(x=or,y=exp))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL),height=0.2,size=0.8)+ 
  geom_point(aes(x=OR,y=exp),size=6,shape=22,fill = "#6a8fce")+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+
  theme(axis.text.x = element_text(size = 20,face = "bold"), 
        axis.text.y = element_text(size = 14,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))


p.dat <- p.dat.all[1:6,]
p.dat$exp <- paste0(p.dat$exp,p.dat$out)
p.dat$exp <- factor(p.dat$exp,levels = rev(p.dat$exp))
OR <- p.dat$or
CI_LL <- p.dat$or_lci95
CI_HL <- p.dat$or_uci95
P <- 
  ggplot(data=p.dat, aes(x=or,y=exp, fill=-Q_pval))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL), color=c(rep(c("#ccbcf4"),6)),height=0.2,size=1.2)+ 
  geom_point(aes(x=OR,y=exp),size=8,shape=23)+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+
  theme(axis.text = element_text(size = 12,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
P  
ggsave(P,filename = "OUM.LC.png", width = 10, height = 3.5/4*18,dpi = 300)


p.dat <- p.dat.all[7:18,]
p.dat$exp <- paste0(p.dat$exp,p.dat$out)
p.dat$exp <- factor(p.dat$exp,levels = rev(p.dat$exp))
OR <- p.dat$or
CI_LL <- p.dat$or_lci95
CI_HL <- p.dat$or_uci95
P <- 
  ggplot(data=p.dat, aes(x=or,y=exp, fill=-Q_pval))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL), color=c(rep(c("#a8c4f4"),6),rep(c("#8cecf4"),6)),height=0.2,size=1.2)+ 
  geom_point(aes(x=OR,y=exp),size=8,shape=23)+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_void()+
  theme(axis.text = element_text(size = 12,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
P  
ggsave(P,filename = "OUM.LC.png", width = 10, height = 3.5/4*18,dpi = 300)




p.dat <- read.csv("./MVMR Prescription Opioid Use;  Tobacco use  Lung Cancer.csv",header = T)
  