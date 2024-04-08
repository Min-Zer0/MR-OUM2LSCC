library(MendelianRandomization)
library(TwoSampleMR)
setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer//")

load("./Exp.list.Rdata")

GWAS.ID <- read.table("./Lung_cancer.ieu_ID.txt",sep = "\t")

MVMR.ivw.data <- data.frame()

for (i in 1:length(GWAS.ID$V1)){
  print(GWAS.ID[i,1])
  Out <- extract_outcome_data(snps =  Exp$SNP,
                              outcomes = GWAS.ID[i,1],
                              access_token = NULL)
  Out <- subset(Out,Out$pval.outcome > 1e-2)
  exp.out <- harmonise_data(exposure_dat =  Exp,outcome_dat = Out)
  duplicated_rows <- duplicated(exp.out)
  exp.out <- exp.out[!duplicated_rows, ]
  colnames(exp.out)
  print(exp.out$outcome[1])
  LungCancer.out <- exp.out[,c("SNP","beta.outcome","se.outcome")]
  MVMR.dat <- merge(Exp.list,LungCancer.out,by="SNP")
  
  nsnp = sum(MVMR.dat$type.Opioid_use)
  print(c("Opioid_use nsnp",nsnp))
  
  nsnp = sum(MVMR.dat$type.CigDay)
  print(c("CigDay nsnp",nsnp))
  
  nsnp = sum(MVMR.dat$type.SmkInit)
  print(c("SmkInit nsnp",nsnp))
  
  nsnp = sum(MVMR.dat$type.SmkCes)
  print(c("SmkCes nsnp",nsnp))
  
  nsnp = sum(MVMR.dat$type.AgeSmk)
  print(c("AgeSmk nsnp",nsnp))
  
  
  MRMVInputObject <- mr_mvinput(bx = cbind(MVMR.dat$beta.Opioid, 
                                           MVMR.dat$beta.CigDay,
                                           MVMR.dat$beta.SmkInit,
                                           MVMR.dat$beta.SmkCes,
                                           MVMR.dat$beta.AgeSmk),
                                bxse = cbind(MVMR.dat$se.Opioid,
                                             MVMR.dat$se.CigDay,
                                             MVMR.dat$se.SmkInit,
                                             MVMR.dat$se.SmkCes,
                                             MVMR.dat$se.AgeSmk),
                                by = MVMR.dat$beta.outcome,
                                byse = MVMR.dat$se.outcome)
  MRMVInputObject
  
  print("IVW")
  ivw.re <- mr_mvivw(MRMVInputObject,
           model = "random",
           robust = T,
           distribution = "normal",
           alpha = 0.05)
  ivw.re.table <- cbind(as.data.frame(ivw.re@Estimate),as.data.frame(ivw.re@StdError))
  rownames(ivw.re.table) =c("Opioid","CigDay","SmkInit","SmkCes","AgeSmk")
  ivw.re.table$OR = exp(ivw.re.table$`ivw.re@Estimate`)
  ivw.re.table$or_lci95 = exp(ivw.re.table$`ivw.re@Estimate` - 1.96*ivw.re.table$`ivw.re@StdError`)
  ivw.re.table$or_uci95 = exp(ivw.re.table$`ivw.re@Estimate` + 1.96*ivw.re.table$`ivw.re@StdError`)
  ivw.re.table$pvalue <- ivw.re@Pvalue
  
  plot.dat <- ivw.re.table[,-c(1:2)]
  plot.dat$exp = paste0(exp.out$outcome[1],rownames(plot.dat))
  MVMR.ivw.data <- rbind(MVMR.ivw.data,plot.dat)
  #MVMR.ivw.data <- rbind(MVMR.ivw.data,Na.dat)
  
  print(ivw.re)
  print(ivw.re.table)
  
  
  print("MR-egger")
  mvegger.re <- mr_mvegger(MRMVInputObject, 
             orientate = 2,
             correl = FALSE,
             distribution = "normal", 
             alpha = 0.05)
  
  mvegger.re.table <- cbind(as.data.frame(mvegger.re@Estimate),as.data.frame(mvegger.re@StdError.Est))
  rownames(mvegger.re.table) =c("Opioid","CigDay","SmkInit","SmkCes","AgeSmk")
  mvegger.re.table$OR = exp(mvegger.re.table$`mvegger.re@Estimate`)
  mvegger.re.table$or_lci95 = exp(mvegger.re.table$`mvegger.re@Estimate` - 1.96*mvegger.re.table$`mvegger.re@StdError.Est`)
  mvegger.re.table$or_uci95 = exp(mvegger.re.table$`mvegger.re@Estimate` + 1.96*mvegger.re.table$`mvegger.re@StdError.Est`)
  mvegger.re.table$pvalue <- mvegger.re@Pvalue.Est
  print(mvegger.re)
  print(mvegger.re.table)
}




p.dat <- MVMR.ivw.data
library(ggplot2)
library(scales)
library(viridis)

p.dat$exp <- factor(p.dat$exp,levels = rev(p.dat$exp))
OR <- p.dat$OR
CI_LL <- p.dat$or_lci95
CI_HL <- p.dat$or_uci95
#P <- 
  ggplot(data=p.dat, aes(x=OR,y=exp, fill=-pvalue))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
 # geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL), color=rep(c("#CFBAF0","#C75146","#AD2E24","#81171B","#540804"),6),height=0.2,size=1.2)+ 
  geom_point(aes(x=OR,y=exp),size=4,shape=23)+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_bw()+
  theme(axis.text = element_text(size = 12,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
P
ggsave(P,filename = "MVMR.png", width = 8, height = 10,dpi = 300)





MVMR.dat.OUM <- subset(MVMR.dat,MVMR.dat$p.Opioid>0.05)

MRMVInputObject <- mr_mvinput(bx = cbind(MVMR.dat.OUM$beta.CigDay,
                                         MVMR.dat.OUM$beta.SmkInit,
                                         MVMR.dat.OUM$beta.SmkCes,
                                         MVMR.dat.OUM$beta.AgeSmk),
                              bxse = cbind(MVMR.dat.OUM$se.CigDay,
                                           MVMR.dat.OUM$se.SmkInit,
                                           MVMR.dat.OUM$se.SmkCes,
                                           MVMR.dat.OUM$se.AgeSmk),
                              by = MVMR.dat.OUM$beta.Opioid,
                              byse = MVMR.dat.OUM$se.Opioid)
MRMVInputObject



ivw.re <- mr_mvivw(MRMVInputObject,
               model = "random",
               robust = T,
               distribution = "normal",
               alpha = 0.05)
ivw.re.table <- cbind(as.data.frame(ivw.re@Estimate),as.data.frame(ivw.re@StdError))
rownames(ivw.re.table) =c("CigDay","SmkInit","SmkCes","AgeSmk")
ivw.re.table$OR = exp(ivw.re.table$`ivw.re@Estimate`)
ivw.re.table$or_lci95 = exp(ivw.re.table$`ivw.re@Estimate` - 1.96*ivw.re.table$`ivw.re@StdError`)
ivw.re.table$or_uci95 = exp(ivw.re.table$`ivw.re@Estimate` + 1.96*ivw.re.table$`ivw.re@StdError`)
ivw.re.table$pvalue <- ivw.re@Pvalue


mvegger.re <- mr_mvegger(MRMVInputObject, 
                 orientate = 2,
                 correl = FALSE,
                 distribution = "normal", 
                 alpha = 0.05)
mvegger.re.table <- cbind(as.data.frame(mvegger.re@Estimate),as.data.frame(mvegger.re@StdError.Est))
rownames(mvegger.re.table) =c("CigDay","SmkInit","SmkCes","AgeSmk")
mvegger.re.table$OR = exp(mvegger.re.table$`mvegger.re@Estimate`)
mvegger.re.table$or_lci95 = exp(mvegger.re.table$`mvegger.re@Estimate` - 1.96*mvegger.re.table$`mvegger.re@StdError.Est`)
mvegger.re.table$or_uci95 = exp(mvegger.re.table$`mvegger.re@Estimate` + 1.96*mvegger.re.table$`mvegger.re@StdError.Est`)
mvegger.re.table$pvalue <- mvegger.re@Pvalue.Est

sum(MVMR.dat.OUM$type.Opioid_use)
sum(MVMR.dat.OUM$type.CigDay)
sum(MVMR.dat.OUM$type.SmkInit)
sum(MVMR.dat.OUM$type.SmkCes)
sum(MVMR.dat.OUM$type.AgeSmk)

ivw.re
ivw.re.table
mvegger.re
mvegger.re.table

p.dat <- ivw.re.table[,-c(1,2)]
p.dat$exp <- rownames(p.dat)
p.dat$exp <- factor(p.dat$exp,levels = rev(p.dat$exp))

library(ggplot2)
library(scales)
library(viridis)

OR <- p.dat$OR
CI_LL <- p.dat$or_lci95
CI_HL <- p.dat$or_uci95
P <- ggplot(data=p.dat, aes(x=OR,y=exp, fill=-pvalue))+ 
  geom_vline(xintercept = 1,linetype='dashed',size=1, color = "grey")+
  geom_errorbarh(aes(xmax=CI_HL, xmin=CI_LL), color=c("#C75146","#AD2E24","#81171B","#540804"),height=0.2,size=1.2)+ 
  geom_point(aes(x=OR,y=exp),size=8,shape=23)+ 
  scale_fill_distiller(palette = "Spectral")+
  theme_bw()+
  theme(axis.text = element_text(size = 12,face = "bold"), 
        axis.title.y = element_blank(),
        legend.position="none",
        title = element_text(size = 16))
P  
ggsave(P,filename = "MVMR.OUM.png", width = 8, height = 3.5,dpi = 300)
