library(TwoSampleMR)

setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer/")

Opioid_use <- read.table("./Opioid_use_measurement/N02A",header = T)
# n = 78,808 European (U.K.)
# 22,982 European ancestry cases
# 55,826 European ancestry controls

hist(-log10(Opioid_use$p),breaks = 100)

##******************##
Pvalue.t = 5e-5
clump_len = 10000
clump_R2 = 0.01 
##******************##

# Opioid_dependence -> Opioid_use ############################################## 
# Exposure ---------------------------------------------------------------------
# 相关性设置 P值筛选  
Exp.Opioid_use <- subset(Opioid_use,Opioid_use$p < Pvalue.t)
Exp.Opioid_use$Phenotype <- "Opioid_use_measurement"
Exp.Opioid_use$ncase <- 22982
Exp.Opioid_use$ncontrol <- 55826

write.csv(Exp.Opioid_use,file = "./Exp.IV/Opioid_use_measurement.csv")
colnames(Exp.Opioid_use)
# 读取 Exposure
Exp.Opioid_use <- read_exposure_data(filename = "./Exp.IV/Opioid_use_measurement.csv",
                                     clump = FALSE,
                                     sep = ",",
                                     phenotype_col = "Phenotype",
                                     snp_col = "SNP", 
                                     beta_col = "b",
                                     se_col = "se",
                                     eaf_col = "freq",
                                     effect_allele_col = "A1",
                                     other_allele_col = "A2",
                                     pval_col = "p",
                                     #units_col = "units",
                                     ncase_col = "ncase",
                                     ncontrol_col = "ncontrol",
                                     samplesize_col = "n")
# 独立性 设置连锁不平衡
Exp.Opioid_use <- clump_data(Exp.Opioid_use,
                             clump_kb = clump_len,
                             clump_r2 = clump_R2,
                             pop = "EUR")
#######################################################################################################

# 统计强度设置
CigDay <- read.table("./Smoking/data/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkInit<- read.table("./Smoking/data/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkCes<- read.table("./Smoking/data/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
AgeSmk<- read.table("./Smoking/data/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)

Exp.Opioid_use.output <- Exp.Opioid_use[,-c(9,10,12:15)]
Exp.Opioid_use.output <- merge(Exp.Opioid_use.output,CigDay,by.x="SNP",by.y="RSID")
Exp.Opioid_use.output <- Exp.Opioid_use.output[,c(10,11,1:9)]
colnames(Exp.Opioid_use.output) = c("Chr","Position","rsID","Effect_Allele","Other_Allele",
                                    "EAF","Beta","SE","P","SampleSize","Exposure")
Exp.Opioid_use.output$R2 <- (Exp.Opioid_use.output$Beta**2)*(1-Exp.Opioid_use.output$EAF)*2*Exp.Opioid_use.output$EAF
Exp.Opioid_use.output$F_statistic <-Exp.Opioid_use.output$R2*(Exp.Opioid_use.output$SampleSize-2)/(1-Exp.Opioid_use.output$R2)
head(Exp.Opioid_use.output)
min(Exp.Opioid_use.output$F_statistic)

  CigDay <- CigDay[,c(3,9)]
  colnames(CigDay) =c("rsID","P.CigDay")
  head(CigDay)
  
  SmkInit <- SmkInit[,c(3,9)]
  colnames(SmkInit) =c("rsID","P.SmkInit")
  head(SmkInit)
  
  SmkCes <- SmkCes[,c(3,9)]
  colnames(SmkCes) =c("rsID","P.SmkCes")
  head(SmkCes)
  
  AgeSmk <- AgeSmk[,c(3,9)]
  colnames(AgeSmk) =c("rsID","P.AgeSmk")
  head(AgeSmk)

  Exp.Opioid_use.output <- merge(Exp.Opioid_use.output,CigDay,by="rsID")
  Exp.Opioid_use.output <- merge(Exp.Opioid_use.output,SmkInit,by="rsID")
  Exp.Opioid_use.output <- merge(Exp.Opioid_use.output,SmkCes,by="rsID")
  Exp.Opioid_use.output <- merge(Exp.Opioid_use.output,AgeSmk,by="rsID")
  
  Exp.Opioid_use.output$SNP_Type <-ifelse(Exp.Opioid_use.output$P.CigDay < 0.05 |
                                          Exp.Opioid_use.output$P.SmkInit < 0.05 |
                                          Exp.Opioid_use.output$P.SmkCes < 0.05 |
                                          Exp.Opioid_use.output$P.AgeSmk < 0.05,    
                                          "OUM&Smoking","OUM-Only")
  
  
  table(Exp.Opioid_use.output$SNP_Type)
  Exp.Opioid_use.output <- Exp.Opioid_use.output[,-c(14:17)]
  head(Exp.Opioid_use.output)
  write.csv(Exp.Opioid_use.output,file = "./Exp.IV/Opioid_use_measurement.csv")
  
    