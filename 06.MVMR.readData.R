library(TwoSampleMR)
setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer/")


Opioid_use <- read.csv("./Exp.IV/Opioid_use_measurement.csv",header = T)
table(Opioid_use$SNP_Type)
# n = 78,808 European (U.K.)
# 22,982 European ancestry cases
# 55,826 European ancestry controls

Exp.Opioid_use <- Opioid_use
Exp.Opioid_use$ncase <- 22982
Exp.Opioid_use$ncontrol <- 55826

head(Exp.Opioid_use)
write.csv(Exp.Opioid_use,file = "./Exp.csv")
colnames(Exp.Opioid_use)
# 读取 Exposure
Exp.Opioid_use <- read_exposure_data(filename = "./Exp.csv",
                                     clump = FALSE,
                                     sep = ",",
                                     phenotype_col = "Exposure",
                                     snp_col = "rsID",
                                     beta_col = "Beta",
                                     se_col = "SE",
                                     eaf_col = "EAF",
                                     effect_allele_col = "Effect_Allele",
                                     other_allele_col = "Other_Allele",
                                     pval_col = "P",
                                     ncase_col = "ncase",
                                     ncontrol_col = "ncontrol",
                                     samplesize_col = "SampleSize")



##******************##
Pvalue.t = 5e-8
MAF.t = 0.01
clump_len = 10000
clump_R2 = 0.001 
##******************##


CigDay <- read.table("./Smoking/data/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
CigDay$MAF=ifelse(CigDay$AF_1000G < 0.5, CigDay$AF_1000G ,1- CigDay$AF_1000G)
head(CigDay)
Exp.somking <- subset(CigDay,CigDay$P < Pvalue.t)
Exp.somking <- subset(Exp.somking,Exp.somking$MAF > MAF.t)
Exp.somking$Phenotype <- "CigDay"

write.csv(Exp.somking,file = "./Exp.somking.csv")
colnames(Exp.somking)
# 读取 Exposure
Exp.somking <- read_exposure_data(filename = "./Exp.somking.csv",
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "Phenotype",
                                  snp_col = "RSID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "AF_1000G",
                                  effect_allele_col = "EFFECT_ALLELE",
                                  other_allele_col = "OTHER_ALLELE",
                                  pval_col = "P",
                                  samplesize_col = "N",
                                  chr_col = "CHR",
                                  pos_col = "POS")
# 独立性 设置连锁不平衡
Exp.somking <- clump_data(Exp.somking,
                          clump_kb = clump_len,
                          clump_r2 = clump_R2,
                          pop = "EUR")
Exp.CigDay <-Exp.somking


SmkInit<- read.table("./Smoking/data/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkInit$MAF=ifelse(SmkInit$AF_1000G < 0.5, SmkInit$AF_1000G ,1- SmkInit$AF_1000G)
head(SmkInit)
Exp.somking <- subset(SmkInit,SmkInit$P < Pvalue.t)
Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
Exp.somking$Phenotype <- "SmkInit"

write.csv(Exp.somking,file = "./Exp.somking.csv")
colnames(Exp.somking)
# 读取 Exposure
Exp.somking <- read_exposure_data(filename = "./Exp.somking.csv",
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "Phenotype",
                                  snp_col = "RSID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "AF_1000G",
                                  effect_allele_col = "EFFECT_ALLELE",
                                  other_allele_col = "OTHER_ALLELE",
                                  pval_col = "P",
                                  samplesize_col = "N",
                                  chr_col = "CHR",
                                  pos_col = "POS")
# 独立性 设置连锁不平衡
Exp.somking <- clump_data(Exp.somking,
                          clump_kb = clump_len,
                          clump_r2 = clump_R2,
                          pop = "EUR")
Exp.SmkInit <-Exp.somking




SmkCes<- read.table("./Smoking/data/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkCes$MAF=ifelse(SmkCes$AF_1000G < 0.5, SmkCes$AF_1000G ,1- SmkCes$AF_1000G)
head(SmkCes)
Exp.somking <- subset(SmkCes,SmkCes$P < Pvalue.t)
Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
Exp.somking$Phenotype <- "SmkCes"

write.csv(Exp.somking,file = "./Exp.somking.csv")
colnames(Exp.somking)
# 读取 Exposure
Exp.somking <- read_exposure_data(filename = "./Exp.somking.csv",
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "Phenotype",
                                  snp_col = "RSID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "AF_1000G",
                                  effect_allele_col = "EFFECT_ALLELE",
                                  other_allele_col = "OTHER_ALLELE",
                                  pval_col = "P",
                                  samplesize_col = "N",
                                  chr_col = "CHR",
                                  pos_col = "POS")
# 独立性 设置连锁不平衡
Exp.somking <- clump_data(Exp.somking,
                          clump_kb = clump_len,
                          clump_r2 = clump_R2,
                          pop = "EUR")
Exp.SmkCes <-Exp.somking



AgeSmk<- read.table("./Smoking/data/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
AgeSmk$MAF=ifelse(AgeSmk$AF_1000G < 0.5, AgeSmk$AF_1000G ,1- AgeSmk$AF_1000G)
head(AgeSmk)
Exp.somking <- subset(AgeSmk,AgeSmk$P < Pvalue.t)
Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
Exp.somking$Phenotype <- "AgeSmk"

write.csv(Exp.somking,file = "./Exp.somking.csv")
colnames(Exp.somking)
# 读取 Exposure
Exp.somking <- read_exposure_data(filename = "./Exp.somking.csv",
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "Phenotype",
                                  snp_col = "RSID",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  eaf_col = "AF_1000G",
                                  effect_allele_col = "EFFECT_ALLELE",
                                  other_allele_col = "OTHER_ALLELE",
                                  pval_col = "P",
                                  samplesize_col = "N",
                                  chr_col = "CHR",
                                  pos_col = "POS")
# 独立性 设置连锁不平衡
Exp.somking <- clump_data(Exp.somking,
                          clump_kb = clump_len,
                          clump_r2 = clump_R2,
                          pop = "EUR")

Exp.AgeSmk <-Exp.somking

save(Exp.CigDay,Exp.SmkInit,Exp.SmkCes,Exp.AgeSmk,Exp.Opioid_use,file = "MVMR.Exp.Rdata")


################################################################################
load("./MVMR.Exp.Rdata")
Exp.Opioid_use$type.Opioid_use = 1
Exp.CigDay$type.CigDay = 1
Exp.SmkInit$type.SmkInit = 1
Exp.SmkCes$type.SmkCes = 1
Exp.AgeSmk$type.AgeSmk = 1


Exp.IV <- merge(Exp.Opioid_use,Exp.CigDay,by="SNP",all =T)
Exp.IV <- merge(Exp.IV,Exp.SmkInit,by="SNP",all =T)
Exp.IV <- merge(Exp.IV,Exp.SmkCes,by="SNP",all =T)
Exp.IV <- merge(Exp.IV,Exp.AgeSmk,by="SNP",all =T)
colnames(Exp.IV)
Exp.IV <- Exp.IV[,c("SNP","type.Opioid_use","type.CigDay","type.SmkInit","type.SmkCes","type.AgeSmk")]
Exp.IV[is.na(Exp.IV)] <- 0


Opioid_use <- read.table("./Opioid_use_measurement/N02A",header = T)
CigDay <- read.table("./Smoking/data/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkInit<- read.table("./Smoking/data/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkCes<- read.table("./Smoking/data/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
AgeSmk<- read.table("./Smoking/data/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)


Opioid_use <- subset(Opioid_use, SNP %in% Exp.IV$SNP)
SNP.ID <- Opioid_use$SNP
CigDay  <- subset(CigDay, RSID %in% SNP.ID)
SmkInit  <- subset(SmkInit, RSID %in% SNP.ID)
SmkCes  <- subset(SmkCes, RSID %in% SNP.ID)
AgeSmk  <- subset(AgeSmk, RSID %in% SNP.ID)

save(CigDay,SmkInit,SmkCes,AgeSmk,Opioid_use,file = "MVMR.Exp.table.Rdata")

load("MVMR.Exp.table.Rdata")
Opioid_use <- Opioid_use[,c(1:7)]
CigDay <- CigDay[,c(3:9)]
SmkInit <- SmkInit[,c(3:9)]
SmkCes <- SmkCes[,c(3:9)]
AgeSmk <- AgeSmk[,c(3:9)]

colnames(Opioid_use)=c("SNP","effect_allele","other_allele","eaf","beta","se","p")

colnames(CigDay)=c("SNP","effect_allele","other_allele","eaf","beta","se","p")
colnames(SmkInit)=c("SNP","effect_allele","other_allele","eaf","beta","se","p")
colnames(SmkCes)=c("SNP","effect_allele","other_allele","eaf","beta","se","p")
colnames(AgeSmk)=c("SNP","effect_allele","other_allele","eaf","beta","se","p")

write.csv(Opioid_use,file = "Exp.csv")
Exp <- read_exposure_data(filename = "./Exp.csv",clump = FALSE,sep = ",",snp_col = "SNP",
                          beta_col = "beta",se_col = "se",eaf_col = "eaf",pval_col = "p",
                          effect_allele_col = "effect_allele",other_allele_col = "other_allele")


##=============================================================================#
write.csv(CigDay,file = "Out.csv")
Out <- read_outcome_data(filename = "./Out.csv",sep = ",",
                         snp_col = "SNP",beta_col = "beta",
                         se_col = "se",eaf_col = "eaf",pval_col = "p",
                         effect_allele_col = "effect_allele",other_allele_col = "other_allele")

exp.out <- harmonise_data(exposure_dat = Exp,
                          outcome_dat = Out)
colnames(exp.out)
Out <- exp.out[,c("SNP","beta.outcome","se.outcome","pval.outcome")]
colnames(Out) <- c("SNP","beta","se","p")
colnames(Out)[2:4] <-paste0(colnames(Out)[2:4],".CigDay") 

CigDay <- Out
head(CigDay)
##=============================================================================#
write.csv(SmkInit,file = "Out.csv")
Out <- read_outcome_data(filename = "./Out.csv",sep = ",",
                         snp_col = "SNP",beta_col = "beta",
                         se_col = "se",eaf_col = "eaf",pval_col = "p",
                         effect_allele_col = "effect_allele",other_allele_col = "other_allele")

exp.out <- harmonise_data(exposure_dat = Exp,
                          outcome_dat = Out)
Out <- exp.out[,c("SNP","beta.outcome","se.outcome","pval.outcome")]
colnames(Out) <- c("SNP","beta","se","p")
colnames(Out)[2:4] <-paste0(colnames(Out)[2:4],".SmkInit") 
SmkInit <- Out

##=============================================================================#
write.csv(SmkCes,file = "Out.csv")
Out <- read_outcome_data(filename = "./Out.csv",sep = ",",
                         snp_col = "SNP",beta_col = "beta",
                         se_col = "se",eaf_col = "eaf",pval_col = "p",
                         effect_allele_col = "effect_allele",other_allele_col = "other_allele")

exp.out <- harmonise_data(exposure_dat = Exp,
                          outcome_dat = Out)
Out <- exp.out[,c("SNP","beta.outcome","se.outcome","pval.outcome")]
colnames(Out) <- c("SNP","beta","se","p")
colnames(Out)[2:4] <-paste0(colnames(Out)[2:4],".SmkCes") 
SmkCes <- Out

##=============================================================================#
write.csv(AgeSmk,file = "Out.csv")
Out <- read_outcome_data(filename = "./Out.csv",sep = ",",
                         snp_col = "SNP",beta_col = "beta",
                         se_col = "se",eaf_col = "eaf",pval_col = "p",
                         effect_allele_col = "effect_allele",other_allele_col = "other_allele")

exp.out <- harmonise_data(exposure_dat = Exp,
                          outcome_dat = Out)
Out <- exp.out[,c("SNP","beta.outcome","se.outcome","pval.outcome")]
colnames(Out) <- c("SNP","beta","se","p")
colnames(Out)[2:4] <-paste0(colnames(Out)[2:4],".AgeSmk") 
AgeSmk <- Out


head(Opioid_use)
Opioid_use <- Opioid_use[,c(1,5,6,7)]
colnames(Opioid_use) <- c("SNP","beta.Opioid","se.Opioid","p.Opioid")


Exp.list <- merge(Opioid_use,CigDay,by="SNP")
Exp.list <- merge(Exp.list,SmkInit,by="SNP")
Exp.list <- merge(Exp.list,SmkCes,by="SNP")
Exp.list <- merge(Exp.list,AgeSmk,by="SNP")

Exp.list <- merge(Exp.list,Exp.IV,by="SNP")

head(Exp.list)

save(Exp.list,Exp,file = "./Exp.list.Rdata")



load("./Exp.list.Rdata")

write.csv(Exp.list,file = "./Exp.list.csv")
