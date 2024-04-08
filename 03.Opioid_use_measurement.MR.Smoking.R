library(TwoSampleMR)

setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer/")

Opioid_use <- read.csv("./Exp.IV/Opioid_use_measurement.csv",header = T)
table(Opioid_use$SNP_Type)
# n = 78,808 European (U.K.)
# 22,982 European ancestry cases
# 55,826 European ancestry controls

CigDay <- read.table("./Smoking/data/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkInit<- read.table("./Smoking/data/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkCes<- read.table("./Smoking/data/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
AgeSmk<- read.table("./Smoking/data/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)

somking_list = list(CigDay, SmkInit,SmkCes,AgeSmk)
somking_list
names(somking_list) = c('CigDay', 'SmkInit','SmkCes','AgeSmk')


# OUM-all -> Lung cancer #############################################
# Exposure ---------------------------------------------------------------------
# 相关性设置 P值筛选  
#Exp.Opioid_use <- Opioid_use
#Exp.Opioid_use <- subset(Opioid_use,Opioid_use$SNP_Type=="OUM&Smoking")    
Exp.Opioid_use <- subset(Opioid_use,Opioid_use$SNP_Type=="OUM-Only")

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

# Outcome Smoking ----------------------------------------------------------------------


Exp <- Exp.Opioid_use
colnames(Exp) <- paste0("exp.",colnames(Exp))

for (i in 1:length(somking_list)){
  #i = 2
  Out.smoking <- merge(Exp,somking_list[[i]],by.x="exp.SNP",by.y="RSID")
  Out.smoking <- Out.smoking[,c(1,16:24)]
  colnames(Out.smoking)[1] <- "SNP"
  Out.smoking$Phenotype = names(somking_list[i])
  head(Out.smoking)

  write.csv(Out.smoking,file = "./Out.somking.csv")
  colnames(Out.smoking)
  
  Out.smoking <- read_outcome_data(filename = "./Out.somking.csv",
                                   sep = ",",
                                   phenotype_col = "Phenotype",
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "AF_1000G",
                                   effect_allele_col = "EFFECT_ALLELE",
                                   other_allele_col = "OTHER_ALLELE",
                                   pval_col = "P",
                                   samplesize_col = "N",
                                   chr_col = "CHR",
                                   pos_col = "POS")
  
    exp.out <- harmonise_data(exposure_dat = Exp.Opioid_use,
                              outcome_dat = Out.smoking)
  

  MR.re = generate_odds_ratios(mr(exp.out))
  MR.re <- MR.re[c(3,2,5,1),-c(1,2,4,7,8,10,11)]
  MR.egger.inter <- mr_pleiotropy_test(exp.out)
  MR.egger.inter <- MR.egger.inter[,-c(1,2,4)]
  MR_heter <- mr_heterogeneity(exp.out)
  MR_heter <- MR_heter[1,-c(1,2,4)]
  print(MR.re)
  print(MR.egger.inter)
  print(MR_heter)
}
