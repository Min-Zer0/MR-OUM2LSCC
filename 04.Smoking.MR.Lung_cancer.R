library(TwoSampleMR)

setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer//")

CigDay <- read.table("./Smoking/data/GSCAN_CigDay_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkInit<- read.table("./Smoking/data/GSCAN_SmkInit_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
SmkCes<- read.table("./Smoking/data/GSCAN_SmkCes_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)
AgeSmk<- read.table("./Smoking/data/GSCAN_AgeSmk_2022_GWAS_SUMMARY_STATS_EUR.txt",header = T)

##******************##
  Pvalue.t = 5e-8
  MAF.t = 0.01
  clump_len = 10000
  clump_R2 = 0.001 
##******************##
  
# Somking -> Lung cancer #######################################################
# Exposure ---------------------------------------------------------------------
  # 相关性设置 P值筛选

################################################################################  
  #CigDay$MAF=ifelse(CigDay$AF_1000G < 0.5, CigDay$AF_1000G ,1- CigDay$AF_1000G)
  #head(CigDay)
  #Exp.somking <- subset(CigDay,CigDay$P < Pvalue.t)
  #Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
  #Exp.somking$Phenotype <- "CigDay"
  
  #SmkInit$MAF=ifelse(SmkInit$AF_1000G < 0.5, SmkInit$AF_1000G ,1- SmkInit$AF_1000G)
  #head(SmkInit)
  #Exp.somking <- subset(SmkInit,SmkInit$P < Pvalue.t)
  #Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
  #Exp.somking$Phenotype <- "SmkInit"
  
  #SmkCes$MAF=ifelse(SmkCes$AF_1000G < 0.5, SmkCes$AF_1000G ,1- SmkCes$AF_1000G)
  #head(SmkCes)
  #Exp.somking <- subset(SmkCes,SmkCes$P < Pvalue.t)
  #Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
  #Exp.somking$Phenotype <- "SmkCes"
  
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
  
  Exp.somking.output <- Exp.somking[,c(3,1,2,4:11)]
  colnames(Exp.somking.output) = c("rsID","Chr","Position","Effect_Allele","Other_Allele",
                                      "EAF","Beta","SE","P","SampleSize","Exposure")
  Exp.somking.output$R2 <- (Exp.somking.output$Beta**2)*(1-Exp.somking.output$EAF)*2*Exp.somking.output$EAF
  Exp.somking.output$F_statistic <-Exp.somking.output$R2*(Exp.somking.output$SampleSize-2)/(1-Exp.somking.output$R2)
  head(Exp.somking.output)
  min(Exp.somking.output$F_statistic)
  Exp.somking.output$SNP_Type <- "-"
  
  
  write.csv(Exp.somking.output, file = "./Exp.IV/AgeSmk.csv")
  # 统计强度设置
    
# Outcome ----------------------------------------------------------------------
  GWAS.ID <- read.table("./Lung_cancer.ieu_ID.txt",sep = "\t")
  
  for (i in 1:length(GWAS.ID$V1)){
    #i = 3
    print(GWAS.ID[i,1])
    Out <- extract_outcome_data(snps =  Exp.somking$SNP,
                                outcomes = GWAS.ID[i,1],
                                access_token = NULL)
    Out <- subset(Out,Out$pval.outcome > 1e-2)
    
    exp.out <- harmonise_data(exposure_dat =  Exp.somking,outcome_dat = Out)
    duplicated_rows <- duplicated(exp.out)
    exp.out <- exp.out[!duplicated_rows, ]
    MR.re = generate_odds_ratios(mr(exp.out))
    MR.re <- as.data.frame(MR.re)
    MR.re <- MR.re[c(3,2,5,1),-c(1,2,4,7,8,10,11)]
    MR.egger.inter <- mr_pleiotropy_test(exp.out)
    MR.egger.inter <- MR.egger.inter[,-c(1,2,4)]
    MR_heter <- mr_heterogeneity(exp.out)
    MR_heter <- MR_heter[1,-c(1,2,4)]
    print(MR.re)
    print(MR.egger.inter)
    print(MR_heter)
  }
    