library(TwoSampleMR)

setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer//")

Opioid_use <- read.table("./Opioid_use_measurement/N02A",header = T)

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


##******************##
  Pvalue.t = 5e-8
  MAF.t = 0.01
  clump_len = 10000
  clump_R2 = 0.001 
##******************##

  
  
  for (i in 1:length(somking_list)){
    print(names(somking_list[i]))
    Exp.somking <- somking_list[[i]]
    Exp.somking$MAF=ifelse(Exp.somking$AF_1000G < 0.5, Exp.somking$AF_1000G ,1- Exp.somking$AF_1000G)
    Exp.somking <- subset(Exp.somking,Exp.somking$P < Pvalue.t)
    Exp.somking <- subset(Exp.somking,Exp.somking$MAF > 0.01)
    Exp.somking$Phenotype <- names(somking_list[i])
  
  write.csv(Exp.somking,file = "./Exp.somking.csv")
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
  
  Out.Opioid_use <- merge(Exp.somking, Opioid_use,by="SNP")
  Out.Opioid_use <- Out.Opioid_use[,c(1,16:22)]
  Out.Opioid_use <- subset(Out.Opioid_use, Out.Opioid_use$p > 0.05)
  Out.Opioid_use$Phenotype = "Opioid_use_measurement"
  
  write.csv(Out.Opioid_use,file = "./Out.Opioid.csv")
  colnames(Out.Opioid_use)
  
  Out.Opioid <- read_outcome_data(filename = "./Out.Opioid.csv",
                                   sep = ",",
                                   phenotype_col = "Phenotype",
                                   snp_col = "SNP",
                                   beta_col = "b",
                                   se_col = "se",
                                   eaf_col = "freq",
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2",
                                   pval_col = "p",
                                   samplesize_col = "n")
  
  exp.out <- harmonise_data(exposure_dat = Exp.somking,
                            outcome_dat = Out.Opioid)
  
  
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
    