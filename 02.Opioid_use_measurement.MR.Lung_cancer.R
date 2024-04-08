library(TwoSampleMR)

setwd("/home/Method/MR/Opioid_use_measurement.smoking.lung_cancer/")

Opioid_use <- read.csv("./Exp.IV/Opioid_use_measurement.csv",header = T)
# n = 78,808 European (U.K.)
# 22,982 European ancestry cases
# 55,826 European ancestry controls
table(Opioid_use$SNP_Type)
# OUM-all -> Lung cancer #############################################
# Exposure ---------------------------------------------------------------------
  # 相关性设置 P值筛选  
    Exp.Opioid_use <- Opioid_use
    #Exp.Opioid_use <- subset(Opioid_use,Opioid_use$SNP_Type=="OUM&Smoking")    
    #Exp.Opioid_use <- subset(Opioid_use,Opioid_use$SNP_Type=="OUM-Only")
    
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
# Outcome  Lung Cancer----------------------------------------------------------------------
    
    GWAS.ID <- read.table("./Lung_cancer.ieu_ID.txt",sep = "\t")
    
    for (i in 1:length(GWAS.ID$V1)){
      print(GWAS.ID[i,1])
      Out <- extract_outcome_data(snps = Exp.Opioid_use$SNP,
                                  outcomes = GWAS.ID[i,1],
                                  access_token = NULL)
      #Out <- subset(Out,Out$pval.outcome > 0.01)
      exp.out <- harmonise_data(exposure_dat = Exp.Opioid_use,outcome_dat = Out)
      duplicated_rows <- duplicated(exp.out)
      exp.out <- exp.out[!duplicated_rows, ]
    #---------------------------------------------------------------------------
      MR.re = generate_odds_ratios(mr(exp.out))
      mr_scatter_plot(mr(exp.out, method_list =c("mr_ivw",
                                                 "mr_weighted_median",
                                                 "mr_egger_regression") ),exp.out)
      
      MR.re <- MR.re[c(3,2,5,1),-c(1,2,4,7,8,10,11)]
      MR.egger.inter <- mr_pleiotropy_test(exp.out)
      MR.egger.inter <- MR.egger.inter[,-c(1,2,4)]
      MR_heter <- mr_heterogeneity(exp.out)
      MR_heter <- MR_heter[1,-c(1,2,4)]
      
      print(MR.re)
      print(MR.egger.inter)
      print(MR_heter)
    }
    
    
#devtools::install_github("rondolab/MR-PRESSO",force = TRUE)
#library(MRPRESSO)
#mr_presso(data = exp.out, 
#          BetaOutcome = "beta.outcome", 
#          BetaExposure = "beta.exposure", 
#          SdOutcome = "se.outcome", 
#          SdExposure = "se.exposure", 
#          OUTLIERtest = TRUE, 
#          DISTORTIONtest = TRUE, 
#          NbDistribution = 1000,  
#          SignifThreshold = 0.05)

    