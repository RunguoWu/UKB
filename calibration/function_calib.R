library(tidyverse)
library(survival)
library(stats)
library(eha)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

source(file.path(prep, "model_preparation.R"))

# cf <- readRDS(file.path(vali_data, "cf_20210508.rds"))

cf <- readRDS(file.path(vali_data, "cf_20210714.rds"))

'%ni%' <- Negate('%in%')


# prepare data ------------------------------------------------------------

# to be cautious, only generate those with all 0
for (i in 1:26) {
  text <- paste0("st",i)
  d4s[[text]] <- rep(0,dim(d4s)[1])
}
for (i in c("wstatina", "wstatinayr1", "wstatinayr2And", "ind6_1", "ind7_1")) {
  d4s[[i]] <- rep(0,dim(d4s)[1])
}

# generate an indicator of old diabetes RW: 21-09-20
d4s$diabetes_baseline_0_10 <- 0

d4s$diabetes_baseline_0_10[d4s$diabetes.fo.pre==1 & 
                     d4s$recruit.date - d4s$diabetes.fo.pre.date<9*365.25] <- 1

d4s$diabetes_baseline_10_20 <- 0

d4s$diabetes_baseline_10_20[d4s$diabetes.fo.pre==1 & 
                        d4s$recruit.date - d4s$diabetes.fo.pre.date>=9*365.25 &
                          d4s$recruit.date - d4s$diabetes.fo.pre.date < 19*365.25] <- 1

d4s$diabetes_baseline_20_inf <- 0
d4s$diabetes_baseline_20_inf[d4s$diabetes.fo.pre==1 & 
                              d4s$recruit.date - d4s$diabetes.fo.pre.date>=19*365.25] <- 1

d4s$diabetes_baseline_10yr <- 0
d4s$diabetes_baseline_10yr[d4s$diabetes_baseline_0_10==1] <- 1
d4s$diabetes_baseline_10yr[d4s$diabetes_baseline_10_20==1] <- 2
d4s$diabetes_baseline_10yr[d4s$diabetes_baseline_20_inf==1] <- 3


# calibration function ----------------------------------------------------

# e belongs to c("MI", "stroke", "CRV", "cancer","VD", "NVD", "diabetes")
# p belongs to c("pp", "sp")
# d belongs to c("exp", "wei, "gom")
# gp_record: TRUE indicates include people with GP record only
# var: cancer familiy, dm family, smoke, PA, mental, townsend
# var must be categorical variables
calib <- function(e, p, d, var=NULL, gp_record=FALSE, rmcvd = FALSE){
  
  # generate data for modeling from d4s ----
  if (e %ni% c("cancer", "diabetes")){
    if (p=="pp") d4m <- subset(d4s, CVhist==0) else d4m <- subset(d4s, CVhist>0)
  } else {
    if (e=="cancer") d4m <- subset(d4s, cancer.baseline.all==0) 
    else if (e=="diabetes") d4m <- subset(d4s, gp_record==1 & diabetes.fo.pre==0)
  }
  
  # activate this below when estimating diabetes coefficients
  # mute them if diabetes covars are not the key.
  # if (!is.null(var) & sum(grepl("dm", var))>0) d4m <- subset(d4m, gp_record==1) # keep those with GP records only for estimating diabetes effects
  
  if(gp_record) d4m <- subset(d4m, gp_record==1)
  
  # ### test for people with baseline dm only
  # # mute it after use
  # d4m <- subset(d4m, diabetes.fo.pre==1)
  # ###
  
  
  surv <- as.formula(paste0("Surv(", e, ".time", ", ", e, ".status", ")~."))
  
  dt <- survSplit(surv, d4m, 
                  cut = 365.25*(1:15), episode = "timegroup")
  
  dt$cir.start <- dt$recruit.date + dt$tstart
  dt$cir.stop <- dt$recruit.date + dt[[paste0(e, ".time")]]
  
  # age
  # first timegroup since recruitment is 1 
  # per 10 years age, centred at 60
  dt$CurrAge_cent <- (dt$age.recruit - 60 + dt$timegroup - 1)/10 
  
  # generate covariate events by year ---- 
  events <- c("MI", "stroke", "CRV") # generate cancer and dm covar separately
  
  for (cov in setdiff(events, e)) {
    dt[[paste0(tolower(cov), "_0_1")]] <- 
      dt[[paste0(tolower(cov), "_1_2")]] <- 
      dt[[paste0(tolower(cov), "_2_3")]] <- 
      dt[[paste0(tolower(cov), "_3_inf")]] <- rep(0, nrow(dt))
    
    dt[[paste0(tolower(cov), "_0_1")]][dt[[paste0(cov, ".b4", e)]]==1 & 
           dt[[paste0(cov, ".b4", e, ".date")]] > dt$cir.start &
           dt[[paste0(cov, ".b4", e, ".date")]] < dt$cir.start + 365.25] <- 1 
    
    dt[[paste0(tolower(cov), "_1_2")]][dt[[paste0(cov, ".b4", e)]]==1 & 
           dt[[paste0(cov, ".b4", e, ".date")]] + 365.25 > dt$cir.start &
           dt[[paste0(cov, ".b4", e, ".date")]] + 365.25 <= dt$cir.start + 365.25] <- 1                      
    
    dt[[paste0(tolower(cov), "_2_3")]][dt[[paste0(cov, ".b4", e)]]==1 & 
           dt[[paste0(cov, ".b4", e, ".date")]] + 730.5 > dt$cir.start &
           dt[[paste0(cov, ".b4", e, ".date")]] + 730.5 <= dt$cir.start + 365.25] <- 1
    
    dt[[paste0(tolower(cov), "_3_inf")]][dt[[paste0(cov, ".b4", e)]]==1 & 
           dt[[paste0(cov, ".b4", e, ".date")]] + 1095.75 <= dt$cir.start + 365.25] <- 1
    
    # for level collapse
    dt[[paste0(tolower(cov), "_2_inf")]] <- dt[[paste0(tolower(cov), "_1_inf")]] <- 
      dt[[paste0(tolower(cov), "_0_inf")]] <- dt[[paste0(tolower(cov), "_none")]] <- 0
    
    dt[[paste0(tolower(cov), "_2_inf")]][dt[[paste0(tolower(cov), "_3_inf")]]==1 | 
                                           dt[[paste0(tolower(cov), "_2_3")]]==1] <- 1
    
    dt[[paste0(tolower(cov), "_1_inf")]][dt[[paste0(tolower(cov), "_2_inf")]]==1 | 
                                           dt[[paste0(tolower(cov), "_1_2")]]==1] <- 1  
    
    dt[[paste0(tolower(cov), "_0_inf")]][dt[[paste0(tolower(cov), "_1_inf")]]==1 | 
                                           dt[[paste0(tolower(cov), "_0_1")]]==1] <- 1
    
    dt[[paste0(tolower(cov), "_none")]] <- ifelse(dt[[paste0(tolower(cov), "_0_inf")]]==1, 0, 1)
  }

  # generate diabetes by 10 year interval ----
  if (e != "diabetes") {
    dt[["dm_0_10"]] <- dt[["dm_10_inf"]] <- rep(0, nrow(dt))
    
    dt[["dm_0_10"]][dt[[paste0("diabetes.b4", e)]]==1 &
                      dt[[paste0("diabetes.b4", e, ".date")]] + 3287.25 > dt$cir.start &
                      dt[[paste0("diabetes.b4", e, ".date")]] < dt$cir.start + 365.25] <- 1
    
    dt[["dm_10_inf"]][dt[[paste0("diabetes.b4", e)]]==1 &
                        dt[[paste0("diabetes.b4", e, ".date")]] + 3652.5 <= dt$cir.start + 365.25] <- 1
    
    
    dt[["dm_10_inf_old"]] <- dt[["dm_10_inf_new"]] <- rep(0, nrow(dt))
    dt[["dm_10_inf_old"]][dt$diabetes_VeryOld==1] <- 1
    dt[["dm_10_inf_new"]][dt[["dm_10_inf"]]==1 & dt$diabetes_VeryOld==0] <- 1
    
    
    # code pre-diabetes
    
    dt[["dmPre_1"]] <- dt[["dmPre_2"]] <- dt[["dmPre_3"]] <- dt[["dmPre_4"]] <- 0
    
    dt[["dmPre_1"]][dt[["dm_0_10"]]==0 & dt[["dm_10_inf"]]==0 & 
                      dt[["hba1c"]]<32] <- 1
    
    dt[["dmPre_2"]][dt[["dm_0_10"]]==0 & dt[["dm_10_inf"]]==0 & 
                      dt[["hba1c"]]>=32 & dt[["hba1c"]]<37] <- 1
    
    dt[["dmPre_3"]][dt[["dm_0_10"]]==0 & dt[["dm_10_inf"]]==0 & 
                      dt[["hba1c"]]>=37 & dt[["hba1c"]]<42] <- 1
    
    dt[["dmPre_4"]][dt[["dm_0_10"]]==0 & dt[["dm_10_inf"]]==0 & 
                      dt[["hba1c"]]>=42] <- 1
    
    
    ##################################################
    ###### check diabetes by baseline and incident
    ##################################################
    
    # generate diabetes events by 5-year interval
    dt[["dm_0_5"]] <- dt[["dm_5_10"]] <- dt[["dm_10_15"]] <- dt[["dm_15_20"]] <- dt[["dm_20_inf"]] <- rep(0, nrow(dt))
    
    dt[["dm_0_5"]][dt[[paste0("diabetes.b4", e)]]==1 & dt[[paste0("diabetes.b4", e, ".date")]] + 1461 > dt$cir.start &
                     dt[[paste0("diabetes.b4", e, ".date")]] < dt$cir.start + 365.25] <- 1
    
    dt[["dm_5_10"]][dt[[paste0("diabetes.b4", e)]]==1 & dt[[paste0("diabetes.b4", e, ".date")]] + 3287.25 > dt$cir.start &
                      dt[[paste0("diabetes.b4", e, ".date")]] + 1826.25 <= dt$cir.start + 365.25] <- 1
    
    dt[["dm_10_15"]][dt[[paste0("diabetes.b4", e)]]==1 & dt[[paste0("diabetes.b4", e, ".date")]] + 5113.5 > dt$cir.start &
                       dt[[paste0("diabetes.b4", e, ".date")]] + 3652.5 <= dt$cir.start + 365.25] <- 1
    
    dt[["dm_15_20"]][dt[[paste0("diabetes.b4", e)]]==1 & dt[[paste0("diabetes.b4", e, ".date")]] + 6939.75 > dt$cir.start &
                       dt[[paste0("diabetes.b4", e, ".date")]] + 5478.75 <= dt$cir.start + 365.25] <- 1
    
    dt[["dm_20_inf"]][dt[[paste0("diabetes.b4", e)]]==1 &
                        dt[[paste0("diabetes.b4", e, ".date")]] + 7305 <= dt$cir.start + 365.25] <- 1
    
    dt[["dm_10_20"]] <- rep(0, nrow(dt))
    dt[["dm_10_20"]][dt[[paste0("diabetes.b4", e)]]==1 & dt[[paste0("diabetes.b4", e, ".date")]] + 6939.75 > dt$cir.start &
                       dt[[paste0("diabetes.b4", e, ".date")]] + 3652.5 <= dt$cir.start + 365.25] <- 1
    
    
    ###################################################################
    dt[["dm_hba1c"]] <- NA
    dt[["dm_hba1c"]][dt[["dmPre_1"]]==1] <- "dmPre_1"
    dt[["dm_hba1c"]][dt[["dmPre_2"]]==1] <- "dmPre_2"
    dt[["dm_hba1c"]][dt[["dmPre_3"]]==1] <- "dmPre_3"
    dt[["dm_hba1c"]][dt[["dmPre_4"]]==1] <- "dmPre_4"
    dt[["dm_hba1c"]][dt[["dm_0_10"]]==1] <- "dm_0_10"
    dt[["dm_hba1c"]][dt[["dm_10_inf"]]==1] <- "dm_10_inf"
    dt[["dm_hba1c"]] <- factor(dt[["dm_hba1c"]],
                               levels=c("dmPre_2", "dmPre_1", "dmPre_3",
                                        "dmPre_4", "dm_0_10", "dm_10_inf"))
    ###################################################################
    
    # # divide T1 diabetes from baseline diabetes
    # dt[["dmT1_0_10"]] <- dt[["dmT1_10_20"]] <- dt[["dmT1_20_inf"]] <- 0
    # dt[["dmT2_0_10"]] <- dt[["dmT2_10_20"]] <- dt[["dmT2_20_inf"]] <- 0
    # 
    # dt[["dmT1_0_10"]][dt[["dm_0_10"]]==1 & dt$dmT1==1] <- 1
    # dt[["dmT1_10_20"]][dt[["dm_10_20"]]==1 & dt$dmT1==1] <- 1
    # dt[["dmT1_20_inf"]][dt[["dm_20_inf"]]==1 & dt$dmT1==1] <- 1
    # 
    # dt[["dmT2_0_10"]][dt[["dm_0_10"]]==1 & dt$dmT1==0] <- 1
    # dt[["dmT2_10_20"]][dt[["dm_10_20"]]==1 & dt$dmT1==0] <- 1
    # dt[["dmT2_20_inf"]][dt[["dm_20_inf"]]==1 & dt$dmT1==0] <- 1
    
  }
  
  # generate cancer ----
  if (e != "cancer") {

    # incident
    dt[["cancer_icd_0_1"]] <- dt[["cancer_icd_1_2"]] <- dt[["cancer_icd_2_3"]] <- dt[["cancer_icd_3_4"]] <- dt[["cancer_icd_4_5"]] <- dt[["cancer_icd_0_5"]] <- dt[["cancer_icd_5_inf"]] <- rep(0, nrow(dt))
    
    dt[["cancer_icd_0_1"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                             dt[[paste0("cancer_icd.b4", e, ".date")]] > dt$cir.start &
                             dt[[paste0("cancer_icd.b4", e, ".date")]] < dt$cir.start + 365.25] <- 1  
    dt[["cancer_icd_1_2"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 365.25 > dt$cir.start &
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 365.25 <= dt$cir.start + 365.25] <- 1   
    dt[["cancer_icd_2_3"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 730.5 > dt$cir.start &
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 730.5 <= dt$cir.start + 365.25] <- 1    
    dt[["cancer_icd_3_4"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 1095.75 > dt$cir.start &
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 1095.75 <= dt$cir.start + 365.25] <- 1  
    dt[["cancer_icd_4_5"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 1461 > dt$cir.start &
                             dt[[paste0("cancer_icd.b4", e, ".date")]] + 1461 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_icd_0_5"]][dt[["cancer_icd_0_1"]]==1 | dt[["cancer_icd_1_2"]]==1 | 
                             dt[["cancer_icd_2_3"]]==1 | dt[["cancer_icd_3_4"]]==1 |
                             dt[["cancer_icd_4_5"]]==1] <- 1
    dt[["cancer_icd_5_inf"]][dt[[paste0("cancer_icd.b4", e)]]==1 & 
                               dt[[paste0("cancer_icd.b4", e, ".date")]] + 1826.25 <= dt$cir.start + 365.25] <- 1
    
    # baseline
    dt[["cancer_bsl_1_2"]] <- dt[["cancer_bsl_2_3"]] <- dt[["cancer_bsl_3_4"]] <- dt[["cancer_bsl_4_5"]] <- dt[["cancer_bsl_0_5"]] <- rep(0, nrow(dt))
    
    dt[["cancer_bsl_1_2"]][dt[["cancer.baseline.all"]]==1 & 
                             dt[["cancer.baseline.all.date"]] + 365.25 > dt$cir.start &
                             dt[["cancer.baseline.all.date"]] + 365.25 <= dt$cir.start + 365.25] <- 1            
    dt[["cancer_bsl_2_3"]][dt[["cancer.baseline.all"]]==1 & 
                             dt[["cancer.baseline.all.date"]] + 730.5 > dt$cir.start &
                             dt[["cancer.baseline.all.date"]] + 730.5 <= dt$cir.start + 365.25] <- 1             
    dt[["cancer_bsl_3_4"]][dt[["cancer.baseline.all"]]==1 & 
                             dt[["cancer.baseline.all.date"]] + 1095.75 > dt$cir.start &
                             dt[["cancer.baseline.all.date"]] + 1095.75 <= dt$cir.start + 365.25] <- 1           
    dt[["cancer_bsl_4_5"]][dt[["cancer.baseline.all"]]==1 & 
                             dt[["cancer.baseline.all.date"]] + 1461 > dt$cir.start &
                             dt[["cancer.baseline.all.date"]] + 1461 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_bsl_0_5"]][dt[["cancer_bsl_1_2"]]==1 | dt[["cancer_bsl_2_3"]]==1 | 
                             dt[["cancer_bsl_3_4"]]==1 | dt[["cancer_bsl_4_5"]]==1] <- 1
    
    # combine incident cancer 5_inf to baseline cancer 5_10
    dt[["cancer_com_5_10"]] <- dt[["cancer_com_10_15"]] <- dt[["cancer_com_15_20"]] <- dt[["cancer_com_20_inf"]] <- dt[["cancer_com_5_inf"]] <- rep(0, nrow(dt))
    
    dt[["cancer_com_5_10"]][dt[["cancer.baseline.all"]]==1 &
                              dt[["cancer.baseline.all.date"]] + 3287.25 > dt$cir.start &
                              dt[["cancer.baseline.all.date"]] + 1826.25 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_com_5_10"]][dt[["cancer_icd_5_inf"]]==1] <- 1
    dt[["cancer_com_10_15"]][dt[["cancer.baseline.all"]]==1 &
                               dt[["cancer.baseline.all.date"]] + 5113.5 > dt$cir.start &
                               dt[["cancer.baseline.all.date"]] + 3652.5 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_com_15_20"]][dt[["cancer.baseline.all"]]==1 &
                               dt[["cancer.baseline.all.date"]] + 6939.75 > dt$cir.start &
                               dt[["cancer.baseline.all.date"]] + 5478.75 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_com_20_inf"]][dt[["cancer.baseline.all"]]==1 &
                                dt[["cancer.baseline.all.date"]] + 7305 <= dt$cir.start + 365.25] <- 1
    dt[["cancer_com_5_inf"]][dt[["cancer_com_5_10"]]==1 | dt[["cancer_com_10_15"]]==1 | 
                               dt[["cancer_com_15_20"]]==1 | dt[["cancer_com_20_inf"]]==1] <- 1
   }
  
  # generate xBeta ----
  
  dic <- list(pp="prim", sp="sec")
  
  coef <- c(cf[[tolower(e)]][[dic[[p]]]][[d]][["cf_b"]], 
            cf[[tolower(e)]][[dic[[p]]]][[d]][["cf_t"]])
  
  # remove intercept, original cancer, baseline diabetes coefficient
  coef <- coef[!grepl("Intercept|cancer|dm|smoker|race", names(coef))] 
  # this means to exclude the interaction with smoker as well
  # to include interaction with smoker, change |smoker| to |^smoker|
  
  if (rmcvd == TRUE) coef <- coef[!grepl("mi_|stroke_", names(coef))] # temporarily remove mi and stroke 
  
  
  ####
  # RW 2021-02-23
  if (!is.null(var)){
    
    # check interaction elements
    inters_int <- strsplit(grep("_int_|_noageint_", var, value = TRUE), "_int_|_noageint_")
    inters_int <- unique(unlist(inters_int))
    # add to var list
    inters <- paste(intersect(unique(c(var, inters_int)), names(coef)), collapse = "|")
    
    # inters =="", grepl(inters, names(coef)) will return all TRUE
    if (inters !="") coef <- coef[!grepl(inters, names(coef))]
  }

  ####
  
  # interaction
  # get the names of variables that interact with age
  int_names <- sapply(strsplit(grep("_int_", names(coef), value = TRUE), "_int_"), "[", 2)
  # generate the interaction terms
  for (int in int_names) {
    dt[[paste0("CurrAge_cent_int_", int)]] <- dt$CurrAge_cent * dt[[int]]
  }
  # in case there is interaction in input var
  int_names <- sapply(strsplit(grep("_int_", var, value = TRUE), "_int_"), "[", 2)
  # generate the interaction terms
  for (int in int_names) {
    dt[[paste0("CurrAge_cent_int_", int)]] <- dt$CurrAge_cent * dt[[int]]
  }
  # in case there is non-age interaction in input var
  int_names <- strsplit(grep("_noageint_", var, value = TRUE), "_noageint_")
  if (length(int_names)>0){
    for (i in 1:length(int_names)) {
      dt[[paste0(int_names[[i]][1],"_noageint_", int_names[[i]][2])]] <- 
        dt[[int_names[[i]][1]]]*dt[[int_names[[i]][2]]]
    } 
  }

  # select variables that exactly match the CTT equation
  # dt2 <- dt %>% select(names(coef)) # avoid to use dplyr because stepAIC conflit this
  dt2 <- subset(dt, select = names(coef))
  
  # convert to matrix for multiplication
  dt2 <- data.matrix(dt2)
  
  # generate Xbeta
  dt$Xbeta <- dt2[, names(coef)] %*% coef
  
  # model -----
  
  # fit survival model using Xbeta
  dt$dummy <- ifelse(dt[[paste0(e, ".status")]]==1, 1, 0)
  dt$tstart <- dt$tstart/365.25
  dt[[paste0(e, ".time")]] <- dt[[paste0(e, ".time")]]/365.25
  
  if (is.null(var)) {
    fm <- as.formula(paste0("Surv(tstart, ", e, ".time, dummy) ~ Xbeta")) 
  } else {
    covar <- paste(var, collapse = "+")
    fm <- as.formula(paste0("Surv(tstart, ", e, ".time, dummy) ~ Xbeta + ", covar))
  }
  
  if (d=="exp") {
    est <- phreg(fm, dist="weibull", shape=1, data=dt)
  } else if (d=="wei") {
    est <- phreg(fm, dist="weibull", data=dt)
  } else if (d=="gom") {
    est <- phreg(fm, dist="gompertz", param = "rate", data=dt)
  } else if (d=="cox")
    est <- coxph(fm, data = dt)
  
  return(est)
}  
  
  
  
  
  
