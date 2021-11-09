library(survival)
library(tidyverse)
library(eha)
library(survminer)

# set path
source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

# prepare data
source(file.path(prep, "model_preparation.R"))

# read in the original coefficient file
cf <- readRDS(file.path(vali_data, "cf_20210714.rds"))

# prepare for model -------------------------------------------------------

# only keep people without baseline diabetes 
d4s <- subset(d4s, gp_record==1 & diabetes.fo.pre==0)

# prepare data by annual cycle
dt <- survSplit(Surv(diabetes.time, diabetes.status)~., d4s, cut = c(365.25, 730.5, 1095.75, 1461, 1826.25, 2191.50, 2556.75, 2922, 3287.25, 3652.5), episode = "timegroup")

# generate a start and stop date for each circle
dt$cir.start <- dt$recruit.date + dt$tstart
dt$cir.stop <- dt$recruit.date + dt$diabetes.time

# age
# first timegroup since recruitment is 1 
# per 10 years age, centred at 60
dt$age <- (dt$age.recruit - 60 + dt$timegroup - 1)/10 

dt <- dt %>% mutate(age_cat = case_when(age < -1.5 ~ "<45", 
                                        age >=-1.5 & age < -1 ~ "45-50",
                                        age >=-1 & age < -0.5 ~ "50-55",
                                        age >=-0.5 & age < 0 ~ "55-60",
                                        age >=0 & age < 0.5 ~ "60-65",
                                        age >=0.5 & age < 1 ~ "65-70",
                                        age >=1 & age < 1.5 ~ "70-75",
                                        age >=1.5 ~ "75+")) 

# prepare event covariables

# generate covariate events by year ---- 
events <- c("MI", "stroke", "CRV") # generate cancer and dm covar separately

for (cov in events) {
  dt[[paste0(tolower(cov), "_0_1")]] <- 
    dt[[paste0(tolower(cov), "_1_2")]] <- 
    dt[[paste0(tolower(cov), "_2_3")]] <- 
    dt[[paste0(tolower(cov), "_3_inf")]] <- rep(0, nrow(dt))
  
  dt[[paste0(tolower(cov), "_0_1")]][dt[[paste0(cov, ".b4diabetes")]]==1 & 
        dt[[paste0(cov, ".b4diabetes.date")]] > dt$cir.start &
        dt[[paste0(cov, ".b4diabetes.date")]] < dt$cir.start + 365.25] <- 1 
  
  dt[[paste0(tolower(cov), "_1_2")]][dt[[paste0(cov, ".b4diabetes")]]==1 & 
        dt[[paste0(cov, ".b4diabetes.date")]] + 365.25 > dt$cir.start &
        dt[[paste0(cov, ".b4diabetes.date")]] + 365.25 <= dt$cir.start + 365.25] <- 1                      
  
  dt[[paste0(tolower(cov), "_2_3")]][dt[[paste0(cov, ".b4diabetes")]]==1 & 
        dt[[paste0(cov, ".b4diabetes.date")]] + 730.5 > dt$cir.start &
        dt[[paste0(cov, ".b4diabetes.date")]] + 730.5 <= dt$cir.start + 365.25] <- 1
  
  dt[[paste0(tolower(cov), "_3_inf")]][dt[[paste0(cov, ".b4diabetes")]]==1 & 
        dt[[paste0(cov, ".b4diabetes.date")]] + 1095.75 <= dt$cir.start + 365.25] <- 1

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

# cancer
dt$cancer_icd_0_5 <- dt$cancer_bsl_0_5 <- dt$cancer_com_5_inf <- rep(0, nrow(dt))
# incident cancer
dt$cancer_icd_0_5[dt$cancer_icd.b4diabetes==1 & 
                  dt$cancer_icd.b4diabetes.date + 1461 > dt$cir.start &
                  dt$cancer_icd.b4diabetes.date < dt$cir.start + 365.25] <- 1
# baseline cancer
dt$cancer_bsl_0_5[dt$cancer.baseline.all==1 & 
                  dt$cancer.baseline.all.date + 1461 > dt$cir.start &
                  dt$cancer.baseline.all.date < dt$cir.start + 365.25] <- 1

# combined cancer 5_inf
dt$cancer_com_5_inf[(dt$cancer_icd.b4diabetes==1 & 
                  dt$cancer_icd.b4diabetes.date + 1826.25 <= dt$cir.start + 365.25) |
                  (dt$cancer.baseline.all==1 & 
                   dt$cancer.baseline.all.date + 1826.25 <= dt$cir.start + 365.25)] <- 1

dt$dummy <- ifelse(dt$diabetes.status==1, 1, 0)


# fit model ---------------------------------------------------------------

# # before PMG in June
# baseline <- "male + smk_ex + smk_cur + race_afro + race_otherna + race_sa +
#              underweight + overweight + obese + obese2 + obese3 + 
#              Txhypen + NEWB_LDL_CL_cent + lhdl + sys_bpS + b_creann + 
#              h_mi_raw_only + cvd_only + othchd_only + h_pad_raw_only + number_events2 + 
#              pa_low + pa_high + pa_mis + severe_mental_illness + town1 + town2 + town4 + town5"
# 
# CRV <- "crv_0_1 + crv_1_inf"
# cancer <- "cancer_icd_0_5 + cancer_bsl_0_5 + cancer_com_5_inf"
# 
# dt$dummy <- ifelse(dt$diabetes.status==1, 1, 0)
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age", CRV, cancer, sep = "+") ))

# the selected specification after remove HbA1c>=48 and ESRD 
baseline <- "male + race_afro + race_otherna + race_sa +
             underweight + overweight + obese + obese2 + obese3 + 
             Txhypen + lhdl + sys_bpS + lnbcreann + hba1c + 
             h_mi_raw_only + cvd_only + othchd_only + h_pad_raw_only + number_events2 + 
             pa_low + pa_high + pa_mis + severe_mental_illness + 
             town1 + town2 + town4 + town5"

MI <- "mi_0_inf"
CRV <- "crv_0_1 + crv_1_inf"
cancer <- "cancer_icd_0_5 + cancer_bsl_0_5 + cancer_com_5_inf"

# # fm_1 <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age_cat", MI, CRV, cancer, sep = "+") ))
# # 
# # cox_dm <- coxph(fm_1, dt, id=id)
# # 
# # saveRDS(cox_dm, file = file.path(vali_data, "cox_dm_agecat.rds"))
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age","age*hba1c", MI, CRV, cancer, sep = "+") ))
# 
# cox_dm_int <- coxph(fm, dt, id=id)
# 
# coxph.csv(cox_dm_int)


# fit parametric models
dt$tstart <- dt$tstart/365.25
dt$diabetes.time <- dt$diabetes.time/365.25

fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age","age*hba1c", MI, CRV, cancer, sep = "+") ))

# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age", MI, CRV, cancer, sep = "+") ))
# gom_dm_noint <- phreg(fm, dist="gompertz", param = "rate", data=dt)

exp_dm_final <- phreg(fm, dist="weibull", shape=1, data=dt)
wei_dm_final <- phreg(fm, dist="weibull", data=dt)
gom_dm_final <- phreg(fm, dist="gompertz", param = "rate", data=dt)

dm_est <- list(exp=exp_dm_final, wei=wei_dm_final, gom=gom_dm_final, cox=cox_dm_int)

saveRDS(dm_est, file = file.path(vali_data, "dm_est_20210815.rds"))

# add coefficients into cf ------------------------------------------------

# dm_est <- readRDS(file.path(vali_data, "dm_est_20210815.rds"))

diabetes <- list("prim"=list())

for (i in c("exp", "wei", "gom")) {
  
  coef <- dm_est[[i]]$coefficients
  
  #baseline coef. var before age
  cf_b <- coef[1:(which(names(coef)=="age") - 1)]
  
  # rename time varying CRV coef.
  cf_t <- c("CurrAge_cent"=as.numeric(coef["age"]), 
            "crv_0_1"= as.numeric(coef["crv_0_1"]), 
            "crv_1_2"= as.numeric(coef["crv_1_inf"]),
            "crv_2_3"= as.numeric(coef["crv_1_inf"]),
            "crv_3_inf"= as.numeric(coef["crv_1_inf"]), 
            "mi_0_1"= as.numeric(coef["mi_0_inf"]), 
            "mi_1_2"= as.numeric(coef["mi_0_inf"]),
            "mi_2_3"= as.numeric(coef["mi_0_inf"]),
            "mi_3_inf"= as.numeric(coef["mi_0_inf"]),
            "cancer_icd_0_1"= as.numeric(coef["cancer_icd_0_5"]),
            "cancer_icd_1_2"= as.numeric(coef["cancer_icd_0_5"]),
            "cancer_icd_2_3"= as.numeric(coef["cancer_icd_0_5"]),
            "cancer_icd_3_4"= as.numeric(coef["cancer_icd_0_5"]),
            "cancer_icd_4_5"= as.numeric(coef["cancer_icd_0_5"]),
            "cancer_bsl_1_2"= as.numeric(coef["cancer_bsl_0_5"]),
            "cancer_bsl_2_3"= as.numeric(coef["cancer_bsl_0_5"]),
            "cancer_bsl_3_4"= as.numeric(coef["cancer_bsl_0_5"]),
            "cancer_bsl_4_5"= as.numeric(coef["cancer_bsl_0_5"]),
            "cancer_icd_5_10"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_icd_10_15"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_icd_15_20"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_icd_20_inf"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_bsl_5_10"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_bsl_10_15"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_bsl_15_20"= as.numeric(coef["cancer_com_5_inf"]),
            "cancer_bsl_20_inf"= as.numeric(coef["cancer_com_5_inf"]), 
            "CurrAge_cent_int_hba1c"= as.numeric(coef["hba1c:age"])
            )
  
  diabetes[["prim"]][[i]][["cf_b"]] <- cf_b
  diabetes[["prim"]][[i]][["cf_t"]] <- cf_t
}

# add intercept and shape
# exp
coef <- dm_est$exp$coefficients
diabetes[["prim"]][["exp"]][["cf_b"]] <- c("Intercept"=as.numeric(-coef["log(scale)"]), diabetes[["prim"]][["exp"]][["cf_b"]])
diabetes[["prim"]][["exp"]][["shape"]] <- NA

# wei
coef <- dm_est$wei$coefficients
diabetes[["prim"]][["wei"]][["cf_b"]] <- c("Intercept"=as.numeric(-exp(coef["log(shape)"])*coef["log(scale)"]), diabetes[["prim"]][["wei"]][["cf_b"]])
diabetes[["prim"]][["wei"]][["shape"]] <- as.numeric(exp(coef["log(shape)"]))

# gom
coef <- dm_est$gom$coefficients
diabetes[["prim"]][["gom"]][["cf_b"]] <- c("Intercept"=as.numeric(coef["log(level)"]), diabetes[["prim"]][["gom"]][["cf_b"]])
diabetes[["prim"]][["gom"]][["shape"]] <- as.numeric(coef["rate"])

# because we combine primary and secondary prevention populations, copy prim to sec
diabetes[["sec"]] <- diabetes[["prim"]]




# finally add the diabetes coef. into cf file
cf[["diabetes"]] <- diabetes

saveRDS(cf, file = file.path(vali_data, "cf_20210714.rds"))







