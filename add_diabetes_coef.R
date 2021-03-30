library(survival)
library(tidyverse)
library(eha)

# set path
source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

# prepare data
source(file.path(prep, "model_preparation.R"))

# read in the original coefficient file
cf <- readRDS(file.path(vali_data, "cf.rds"))

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


# selecting model element -------------------------------------------------

# dt$dummy <- ifelse(dt$diabetes.status==1, 1, 0)
# 
# baseline <- "male + smk_ex + smk_cur + race_afro + race_otherna + race_sa + 
#              underweight + overweight + obese + obese2 + obese3 + 
#              Txhypen + NEWB_LDL_CL_cent + lhdl + sys_bpS + dias_bpS + b_creann + 
#              h_mi_raw_only + cvd_only + othchd_only + h_pad_raw_only + number_events2 + 
#              pa_low + pa_high + pa_mis + severe_mental_illness + town1 + town2 + town4 + town5"
# 
# cancer <- "cancer_icd_0_5 + cancer_bsl_0_5 + cancer_com_5_inf"
# 
# # full model
# 
# mi <- "mi_0_1 + mi_1_2 + mi_2_3 + mi_3_inf"
# stroke <- "stroke_0_1 + stroke_1_2 + stroke_2_3 + stroke_3_inf"
# crv <- "crv_0_1 + crv_1_2 + crv_2_3 + crv_3_inf"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#           paste(baseline,"age", mi, stroke, crv, cancer, sep = "+")))
# 
# cox_dm_all <- coxph(fm, dt, id=id)
# 
# # 3 v 2
# mi <- "mi_0_1 + mi_1_2 + mi_none + mi_3_inf"
# stroke <- "stroke_0_1 + stroke_1_2 + stroke_none + stroke_3_inf"
# crv <- "crv_0_1 + crv_1_2 + crv_none + crv_3_inf"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                        paste(baseline,"age", mi, stroke, crv, cancer, sep = "+")))
# 
# cox_dm_3v2 <- coxph(fm, dt, id=id)
# 
# # 23 v 1
# mi <- "mi_0_1 + mi_none + mi_2_inf"
# stroke <- "stroke_0_1 + stroke_none + stroke_2_inf"
# crv <- "crv_0_1 + crv_none + crv_2_inf"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                        paste(baseline,"age", mi, stroke, crv, cancer, sep = "+")))
# 
# cox_dm_23v1 <- coxph(fm, dt, id=id)
# 
# # 123 v 0
# mi <- "mi_none + mi_1_inf"
# stroke <- "stroke_none + stroke_1_inf"
# crv <- "crv_none + crv_1_inf"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                        paste(baseline,"age", mi, stroke, crv, cancer, sep = "+")))
# 
# cox_dm_123v0 <- coxph(fm, dt, id=id)
# 
# # yes v no
# mi <- "mi_0_inf"
# stroke <- "stroke_0_inf"
# crv <- "crv_0_inf"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                        paste(baseline,"age", mi, stroke, crv, cancer, sep = "+")))
# 
# cox_dm_yvn <- coxph(fm, dt, id=id)

# the result of the level collapse above is 
# crv_0_1 + crv_1_inf
# mi, stroke is not signficant in any case

# select model variables --------------------------------------------------

# # for auto selection
# dt[["cancer_0_5vs5_inf"]] <- ifelse(dt[["cancer_icd_0_5"]]==1, "icd_0_5", 
#                                     ifelse(dt[["cancer_bsl_0_5"]]==1, "bsl_0_5", 
#                                            ifelse(dt[["cancer_com_5_inf"]]==1, "5_inf", "none")))
# dt[["cancer_0_5vs5_inf"]] <- as.factor(dt[["cancer_0_5vs5_inf"]])
# dt[["cancer_0_5vs5_inf"]] <- relevel(dt[["cancer_0_5vs5_inf"]], ref = "none")
# 
# dt[["crv_0_1vs1_inf"]] <- ifelse(dt[["crv_0_1"]]==1, "0_1", 
#                                  ifelse(dt[["crv_1_inf"]]==1, "1_inf", "none"))
# dt[["crv_0_1vs1_inf"]] <- as.factor(dt[["crv_0_1vs1_inf"]])
# dt[["crv_0_1vs1_inf"]] <- relevel(dt[["crv_0_1vs1_inf"]], ref = "none")
# 
# baseline <- "male + smoke.imp + ethn3 + BMI_cat + PA_group + townsend.Q5 +
#              NEWB_LDL_CL_cent + lhdl + sys_bpS + dias_bpS + b_creann + 
#              Txhypen + severe_mental_illness + CVD"
# 
# cancer <- "cancer_0_5vs5_inf"
# mi <- "mi_0_inf"
# stroke <- "stroke_0_inf"
# crv <- "crv_0_1vs1_inf"
# 
# dt$dummy <- ifelse(dt$diabetes.status==1, 1, 0)
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                     paste(baseline,"age", mi, stroke, crv, cancer, sep = "+") ))
# 
# model0 <- coxph(fm, dt, id=id)
# 
# library(MASS)
# stepAIC(model0, test="Chisq", direction = "backward", scope = list(lower= ~ male + age))
# 
# # result
# # - stroke_0_inf           1 151939    0.25  0.616661    
# # - dias_bpS               1 151939    0.81  0.368503 
# # - mi_0_inf               1 151941    2.64  0.104102 
# # as stepAIC keep mi_0_inf, but manually remove it and the two above and stepAIC again
# baseline <- "male + smoke.imp + ethn3 + BMI_cat + PA_group + townsend.Q5 +
#              NEWB_LDL_CL_cent + lhdl + sys_bpS + b_creann + 
#              Txhypen + severe_mental_illness + CVD"
# 
# fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", 
#                        paste(baseline,"age", crv, cancer, sep = "+") ))
# 
# model0 <- coxph(fm, dt, id=id)
# 
# stepAIC(model0, test="Chisq", direction = "backward", scope = list(lower= ~ male + age))
# 
# # result
# # no new variable to be removed

# fit model ---------------------------------------------------------------

baseline <- "male + smk_ex + smk_cur + race_afro + race_otherna + race_sa + 
             underweight + overweight + obese + obese2 + obese3 + 
             Txhypen + NEWB_LDL_CL_cent + lhdl + sys_bpS + b_creann + 
             h_mi_raw_only + cvd_only + othchd_only + h_pad_raw_only + number_events2 + 
             pa_low + pa_high + pa_mis + severe_mental_illness + town1 + town2 + town4 + town5"

CRV <- "crv_0_1 + crv_1_inf"
cancer <- "cancer_icd_0_5 + cancer_bsl_0_5 + cancer_com_5_inf"

dt$dummy <- ifelse(dt$diabetes.status==1, 1, 0)
dt$tstart <- dt$tstart/365.25
dt$diabetes.time <- dt$diabetes.time/365.25

fm <- as.formula(paste("Surv(tstart, diabetes.time, dummy) ~", paste(baseline,"age", CRV, cancer, sep = "+") ))

cox_dm_final <- coxph(fm, dt, id=id)

# coxph.csv(cox_dm_final)
dm_est[["cox"]] <- cox_dm_final

exp_dm_final <- phreg(fm, dist="weibull", shape=1, data=dt)
wei_dm_final <- phreg(fm, dist="weibull", data=dt)
gom_dm_final <- phreg(fm, dist="gompertz", param = "rate", data=dt)

dm_est <- list(exp=exp_dm_final, wei=wei_dm_final, gom=gom_dm_final)

save(dm_est, file = file.path(vali_data, "dm_est.Rdata"))

# add coefficients into cf ------------------------------------------------

load(file.path(vali_data, "dm_est.Rdata"))

diabetes <- list("prim"=list())

for (i in c("exp", "wei", "gom")) {
  
  coef <- dm_est[[i]]$coefficients
  
  #baseline coef.
  cf_b <- coef[1:29]
  
  # rename time varying CRV coef.
  cf_t <- c("CurrAge_cent"=as.numeric(coef["age"]), 
            "crv_0_1"= as.numeric(coef["crv_0_1"]), 
            "crv_1_2"= as.numeric(coef["crv_1_inf"]),
            "crv_2_3"= as.numeric(coef["crv_1_inf"]),
            "crv_3_inf"= as.numeric(coef["crv_1_inf"]), 
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
            "cancer_bsl_20_inf"= as.numeric(coef["cancer_com_5_inf"])
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

saveRDS(cf, file = file.path(vali_data, "cf_cali.rds"))

# check AIC/BIC -----------------------------------------------------------

aic_bic_dm <- c()

for (i in c("exp", "wei", "gom")) {
  est <- dm_est[[i]]
  aic <- extractAIC(est)[2]
  bic <- extractAIC(est, k=log(est$n))[2]
  x <- c("AIC" = aic, "BIC" = bic)
  aic_bic_dm <- cbind(aic_bic_dm, x)
}
colnames(aic_bic_dm) <- c("exp", "wei", "gom")  

aic_bic_dm

# result: gom






