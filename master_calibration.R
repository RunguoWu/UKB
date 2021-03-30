rm(list = ls())

library(foreach)
library(doSNOW)
library(snowfall)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

source(file.path(cali,  "function_calib.R"))

# run add diabetes_coef.R first
# cf <- readRDS(file.path(vali_data, "cf_cali.rds"))

# calibration -------------------------------------------------------------

# don't calibrate diabetes model, because it is solely based on UKB

# model specification

# var to be calibrated based on the result of auto_sel_newvar.R/Rmd

pa <- c("pa_low", "pa_high", "pa_mis")
smk <- c("smk_ex", "smk_cur")
mental <- c("severe_mental_illness")
town <- c("town1", "town2", "town4", "town5")
race <- c("race_afro", "race_sa", "race_otherna")
dm <- c("dm_0_10", "dm_10_inf")
can <- c("cancer_icd_0_5", "cancer_bsl_0_5", "cancer_com_5_inf")
can_nvd <- c("cancer_icd_0_1", "cancer_icd_1_2", "cancer_icd_2_3", "cancer_icd_3_4", "cancer_icd_4_5", 
             "cancer_bsl_1_2", "cancer_bsl_2_3", "cancer_bsl_3_4", "cancer_bsl_4_5",
             "cancer_com_5_10", "cancer_com_10_15","cancer_com_15_20","cancer_com_20_inf")
# RW 2021-02-23
sex <- "male"
age <- "CurrAge_cent"

# calibrate sex for sp, and also age for CRV sp and cancer
# specification 
var_mi_pp <- c(pa, smk, mental, race, dm)
var_mi_sp <- c(sex, smk, race, dm)

var_stroke_pp <- c(smk, mental, town, dm, can)
var_stroke_sp <- c(sex, pa, smk, mental, town, race, dm, can)

var_crv_pp <- c(pa, smk, race, dm, can)
var_crv_sp <- c(sex, age, smk, race, dm)

var_vd_pp <- c(smk, mental, town, dm)
var_vd_sp <- c(sex, pa, smk, town, race, dm, can)

var_nvd_pp <- c(pa, smk, mental, town, race, dm, can_nvd)
var_nvd_sp <- c(sex, pa, smk, town, race, dm, can_nvd)

var_cancer <- c(age, sex, pa, smk, race, dm)
# add interactions
var_cancer <- c(var_cancer, "CurrAge_cent*male", 
                "CurrAge_cent*smk_cur")


# mod <- cbind(rep(c(rep(c("MI", "stroke", "CRV", "VD", "NVD"),2),"cancer"),3),
#              rep(c(rep("pp", 5), rep("sp", 6)), 3),
#              c(rep("exp", 11), rep("wei", 11), rep("gom", 11)), 
#              rep(c("var_mi_pp", "var_stroke_pp", "var_crv_pp", "var_vd_pp", "var_nvd_pp", "var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp", "var_cancer"), 3))

# only calibrate sp
mod <- cbind(rep(c("MI", "stroke", "CRV", "VD", "NVD"),3),
             rep(rep("sp", 5), 3),
             c(rep("exp", 5), rep("wei", 5), rep("gom", 5)),
             rep(c("var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp"), 3))

# just for cancer
# mod <- cbind(rep("cancer",3),
#              rep("sp", 3),
#              c("exp","wei", "gom"), 
#              rep("var_cancer", 3))

sfInit(cpus = 8, parallel = T)
cl <- sfGetCluster()
# clusterExport(cl=cl, list("var_mi_pp", "var_stroke_pp", "var_crv_pp", "var_vd_pp", "var_nvd_pp", "var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp", "var_cancer"))
clusterExport(cl=cl, list("var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp"))
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(eha))
registerDoSNOW(cl)

# the machine memory may not be large enough to allow 33 to run at the same time
# each time run 8-9 model
retval <- foreach (i = 9:15) %dopar% {
  e <- mod[i, 1]
  p <- mod[i, 2]
  d <- mod[i, 3]
  v <- mod[i, 4]
  
  var <- get(v)
  
  est <- calib(e, p, d, var)
  coef <- est$coefficients
  aic <- extractAIC(est)[2]
  bic <- extractAIC(est, k=log(est$n))[2]
  nam <- paste0(e, "_", p, "_", d)
  result <- list(c(coef, aic=aic, bic=bic))
  names(result) <- nam
  return(result)
}

sfStop()

# est_list <- list()
# save
for (i in 1:length(retval)) {
  est_list[[names(retval[[i]])]] <- retval[[i]][[1]]
}

# # rename the interactions in cancer equations
# names(est_list$cancer_sp_exp)[14] <- "CurrAge_cent_int_male"
# names(est_list$cancer_sp_exp)[15] <- "CurrAge_cent_int_smk_cur"
# names(est_list$cancer_sp_wei)[14] <- "CurrAge_cent_int_male"
# names(est_list$cancer_sp_wei)[15] <- "CurrAge_cent_int_smk_cur"
# names(est_list$cancer_sp_gom)[14] <- "CurrAge_cent_int_male"
# names(est_list$cancer_sp_gom)[15] <- "CurrAge_cent_int_smk_cur"

# save(est_list, file = file.path(vali_data, "est_list.Rdata"))

save(est_list, file = file.path(vali_data, "est_list_sp_recali.Rdata"))

# save(est_list, file = file.path(vali_data, "est_list_cancer_recali.Rdata"))


# Incorporate into cf file ------------------------------------------------

# load(file.path(vali_data, "est_list_sp_recali.Rdata"))

# var to be calibrated
pa <- c("pa_low", "pa_high", "pa_mis")
smk <- c("smk_ex", "smk_cur")
mental <- c("severe_mental_illness")
town <- c("town1", "town2", "town4", "town5")
race <- c("race_afro", "race_sa", "race_otherna")
dm <- c("dm_0_10", "dm_10_inf")
# RW 2021-02-23
sex <- "male"
age <- "CurrAge_cent"

can_0_5 <- c("cancer_icd_0_1", "cancer_icd_1_2", "cancer_icd_2_3", "cancer_icd_3_4", "cancer_icd_4_5", 
             "cancer_bsl_1_2", "cancer_bsl_2_3", "cancer_bsl_3_4", "cancer_bsl_4_5")

can_5_inf <- c("5_10", "10_15", "15_20", "20_inf")


# specification 
var_mi_pp <- c(pa, smk, mental, race)
var_mi_sp <- c(sex, smk, race)

var_stroke_pp <- c(smk, mental, town)
var_stroke_sp <- c(sex, pa, smk, mental, town, race)

var_crv_pp <- c(pa, smk, race)
var_crv_sp <- c(sex, age, smk, race)

var_vd_pp <- c(smk, mental, town)
var_vd_sp <- c(sex, pa, smk, town, race)

var_nvd_pp <- c(pa, smk, mental, town, race)
var_nvd_sp <- c(sex, pa, smk, town, race)

var_cancer <- c(age, sex, pa, smk, race)
# add interactions
var_cancer <- c(var_cancer, "CurrAge_cent_int_male", 
                "CurrAge_cent_int_smk_cur")


# mod <- cbind(rep(c(rep(c("MI", "stroke", "CRV", "VD", "NVD"),2),"cancer"),3),
#              rep(c(rep("pp", 5), rep("sp", 6)), 3),
#              c(rep("exp", 11), rep("wei", 11), rep("gom", 11)), 
#              rep(c("var_mi_pp", "var_stroke_pp", "var_crv_pp", "var_vd_pp", "var_nvd_pp", "var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp", "var_cancer"), 3))

# only calibrate sp
# mod <- cbind(rep(c("MI", "stroke", "CRV", "VD", "NVD","cancer"),3),
#              rep(rep("sp", 6), 3),
#              c(rep("exp", 6), rep("wei", 6), rep("gom", 6)), 
#              rep(c("var_mi_sp", "var_stroke_sp", "var_crv_sp", "var_vd_sp", "var_nvd_sp", "var_cancer"), 3))

mod <- cbind(rep("cancer",3),
             rep("sp", 3),
             c("exp","wei", "gom"), 
             rep("var_cancer", 3))


dic <- list(pp="prim", sp="sec")

for (i in 1:dim(mod)[1]) {
  e <- mod[i, 1]
  p <- mod[i, 2]
  d <- mod[i, 3]
  v <- mod[i, 4]
  
  bsl <- get(v)
  
  obj <- paste0(e, "_", p, "_", d) 
  
  cali <- est_list[[obj]]
  
  # copy coef----     
  alpha <- cf[[tolower(e)]][[dic[[p]]]][[d]]
  
  # remove intercept, baseline diabetes and its interactions, smoker, race, cancer----  
  # cf_b
  alpha[["cf_b"]] <- alpha[["cf_b"]][!grepl("Intercept|dm|smoker|race", names(alpha[["cf_b"]]))]
  # for recalibrate sp, further remove male
  if (p=="sp") alpha[["cf_b"]] <- alpha[["cf_b"]][!grepl("male", names(alpha[["cf_b"]]))]
  
  # cf_t
  alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("cancer|dm|smoker", names(alpha[["cf_t"]]))] 
  
  # cancer remove age and age interactions
  if (e=="cancer") alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("CurrAge_cent", names(alpha[["cf_t"]]))] 
  
  # CRV sp remove age
  if (e=="CRV" & p=="sp") alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("CurrAge_cent", names(alpha[["cf_t"]]))]
  
  # add coef of calibrated variables ----
  # all add baseline covar, but add age or age interactions, if any, to cf_t
  for (i in bsl){
    if (!grepl("CurrAge_cent", i)) alpha[["cf_b"]][i] <- cali[i]/cali["Xbeta"] else
      alpha[["cf_t"]][i] <- cali[i]/cali["Xbeta"]
  }

  # all add diabetes covar
  for (i in dm) {
    alpha[["cf_t"]][i] <- cali[i]/cali["Xbeta"]
  }
  
  # add cancer covar depends on
  if (e=="NVD"){
    for (i in can_0_5) {
      alpha[["cf_t"]][i] <- cali[i]/cali["Xbeta"]
    }
    # separate icd and bsl, although share the same value
    for (j in can_5_inf){
      alpha[["cf_t"]][paste0("cancer_icd_", j)] <- cali[paste0("cancer_com_", j)]/cali["Xbeta"]
      alpha[["cf_t"]][paste0("cancer_bsl_", j)] <- cali[paste0("cancer_com_", j)]/cali["Xbeta"]
    }
  } else if (v %in% c( "var_stroke_pp", "var_crv_pp", "var_stroke_sp", "var_vd_sp")){
    for (i in can_0_5) {
      if (grepl("icd", i)){
        alpha[["cf_t"]][i] <- cali["cancer_icd_0_5"]/cali["Xbeta"]
      } else if (grepl("bsl", i)){
        alpha[["cf_t"]][i] <- cali["cancer_bsl_0_5"]/cali["Xbeta"]
      }
    }
    for (j in can_5_inf){
      alpha[["cf_t"]][paste0("cancer_icd_", j)] <- cali["cancer_com_5_inf"]/cali["Xbeta"]
      alpha[["cf_t"]][paste0("cancer_bsl_", j)] <- cali["cancer_com_5_inf"]/cali["Xbeta"]
    }  
  }
  
  # finally add intercept, shape and scale ----
  if(d=="exp"){
    # shape
    alpha[["shape"]] <- NA
    # new intercept
    alpha[["intercept_calibrated"]] <- as.numeric(-cali["log(scale)"])
    # scale
    alpha[["scaleby"]] <- as.numeric(cali["Xbeta"])
  } else if(d=="wei"){
    # shape
    alpha[["shape"]] <- as.numeric(exp(cali["log(shape)"]))
    # new intercept
    alpha[["intercept_calibrated"]] <- as.numeric(-exp(cali["log(shape)"])*cali["log(scale)"])
    # scale
    alpha[["scaleby"]] <- as.numeric(cali["Xbeta"])
  } else if(d=="gom") {
    # shape
    alpha[["shape"]] <- as.numeric(cali["rate"])
    # new intercept
    alpha[["intercept_calibrated"]] <- as.numeric(cali["log(level)"])
    # scale
    alpha[["scaleby"]] <- as.numeric(cali["Xbeta"]) 
  }
  
  cf[[tolower(e)]][["calibrated"]][[dic[[p]]]][[d]] <- alpha
}

# copy a cancer to prim/sec
cf$cancer$calibrated$prim <- cf$cancer$calibrated$sec

# saveRDS(cf, file = file.path(vali_data, "cf_cali.rds"), compress = F)

saveRDS(cf, file = file.path(vali_data, "cf_cali_sp_recali.rds"), compress = F)


# select based on AIC/BIC -------------------------------------------------

# load(file.path(vali_data, "est_list.Rdata"))

# mod <- c("MI_pp","MI_sp", "stroke_pp","stroke_sp", "CRV_pp","CRV_sp", 
#          "VD_pp","VD_sp", "NVD_pp","NVD_sp", "cancer_sp")

mod <- c("MI_sp", "stroke_sp", "CRV_sp", "VD_sp", "NVD_sp")
mod <- "cancer_sp"

aic <- c()
bic <- c()
nam <- c()  
for (i in mod) {
  # AIC
  exp <- est_list[[paste0(i, "_exp")]][["aic"]]
  wei <- est_list[[paste0(i, "_wei")]][["aic"]]
  gom <- est_list[[paste0(i, "_gom")]][["aic"]]
  bd <- c(exp, wei, gom)
  aic <- cbind(aic, bd)
  # BIC
  exp <- est_list[[paste0(i, "_exp")]][["bic"]]
  wei <- est_list[[paste0(i, "_wei")]][["bic"]]
  gom <- est_list[[paste0(i, "_gom")]][["bic"]]
  bd <- c(exp, wei, gom)
  bic <- cbind(bic, bd)
  # column names
  nam <- c(nam, i)
}
aic_bic <- rbind(aic, bic)
colnames(aic_bic) <- nam
rownames(aic_bic) <- c("exp_aic", "wei_aic", "gom_aic", "exp_bic", "wei_bic", "gom_bic")

write.csv(aic_bic, file = file.path(output, "calib_aic_bic.csv"))




