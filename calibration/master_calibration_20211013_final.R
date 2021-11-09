rm(list = ls())

library(foreach)
library(doSNOW)
library(snowfall)
library(parallel)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

source(file.path(cali,  "function_calib.R"))

# run add diabetes_coef.R first

# re-calibrate as Claire update cf file

# as diabetes and qol coefficients have no change, copy them into cf_new beforehand 

# use the cf that is used in function_calib.R

# calibration -------------------------------------------------------------

# don't calibrate diabetes model, because it is solely based on UKB

# model specification

# var to be calibrated based on the result of auto_sel_newvar.R/Rmd

pa <- c("pa_low", "pa_high", "pa_mis")
smk <- c("smk_ex", "smk_cur")
mental <- c("severe_mental_illness")
town <- c("town1", "town2", "town4", "town5")
race <- c("race_afro", "race_sa", "race_otherna")
dm <- c("dmPre_1", "dmPre_3", "dmPre_4", "dm_0_10", "dm_10_inf")
can <- c("cancer_icd_0_5", "cancer_bsl_0_5", "cancer_com_5_inf")
can_nvd <- c("cancer_icd_0_1", "cancer_icd_1_2", "cancer_icd_2_3", "cancer_icd_3_4", "cancer_icd_4_5", 
             "cancer_bsl_1_2", "cancer_bsl_2_3", "cancer_bsl_3_4", "cancer_bsl_4_5",
             "cancer_com_5_10", "cancer_com_10_15","cancer_com_15_20","cancer_com_20_inf")
sex <- "male"
age <- "CurrAge_cent"
t1dm <- "dmT1"
diet <- "unhealthy_diet"
pCVD <- c("h_mi_raw_only", "cvd_only", "h_pad_raw_only", "number_events2")

# calibrate sex for sp, and also age for CRV sp and cancer
# specification 
# following the table 3 on page 136 in document

var_mi_pp <- c(age, pa, smk, mental, race, t1dm, diet, dm, can, "NEWB_LDL_CL_cent",  
               "CurrAge_cent_int_dm_0_10", "CurrAge_cent_int_dm_10_inf")
# var_mi_sp <- c(age, smk, race, t1dm, dm, can) 
var_mi_sp <- c(age, sex, smk, race, t1dm, dm, can) # calib_20211018_rev.rds
# var_mi_sp <- c(age, smk, race, t1dm, dm, can, pCVD)

var_stroke_pp <- c(smk, mental, town, t1dm, diet, dm, can)
var_stroke_sp <- c(pa, smk, mental, town, t1dm, dm, can)

var_crv_pp <- c(pa, smk, race, dm, can)
# var_crv_sp <- c(age, smk, race, dm, can)
var_crv_sp <- c(age, sex, smk, race, dm, can) # calib_20211018_rev.rds
# var_crv_sp <- c(age, smk, race, dm, can, pCVD)


var_vd_pp <- c(age, smk, mental, town, diet, dm,
               "mi_0_1", "mi_1_2", "mi_2_inf",
               "stroke_0_1", "stroke_1_inf",
               "CurrAge_cent_int_dm_0_10",
               "CurrAge_cent_int_dm_10_inf",
               "mi_0_1_noageint_dm_0_10", "mi_0_1_noageint_dm_10_inf",
               "stroke_0_1_noageint_dm_0_10", "stroke_0_1_noageint_dm_10_inf")
# var_vd_pp <- c(age, smk, mental, town, diet, dm,
#                "mi_0_1", "mi_1_2", "mi_2_inf",
#                "stroke_0_1", "stroke_1_inf") # calib_20211018_rev.rds
# var_vd_pp <- c(age, smk, mental, town, diet, dm, "lnbcreann",
#                "lnbcreann_noageint_dm_0_10", "lnbcreann_noageint_dm_10_inf",
#                "CurrAge_cent_int_dm_0_10", "CurrAge_cent_int_dm_10_inf") # calib_20211020_vdpp.rds

var_vd_sp <- c(age, pa, smk, town, dm, can, "NEWB_LDL_CL_cent", 
               "CurrAge_cent_int_NEWB_LDL_CL_cent", 
               "CurrAge_cent_int_dm_0_10", "CurrAge_cent_int_dm_10_inf")

# var_nvd_pp <- c(age,sex, pa, smk, mental, town, race, t1dm, diet, dm, can_nvd, 
#                 "CurrAge_cent_int_dm_0_10", "CurrAge_cent_int_dm_10_inf", 
#                 "CurrAge_cent_int_male")
var_nvd_pp <- c(age,sex, pa, smk, mental, town, race, t1dm, diet, dm, can_nvd, 
                "CurrAge_cent_int_male") # calib_20211018_rev.rds

var_nvd_sp <- c(pa, smk, town, race, diet, dm, can_nvd)

var_cancer <- c(age, sex, pa, smk, race, t1dm, diet, dm, "CurrAge_cent_int_male",
                "CurrAge_cent_int_smk_cur", "CurrAge_cent_int_smk_ex")

mod <- cbind(rep(c(rep(c("MI", "stroke", "CRV", "VD", "NVD"),2),"cancer"),3),
             rep(c(rep("pp", 5), rep("sp", 6)), 3),
             c(rep("exp", 11), rep("wei", 11), rep("gom", 11)),
             rep(c("var_mi_pp", "var_stroke_pp", "var_crv_pp", "var_vd_pp",
                   "var_nvd_pp", "var_mi_sp", "var_stroke_sp", "var_crv_sp",
                   "var_vd_sp", "var_nvd_sp", "var_cancer"), 3))

# stk pp, crv pp and cancer
# so fit model for the five first
# # in mod, that is i = 2,13,24,  3,14,25, 11,22,33
# haveHad <- c(2,13,24,  3,14,25, 11,22,33)
# toRun <- setdiff(1:33, haveHad)
# # correct nvd pp and vd sp
# nvdppvdsp <- c(5, 16, 27, 9, 20, 31)
# 
# # for revised mi sp crv sp vd pp and nvd pp
# mispcrvspvdppnvdpp <- c(6, 17, 28, 4, 15, 26, 8, 19, 30, 5, 16, 27)
# 
# # only vd pp
# vdpp <- c(4, 15, 26)
# 
# mispcrvsp <- c(6, 17, 28, 8, 19, 30) 


est_list <- list()

sfInit(cpus = 12, parallel = T)
cl <- sfGetCluster()
clusterExport(cl=cl, list("var_mi_pp", "var_stroke_pp", "var_crv_pp", "var_vd_pp", 
                          "var_nvd_pp", "var_mi_sp", "var_stroke_sp", "var_crv_sp", 
                          "var_vd_sp", "var_nvd_sp", "var_cancer"))

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(stats))
clusterEvalQ(cl, library(eha))
registerDoSNOW(cl)

retval <- foreach (i = 20) %dopar% {
  
  
  e <- mod[i, 1]
  p <- mod[i, 2]
  d <- mod[i, 3]
  v <- mod[i, 4]
  
  var <- get(v)
  
  if (e=="VD" & p=="pp") est <- calib(e, p, d, var, rmcvd=TRUE) else
  est <- calib(e, p, d, var)
  
  # keep AIC/BIC and coef
  # coef <- est$coefficients
  # aic <- extractAIC(est)[2]
  # bic <- extractAIC(est, k=log(est$n))[2]
  # nam <- paste0(e, "_", p, "_", d)
  # result <- list(c(coef, aic=aic, bic=bic))
  # names(result) <- nam
  # return(result)
  
  # get whole regression object
  result <- list(est)
  nam <- paste0(e, "_", p, "_", d)
  names(result) <- nam
  return(result)
  
  
  # keep coef and CI
  # HR <- round(summary(est)$conf.int[,1], digits = 2)
  # ci1 <- round(summary(est)$conf.int[,3], digits = 2)
  # ci2 <- round(summary(est)$conf.int[,4], digits = 2)
  # CI <- paste0("(", ci1, "-", ci2, ")")
  # nam <- paste0(e, "_", p, "_", d)
  # HR <- paste(HR, CI, sep = " ")
  # HR <- data.frame(Var=rownames(summary(est)$conf.int), HR=HR)
  # result <- list(HR)
  # names(result) <- nam
  # return(result)
}

sfStop()

for (i in 1:length(retval)) {
  est_list[[names(retval[[i]])]] <- retval[[i]][[1]]
}
 

# save
# saveRDS(est_list, file = file.path(vali_data, "calib_20211013.rds"))
# saveRDS(est_list, file = file.path(vali_data, "calib_20211018_rev.rds"))
# saveRDS(est_list, file = file.path(vali_data, "calib_20211020_vdpp.rds"))
# saveRDS(est_list, file = file.path(vali_data, "calib_20211022_micrvsp.rds"))
 
# # for cox models
# saveRDS(est_list, file = file.path(vali_data, "calib_20210926_newDM_cox.rds"))
# for (i in 1:length(est_list)) {
#   
#   write.csv(est_list[[i]], file.path(output,paste0(names(est_list[i]), ".csv")))
#   
# }


# Incorporate into cf file ------------------------------------------------

# est_list <- readRDS(file.path(vali_data, "calib_20211020_micrvsp.rds"))

# edit based on the original cf from Claire, but with diabetes and qol
# so remove calibrated in case confusion
cf <- readRDS(file.path(vali_data, "cf_20210714.rds"))
cf$mi$calibrated <- NULL
cf$stroke$calibrated <- NULL
cf$crv$calibrated <- NULL
cf$cancer$calibrated <- NULL
cf$vd$calibrated <- NULL
cf$nvd$calibrated <- NULL


# define cancer levels in simulation coef.
can_0_5 <- c("cancer_icd_0_1", "cancer_icd_1_2", "cancer_icd_2_3", "cancer_icd_3_4", "cancer_icd_4_5", 
             "cancer_bsl_1_2", "cancer_bsl_2_3", "cancer_bsl_3_4", "cancer_bsl_4_5")

can_5_inf <- c("5_10", "10_15", "15_20", "20_inf")


dic <- list(pp="prim", sp="sec")

for (m in 1:dim(mod)[1]) {
  e <- mod[m, 1]
  p <- mod[m, 2]
  d <- mod[m, 3]
  v <- mod[m, 4]
  
  bsl <- get(v)
  
  # detect if cancer in var for non-NVD events
  detect_cancer <- 0
  if (e != "NVD" & sum(grepl("cancer", bsl))> 0) detect_cancer <- 1
  
  # remove cancer because we will respecify it later 
  bsl <- bsl[!grepl("cancer", bsl)]
  
  obj <- paste0(e, "_", p, "_", d) 
  
  cali <- est_list[[obj]]
  
  # copy coef----     
  alpha <- cf[[tolower(e)]][[dic[[p]]]][[d]]
  
  # remove intercept, baseline diabetes and its interactions, smoker, race, cancer----  
  # cf_b
  alpha[["cf_b"]] <- alpha[["cf_b"]][!grepl("Intercept|dm|smoker|race", names(alpha[["cf_b"]]))]

  # cf_t
  alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("cancer|dm|smoker", names(alpha[["cf_t"]]))] 
  
  # for VD pp, remove all mi and stroke
  # if (e=="VD" & p=="pp") alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("mi_|stroke_", names(alpha[["cf_t"]]))] 
  
  # remove interset variables if any remains
  inters <- paste(intersect(bsl, c(names(alpha[["cf_b"]]),names(alpha[["cf_t"]]))), collapse = "|")
  # inters =="", grepl(inters, names(coef)) will return all TRUE
  if (inters !="") {
    alpha[["cf_b"]] <- alpha[["cf_b"]][!grepl(inters, names(alpha[["cf_b"]]))]
    alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl(inters, names(alpha[["cf_t"]]))]
  }
  
  # add coef of calibrated variables ----
  # all add baseline covar, but add age or age interactions, if any, to cf_t
  for (i in bsl){
    if (!grepl("CurrAge_cent|dm_|dmP|^mi_|stroke_", i)) alpha[["cf_b"]][i] <- cali[i]/cali["Xbeta"] else
      alpha[["cf_t"]][i] <- cali[i]/cali["Xbeta"]
  }

  # for VD pp, edit mi and stroke
  # if (e=="VD" & p=="pp") {
  #   alpha[["cf_t"]]["mi_2_3"] <- alpha[["cf_t"]]["mi_3_inf"] <- alpha[["cf_t"]]["mi_2_inf"]
  #   alpha[["cf_t"]]["stroke_1_2"] <- alpha[["cf_t"]]["stroke_2_3"] <- 
  #     alpha[["cf_t"]]["stroke_3_inf"] <- alpha[["cf_t"]]["stroke_1_inf"]
  #   alpha[["cf_t"]] <- alpha[["cf_t"]][!grepl("mi_2_inf|stroke_1_inf", names(alpha[["cf_t"]]))] 
  # }
  
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
  } else if (detect_cancer==1){
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

# saveRDS(cf, file = file.path(vali_data, "cf_20210714_cali1016.rds"), compress = F)
# saveRDS(cf, file = file.path(vali_data, "cf_20210714_cali1018_rev.rds"), compress = F)
# saveRDS(cf, file = file.path(vali_data, "cf_20210714_cali1020_vdpp.rds"), compress = F)
saveRDS(cf, file = file.path(vali_data, "cf_20210714_cali1022_micrvsp.rds"), compress = F)

# select based on AIC/BIC -------------------------------------------------

# est_list <- readRDS(file.path(vali_data, "calib_20211013.rds"))

mod <- c("MI_pp","MI_sp", "stroke_pp","stroke_sp", "CRV_pp","CRV_sp",
         "VD_pp","VD_sp", "NVD_pp","NVD_sp", "cancer_sp")

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
aic_bic <- round(aic_bic, digits = 0)
colnames(aic_bic) <- nam
rownames(aic_bic) <- c("exp_aic", "wei_aic", "gom_aic", "exp_bic", "wei_bic", "gom_bic")

write.csv(aic_bic, file = file.path(output, "aic_bic211022_micrvsp.csv"))



# finally pick up what we want --------------------------------------------

# based on cf_20210714_cali1018
cf <- cf_20210714_cali1018_rev

# change VDpp with those from cf_20210714_cali1016
cf$vd$calibrated$prim <- cf_20210714_cali1016$vd$calibrated$prim

saveRDS(cf, file = file.path(vali_data, "cf_20210714_cali1022_final.rds"), compress = F)

