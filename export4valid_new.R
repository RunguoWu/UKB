rm(list=ls())

library(tidyverse)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

source(file.path(prep, "model_preparation.R"))

# TODO: 
# load(file.path(vali_data, "cf_cali_cancer.Rdata"))


# read all variables in cf file -------------------------------------------

cf <- readRDS(file.path(vali_data, "cf_cali.rds"))

.es <- c("mi", "stroke", "crv", "cancer", "vd", "nvd")

all_p_b <- all_p_t <- all_s_b <- all_s_t <- c()

for (e in .es) { # exp, wei, gom have the same variables
  .prim_b <- names(cf[[e]][["calibrated"]][["prim"]][["exp"]][["cf_b"]])
  all_p_b <- c(all_p_b, .prim_b)
  
  .prim_t <- names(cf[[e]][["calibrated"]][["prim"]][["exp"]][["cf_t"]])
  all_p_t <- c(all_p_t, .prim_t)
  
  .sec_b <- names(cf[[e]][["calibrated"]][["sec"]][["exp"]][["cf_b"]])
  all_s_b <- c(all_s_b, .sec_b)
  
  .sec_t <- names(cf[[e]][["calibrated"]][["sec"]][["exp"]][["cf_t"]])
  all_s_t <- c(all_s_t, .sec_t)
}

all_p_b <- c(all_p_b, names(cf$diabetes$prim$exp$cf_b))
all_p_t <- c(all_p_t, names(cf$diabetes$prim$exp$cf_t))
all_s_b <- c(all_s_b, names(cf$diabetes$prim$exp$cf_b))
all_s_t <- c(all_s_t, names(cf$diabetes$prim$exp$cf_t))

all_p_b <- unique(all_p_b)
all_p_t <- unique(all_p_t)
all_s_b <- unique(all_s_b)
all_s_t <- unique(all_s_t)


# export files for graphing -----------------------------------------------

.es <- c("MI", "stroke", "CRV", "cancer", "diabetes", "VD", "NVD")

.ps <- c("prim", "sec")

for (e in .es) {
  for (p in .ps) {

    if (p == "prim") dt <- subset(d4s, CVhist==0) else
      dt <- subset(d4s, CVhist>0)
    
    # generate ids for pp and sp separately
    dt$id_new <- 1:dim(dt)[1]

    if (e == "cancer") dt <- subset(dt, cancer.baseline.all==0)

    # TODO: whether limited to those with GP record only
    if (e == "diabetes") dt <- subset(dt, gp_record==1 & diabetes.fo.pre==0)

    fu_end <- dt[[paste0(e, ".time")]]/365.25
    
    event_status <- as.numeric(dt[[paste0(e, ".status")]])-1
    
    temp <- cbind(id_new = dt$id_new, fu_end = fu_end, event_status = event_status)

    temp <- as.data.frame(temp)

    if (e == "cancer") nam <- "cancer_icd" else if (e == "diabetes") nam <- "dm" else
    nam <- tolower(e)
    
    saveRDS(temp, file = file.path(vali_data, paste0("ukb_", nam, "_", p, ".rds")), compress = F)
  }
}


# baseline variables  ----------------------------------------------------

.ctt_var <- c("Intercept", "male", "smoker","race_afro", "race_otherna",   
             "cvd_only", "h_mi_raw_only", "h_pad_raw_only", "othchd_only", "number_events2", 
             "b_creann", "dias_bpS", "sys_bpS", "lhdl", "NEWB_LDL_CL_cent",
             "underweight", "overweight", "obese", "obese2", "obese3",
             "Txhypen", "CAN", "dm", "male_int_dm")
.graph <- c("over_65", "RG5")

.calib <- c("race_sa", "smk_ex", "smk_cur", "town1", "town2", "town4", "town5",
            "severe_mental_illness", "pa_low", "pa_high", "pa_mis") 
             
ukb_b_prim <- subset(d4s, CVhist==0, select = c(.ctt_var, .graph, .calib))
ukb_b_sec <- subset(d4s, CVhist>0, select = c(.ctt_var, .graph, .calib))

ukb_b_prim$id_new<- 1:dim(ukb_b_prim)[1]
ukb_b_sec$id_new<- 1:dim(ukb_b_sec)[1]

ukb_b_prim <- ukb_b_prim %>% relocate(id_new)
ukb_b_sec <- ukb_b_sec %>% relocate(id_new)

rownames(ukb_b_prim) <- NULL
rownames(ukb_b_sec) <- NULL

# generate the cf will call but ukb do not have as 0
for (i in setdiff(all_p_b, colnames(ukb_b_prim))) {
  ukb_b_prim[,i] <- rep(0,dim(ukb_b_prim)[1])
}

for (i in setdiff(all_s_b, colnames(ukb_b_sec))) {
  ukb_b_sec[,i] <- rep(0, dim(ukb_b_sec)[1])
}

# the base is for graphing
ukb_base_prim <- ukb_b_prim %>% 
  select("id_new", "male", "over_65", "RG5", "dm", "CAN") %>% rename("dm_baseline" = "dm") 

ukb_base_sec <- ukb_b_sec %>% 
  select("id_new", "male", "over_65", "RG5", "dm", "CAN") %>% rename("dm_baseline" = "dm") 

ukb_base_prim$RG5[ukb_base_prim$RG5=="<5"] <- 1
ukb_base_prim$RG5[ukb_base_prim$RG5=="[5,10)"] <- 2
ukb_base_prim$RG5[ukb_base_prim$RG5=="[10,15)"] <- 3
ukb_base_prim$RG5[ukb_base_prim$RG5=="[15,20)"] <- 4
ukb_base_prim$RG5[ukb_base_prim$RG5==">=20"] <- 5        

ukb_base_sec$RG5[ukb_base_sec$RG5=="<5"] <- 1
ukb_base_sec$RG5[ukb_base_sec$RG5=="[5,10)"] <- 2
ukb_base_sec$RG5[ukb_base_sec$RG5=="[10,15)"] <- 3
ukb_base_sec$RG5[ukb_base_sec$RG5=="[15,20)"] <- 4
ukb_base_sec$RG5[ukb_base_sec$RG5==">=20"] <- 5 

# events to predict

for (prim in c("prim", "sec")) {
  
  dt <- get(paste0("ukb_b_", prim))
  
  events_to_predict <- dt %>% select(id_new, CAN, dm) %>% 
    mutate(mi = 1, stroke = 1, crv = 1) %>% rename(cancer_icd = CAN) %>% 
    mutate(cancer_icd = 1 - cancer_icd) %>% mutate(dm = 1 - dm) %>% as.matrix()
  
  saveRDS(events_to_predict, 
          file = file.path(vali_data, paste0("events_to_predict_", prim, ".rds")), compress = F)
}


# delete old or useless vars.
ukb_b_prim <- ukb_b_prim %>% select(-over_65, -RG5, -CAN, -dm, -male_int_dm, -smoker)
ukb_b_sec <- ukb_b_sec %>% select(-over_65, -RG5, -CAN, -dm, -male_int_dm, -smoker)

# save ukb_b as matrix for simulation
ukb_b_prim <- as.matrix(ukb_b_prim)
saveRDS(ukb_b_prim, file = file.path(vali_data, "ukb_b_prim.rds"), compress = F)

ukb_b_sec <- as.matrix(ukb_b_sec)
saveRDS(ukb_b_sec, file = file.path(vali_data, "ukb_b_sec.rds"), compress = F)

# save the graph base as dt
saveRDS(ukb_base_prim, file = file.path(vali_data, "ukb_base_prim.rds"), compress = F)
saveRDS(ukb_base_sec, file = file.path(vali_data, "ukb_base_sec.rds"), compress = F)


# time varying variables --------------------------------------------------

ukb_t <- data.frame(cycle=0, 
                    CVhist=d4s$CVhist, 
                    CurrAge_cent=(d4s$age.recruit-60)/10)

# ukb_t$CurrAge_cent_int_dm <- ukb_t$CurrAge_cent*d4s$diabetes.fo.pre 

ukb_t$CurrAge_cent_int_NEWB_LDL_CL_cent <- ukb_t$CurrAge_cent*d4s$NEWB_LDL_CL_cent
ukb_t$CurrAge_cent_int_sys_bpS <- ukb_t$CurrAge_cent*d4s$sys_bpS
ukb_t$CurrAge_cent_int_smk_cur <- ukb_t$CurrAge_cent*d4s$smk_cur
ukb_t$CurrAge_cent_int_male <- ukb_t$CurrAge_cent*d4s$male

ukb_t$wstatina <- ukb_t$wstatinayr1 <- ukb_t$wstatinayr2And <- ukb_t$ind6_1 <- ukb_t$ind7_1 <- 0

# generate empty incident variable
ukb_t$crv_0_1 <- ukb_t$crv_1_2 <- ukb_t$crv_2_3 <- ukb_t$crv_3_inf <- 0
ukb_t$mi_0_1 <- ukb_t$mi_1_2 <- ukb_t$mi_2_3 <- ukb_t$mi_3_inf <- 0
ukb_t$stroke_0_1 <- ukb_t$stroke_1_2 <- ukb_t$stroke_2_3 <- ukb_t$stroke_3_inf <- 0

# TODO: cancer and diabetes
# cancer
## incident
ukb_t$cancer_icd_0_1 <- ukb_t$cancer_icd_1_2 <- ukb_t$cancer_icd_2_3 <- 
  ukb_t$cancer_icd_3_4 <- ukb_t$cancer_icd_4_5 <- ukb_t$cancer_icd_5_10 <- 
    ukb_t$cancer_icd_10_15 <- ukb_t$cancer_icd_15_20 <- ukb_t$cancer_icd_20_inf <- 0 
## baseline
ukb_t$cancer_bsl_1_2 <- ukb_t$cancer_bsl_2_3 <- 
  ukb_t$cancer_bsl_3_4 <- ukb_t$cancer_bsl_4_5 <- ukb_t$cancer_bsl_5_10 <- 
  ukb_t$cancer_bsl_10_15 <- ukb_t$cancer_bsl_15_20 <- ukb_t$cancer_bsl_20_inf <- 0 

# fill the baseline cancer event for the beginning of the cycle
ukb_t[["cancer_bsl_1_2"]][d4s[["cancer.baseline.all"]]==1 & 
                         d4s[["cancer.baseline.all.date"]] + 365.25 > d4s$recruit.date &
                         d4s[["cancer.baseline.all.date"]] + 365.25 <= d4s$recruit.date + 365.25] <- 1            
ukb_t[["cancer_bsl_2_3"]][d4s[["cancer.baseline.all"]]==1 & 
                         d4s[["cancer.baseline.all.date"]] + 730.5 > d4s$recruit.date &
                         d4s[["cancer.baseline.all.date"]] + 730.5 <= d4s$recruit.date + 365.25] <- 1             
ukb_t[["cancer_bsl_3_4"]][d4s[["cancer.baseline.all"]]==1 & 
                         d4s[["cancer.baseline.all.date"]] + 1095.75 > d4s$recruit.date &
                         d4s[["cancer.baseline.all.date"]] + 1095.75 <= d4s$recruit.date + 365.25] <- 1           
ukb_t[["cancer_bsl_4_5"]][d4s[["cancer.baseline.all"]]==1 & 
                         d4s[["cancer.baseline.all.date"]] + 1461 > d4s$recruit.date &
                         d4s[["cancer.baseline.all.date"]] + 1461 <= d4s$recruit.date + 365.25] <- 1
ukb_t[["cancer_bsl_5_10"]][d4s[["cancer.baseline.all"]]==1 &
                          d4s[["cancer.baseline.all.date"]] + 3287.25 > d4s$recruit.date &
                          d4s[["cancer.baseline.all.date"]] + 1826.25 <= d4s$recruit.date + 365.25] <- 1
ukb_t[["cancer_bsl_10_15"]][d4s[["cancer.baseline.all"]]==1 &
                           d4s[["cancer.baseline.all.date"]] + 5113.5 > d4s$recruit.date &
                           d4s[["cancer.baseline.all.date"]] + 3652.5 <= d4s$recruit.date + 365.25] <- 1
ukb_t[["cancer_bsl_15_20"]][d4s[["cancer.baseline.all"]]==1 &
                           d4s[["cancer.baseline.all.date"]] + 6939.75 > d4s$recruit.date &
                           d4s[["cancer.baseline.all.date"]] + 5478.75 <= d4s$recruit.date + 365.25] <- 1
ukb_t[["cancer_bsl_20_inf"]][d4s[["cancer.baseline.all"]]==1 &
                            d4s[["cancer.baseline.all.date"]] + 7305 <= d4s$recruit.date + 365.25] <- 1

# diabetes
ukb_t$dm_0_10 <- ukb_t$dm_10_inf <- 0

# fill the baseline diabetes event for the beginning of the cycle
ukb_t[["dm_0_10"]][d4s[["diabetes.fo.pre"]]==1 & 
                     d4s[["diabetes.fo.pre.date"]] + 3287.25 > d4s$recruit.date &
                     d4s[["diabetes.fo.pre.date"]] < d4s$recruit.date + 365.25] <- 1

ukb_t[["dm_10_inf"]][d4s[["diabetes.fo.pre"]]==1 &
                       d4s[["diabetes.fo.pre.date"]] + 3652.5 <= d4s$recruit.date + 365.25] <- 1

# death
ukb_t$vd <- ukb_t$nvd <- 0

ukb_t$d_without_crv <- ukb_t$d_without_mi <- ukb_t$d_without_stroke <- 
  ukb_t$d_without_cancer_icd <- ukb_t$d_without_dm <- 0

# Pre-define baseline time varying events ---------------------------------

# NA = no baseline event

# cancer
ukb_t$t_cancer <- NA 
# fill t at the beginning of the cycle for those with baseline disease
# because we define new simulated event at 0-1 first, and t start at 0
# so for bsl cancer start at 1-2 cycle, t start at 1, ex: 
# mx_cancerb$t[ukb_t[["cancer_bsl_1_2"]]==1] <- 1
# mx_cancerb$t[ukb_t[["cancer_bsl_2_3"]]==1] <- 2

for (i in 1:20){
  ukb_t$t_cancer[d4s[["cancer.baseline.all"]]==1 & 
                   d4s[["cancer.baseline.all.date"]] + 365.25*i > d4s$recruit.date &
                   d4s[["cancer.baseline.all.date"]] + 365.25*i <= d4s$recruit.date + 365.25] <- i  
}
ukb_t$t_cancer[d4s[["cancer.baseline.all"]]==1 & 
                 d4s[["cancer.baseline.all.date"]] + 7670.25 <= d4s$recruit.date + 365.25] <- 21

# diabetes
ukb_t$t_dm <- NA

for (i in 1:10){
  ukb_t$t_dm[d4s[["diabetes.fo.pre"]]==1 & 
           d4s[["diabetes.fo.pre.date"]] + 365.25*i > d4s$recruit.date &
           d4s[["diabetes.fo.pre.date"]] + 365.25*i <= d4s$recruit.date + 365.25] <- i  
}
ukb_t$t_dm[d4s[["diabetes.fo.pre"]]==1 & 
           d4s[["diabetes.fo.pre.date"]] + 4017.75 <= d4s$recruit.date + 365.25] <- 11

# above end and save -------------------------------------------------------
t_since_eb <- list()

# primary
n <- dim(subset(ukb_t, CVhist==0))[1]

ukb_t_prim <- ukb_t %>% filter(CVhist==0) %>% mutate(id_new=1:n) 

t_since_eb[["cancer_bsl"]] <- ukb_t_prim %>% 
  mutate(ids=id_new, t=t_cancer) %>% select(ids, t) %>% as.matrix()

t_since_eb[["dm"]] <- ukb_t_prim %>% 
  mutate(ids=id_new, t=t_dm) %>% select(ids, t) %>% as.matrix()

saveRDS(t_since_eb, file = file.path(vali_data, "t_since_eb_prim.rds"), compress = F)  

ukb_t_prim <- ukb_t_prim %>% select(-CVhist, -t_cancer, -t_dm) %>% relocate(id_new) %>% as.matrix()

saveRDS(ukb_t_prim, file = file.path(vali_data, "ukb_t_prim.rds"), compress = F)

# secondary
t_since_eb <- list()

n <- dim(subset(ukb_t, CVhist>0))[1]

ukb_t_sec <- ukb_t %>% filter(CVhist>0) %>% mutate(id_new=1:n)

t_since_eb[["cancer_bsl"]] <- ukb_t_sec %>% 
  mutate(ids=id_new, t=t_cancer) %>% select(ids, t) %>% as.matrix()

t_since_eb[["dm"]] <- ukb_t_sec %>% 
  mutate(ids=id_new, t=t_dm) %>% select(ids, t) %>% as.matrix()

saveRDS(t_since_eb, file = file.path(vali_data, "t_since_eb_sec.rds"), compress = F) 

ukb_t_sec <- ukb_t_sec %>% select(-CVhist, -t_cancer, -t_dm) %>% relocate(id_new) %>% as.matrix()

saveRDS(ukb_t_sec, file = file.path(vali_data, "ukb_t_sec.rds"), compress = F)


# MI/CRV calibration ------------------------------------------------------

# generate the MI/CRV calibration file
# according to PMG 2020-11 (4)
# CRV and MI in the same cycle
CRVb4MI <- 55
MIb4CRV <- 2967
p_prim <- 2*CRVb4MI/(CRVb4MI + MIb4CRV)

CRVb4MI <- 40
MIb4CRV <- 999
p_sec <- 2*CRVb4MI/(CRVb4MI + MIb4CRV)

p_crv <- list(prim = p_prim, sec = p_sec)
saveRDS(p_crv, file = file.path(vali_data, "p_crv.rds"), compress = F)


# Time-updated pre-defined check points -----------------------------------

vars_t <- list(
  # time-updated baseline
  b = list(cancer_bsl = c(2, 3, 4, 5, 10, 15, 20, Inf), 
           dm = c(10, Inf)),
  # time-updated within-simulation
  e = list(
    crv = c(1, 2, 3, Inf),
    mi = c(1, 2, 3, Inf),
    stroke = c(1, 2, 3, Inf),
    dm = c(10, Inf),
    cancer_icd = c(1, 2, 3, 4, 5, 10, 15, 20, Inf)
  )
)

saveRDS(vars_t, file =  file.path(vali_data, "vars_t.rds"), compress = F)

