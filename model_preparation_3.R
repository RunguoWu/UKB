library(tidyverse)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

# based on ctt2, the imputed data
load(file.path(work_data, "ctt2.Rdata"))


# only keep useful variables from ctt2 ------------------------------------

d4s <- ctt2 %>% 
  select("id", "birth.year", "sex", "age.recruit", "age.group", 
         "recruit.date", "smoke.imp", "LDL.imp", "HDL.imp", "creatinine.imp", 
         "SBP.imp","DBP.imp", "ethn2","ethn3", "statin", "BMI_cat.imp","HBPtx", "othCHD", 
         "PVD", "qrisk_cat", "QRISK3_2017", "IMD.Q5", "townsend.imp",
         "MI.baseline", "MI.all.post", "MI.all.date", "stroke.baseline", 
         "stroke.all.post", "stroke.all.date", "CRV.all.post", "CRV.all.date", 
         "cancer.baseline.all", "cancer.baseline.all.date", "cancer.incident.only", 
         "cancer.incident.only.date",
         "diabetes.fo", "diabetes.fo.date", "diabetes.fo.pre", "diabetes.fo.pre.date", 
         "diabetes.fo.post", "diabetes.fo.post.date", "diabetes.fo.pre.duration",
         "diabetes.T1.fo.date", "diabetes.T2.fo.date","diabetes.T1.fo.pre", "diabetes.T2.fo.pre",
         "diabetes.T1.fo.pre.date", "diabetes.T2.fo.pre.date",
         "death.date", "death.vascular", "townsend.Q5", 
         "death.vascular.date", "death.nonvascular", "death.nonvascular.date", "CVhist",
         "PA_group", "neuroticism", "smoke.qrisk", "IMD.Q5", "gp_record", "severe_mental_illness"
  )


# generate endpoints and covar --------------------------------------------

.endpoints <- c("MI", "stroke", "CRV", "cancer", "diabetes", "VD", "NVD")

for (e in .endpoints){# generate event status, stop time and duration

  # a dictionary to indicate corresponding variables
  .es <- c("MI", "stroke", "CRV", "cancer", "diabetes", "VD", "NVD")
  .vars <- c("MI.all", "stroke.all", "CRV.all", "cancer.incident.only", 
             "diabetes.fo.post", "death.vascular", "death.nonvascular") 
  dic <- cbind(.es, .vars)
  
  # use 2016-03-31 the official reported censoring date as cut off
  cut <- "2016-03-31"
  
  d4s[[paste0(e, ".stop")]] <- d4s[[paste0(dic[.es==e, ".vars"], ".date")]]
  d4s[[paste0(e, ".stop")]][d4s[[paste0(e, ".stop")]]>cut] <- NA
  d4s[[paste0(e, ".status")]] <- ifelse(is.na(d4s[[paste0(e, ".stop")]]), 0, 1)
  
  # check death, if death happens before cutoff date and event, use death date as endpoint
  d4s$death.end <- 0
  d4s$death.end[d4s[[paste0(e, ".status")]]==0 & !is.na(d4s$death.date) & d4s$death.date <= cut] <- 1
  d4s[[paste0(e, ".stop")]][d4s$death.end==1] <- d4s$death.date[d4s$death.end==1]
  d4s[[paste0(e, ".status")]][d4s$death.end==1] <- 2
  d4s[[paste0(e, ".status")]] <- factor(d4s[[paste0(e, ".status")]])
  
  # code the rest NA date as cutoff date
  d4s[[paste0(e, ".stop")]][is.na(d4s[[paste0(e, ".stop")]])] <- cut
  
  # maybe a few events happened on the same day as recruitment, 
  # survSplit does not recognise time zero, therefore add 0.5 to them.
  d4s[[paste0(e, ".stop")]][d4s[[paste0(e, ".stop")]]==d4s$recruit.date] <- 
    d4s[[paste0(e, ".stop")]][d4s[[paste0(e, ".stop")]]==d4s$recruit.date] + 0.5
  
  # recruitment date as day 0 as the min is 0.5 
  d4s[[paste0(e, ".time")]] <- as.numeric(d4s[[paste0(e, ".stop")]] - d4s$recruit.date)
}

.endpoints <- c("MI", "stroke", "CRV", "cancer", "diabetes", "VD", "NVD")
.covariates <- c("MI", "stroke", "CRV", "cancer", "diabetes")

for (cov in .covariates) { # generate covariates before the endpoint
  for (e in setdiff(.endpoints, cov)) {
    
    # to use the same name as other previous scripts
    if (cov == "cancer") cov <- "cancer_icd"
    
    # a dictionary to indicate corresponding variables
    .covs <- c("MI", "stroke", "CRV", "cancer_icd", "diabetes", "VD", "NVD")
    .vars <- c("MI.all", "stroke.all", "CRV.all", "cancer.incident.only", 
               "diabetes.fo.post", "death.vascular", "death.nonvascular") 
    dic <- cbind(.covs, .vars)
    
    cov.date <- d4s[[paste0(dic[.covs==cov, ".vars"], ".date")]] 
    # don't include events on the same day of endpoint or censoring
    d4s[[paste0(cov, ".b4", e)]] <- ifelse(!is.na(cov.date) & 
                                             cov.date < d4s[[paste0(e, ".stop")]], 1, 0)
    
    if((cov=="MI" & e=="CRV") | (e=="VD" | e=="NVD") ){ 
      # include MI on the same day of CRV
      # for death endpoint, all events on the same day of endpoint are included
      # all events on the same day of censoring (including death not as the endpoint) are excluded
      d4s[[paste0(cov, ".b4", e)]][cov.date == d4s[[paste0(e, ".stop")]] & 
                                     d4s[[paste0(e, ".status")]]==1] <- 1
    }
    
    d4s[[paste0(cov, ".b4", e, ".date")]] <- cov.date
    d4s[[paste0(cov, ".b4", e, ".date")]][d4s[[paste0(cov, ".b4", e)]]==0] <- NA
    
    if(cov=="diabetes") {# incorporate baseline event into cycles according to their date
      d4s[[paste0(cov, ".b4", e, ".date")]][d4s$diabetes.fo.pre==1] <- 
        d4s$diabetes.fo.pre.date[d4s$diabetes.fo.pre==1]
      d4s[[paste0(cov, ".b4", e)]]<- ifelse(is.na(d4s[[paste0(cov, ".b4", e, ".date")]]), 0, 1)
    }
  }
}



# Recode some variables ---------------------------------------------------

d4s$CVD <- ifelse(d4s$CVhist==0, "None",
                  ifelse(d4s$CVhist==1 & d4s$MI.baseline==1, "MI only",  
                         ifelse(d4s$CVhist==1 & d4s$stroke.baseline==1, "Stroke only", 
                                ifelse(d4s$CVhist==1 & d4s$PVD==1, "PVD only", 
                                       ifelse(d4s$CVhist==1 & d4s$othCHD==1, "other CHD only", "Two or more")))))
d4s$CVD <- relevel(as.factor(d4s$CVD), ref = "None")

# BMI reference = 18.5-25
d4s$BMI_cat <- relevel(as.factor(d4s$BMI_cat.imp), ref = "18.5-25")

# ethn ref = white, only black and others
d4s$ethn2 <- relevel(as.factor(d4s$ethn2), ref = "White")

d4s <- d4s %>% 
  mutate(Intercept = rep(1, dim(d4s)[1]),
         id_new = 1:dim(d4s)[1],
         over_65 = ifelse(age.group=="65+", 1, 0),
         CAN=cancer.baseline.all,
         male = ifelse(sex=="Male", 1, 0), 
         smoker= ifelse(smoke.imp=="Current", 1, 0),
         underweight = ifelse(BMI_cat=="<18.5", 1, 0),
         overweight= ifelse(BMI_cat=="25-30", 1, 0),
         obese= ifelse(BMI_cat=="30-35", 1, 0),
         obese2= ifelse(BMI_cat=="35-40", 1, 0),
         obese3= ifelse(BMI_cat=="40+", 1, 0),
         race_afro= ifelse(ethn3=="Black", 1, 0),
         race_otherna= ifelse(ethn3=="Mixed or others", 1, 0),
         race_sa = ifelse(ethn3=="South Asian", 1, 0),
         cvd_only = ifelse(CVD=="Stroke only", 1, 0),
         h_mi_raw_only = ifelse(CVD=="MI only", 1, 0),
         h_pad_raw_only = ifelse(CVD=="PVD only", 1, 0),
         othchd_only = ifelse(CVD=="other CHD only", 1, 0),
         number_events2 = ifelse(CVD=="Two or more", 1, 0),
         NEWB_LDL_CL_cent=LDL.imp-3.6,
         lhdl=log(HDL.imp),
         Txhypen=HBPtx,
         sys_bpS=(SBP.imp-140)/20,
         dias_bpS=(DBP.imp-80)/10,
         b_creann=(creatinine.imp-80)/20,
         dm=diabetes.fo.pre,
         male_int_dm=male*dm,
         town1=ifelse(townsend.Q5==1, 1, 0),
         town2=ifelse(townsend.Q5==2, 1, 0),
         town4=ifelse(townsend.Q5==4, 1, 0),
         town5=ifelse(townsend.Q5==5, 1, 0),
         smk_ex=ifelse(smoke.imp=="Previous", 1, 0),         
         smk_cur=ifelse(smoke.imp=="Current", 1, 0),
         pa_low=ifelse(PA_group=="low", 1, 0),
         pa_mid=ifelse(PA_group=="moderate", 1, 0),
         pa_high=ifelse(PA_group=="high", 1, 0),
         pa_mis=ifelse(PA_group=="missing", 1, 0),
         RG5 = qrisk_cat
  )




# functions for model output ----------------------------------------------

# 95% CI
coxph.csv <- function(fit){ # input = output of coxph() 
  HR <- round(summary(fit)$conf.int[,1], digits = 2)
  CI <- paste(round(summary(fit)$conf.int[,3],digits = 2), round(summary(fit)$conf.int[,4],digits = 2), sep = "-")
  CI <- paste0("(", CI, ")")
  HR <- paste(HR, CI, sep = " ")
  HR <- data.frame(Var=rownames(summary(fit)$conf.int), HR=HR)
  write.csv(HR, paste0("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB", "/output/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

# 99% CI
coxph99.csv <- function(fit){ # input = output of coxph() 
  sm <- summary(fit, conf.int=0.99)
  HR <- round(sm$conf.int[,1], digits = 2)
  CI <- paste(round(sm$conf.int[,3],digits = 2), round(sm$conf.int[,4],digits = 2), sep = "-")
  CI <- paste0("(", CI, ")")
  HR <- paste(HR, CI, sep = " ")
  HR <- data.frame(Var=rownames(sm$conf.int), HR=HR)
  write.csv(HR, paste0("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB", "/output/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

# for status =0/1/2
coxph012.csv <- function(fit){ # input = output of coxph() 
  HR <- round(summary(fit)$conf.int[1:(length(summary(fit)$conf.int[,1])/2), 1], digits = 2)
  CI <- paste(round(summary(fit)$conf.int[1:(length(summary(fit)$conf.int[,1])/2), 3],digits = 2), round(summary(fit)$conf.int[1:(length(summary(fit)$conf.int[,1])/2), 4],digits = 2), sep = "-")
  CI <- paste0("(", CI, ")")
  HR <- paste(HR, CI, sep = " ")
  HR <- data.frame(Var=rownames(summary(fit)$cmap)[-1], HR=HR)
  write.csv(HR, paste0("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB", "/output/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

# for status =0/1/2
coxph99012.csv <- function(fit){ # input = output of coxph() 
  sm <- summary(fit, conf.int=0.99)
  HR <- round(sm$conf.int[1:(length(sm$conf.int[,1])/2), 1], digits = 2)
  CI <- paste(round(sm$conf.int[1:(length(sm$conf.int[,1])/2), 3],digits = 2), round(sm$conf.int[1:(length(sm$conf.int[,1])/2), 4],digits = 2), sep = "-")
  CI <- paste0("(", CI, ")")
  HR <- paste(HR, CI, sep = " ")
  HR <- data.frame(Var=rownames(sm$cmap)[-1], HR=HR)
  write.csv(HR, paste0("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB", "/output/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

# for eha model
phreg.exp <- function(fit){
  coef <- fit$coefficients
  se <- sqrt(diag(fit$var))
  up <- round(exp(coef + 1.96*se), digits = 2)
  low <- round(exp(coef - 1.96*se), digits = 2)
  CI <- paste0("(", low, " to ", up, ")")
  HR <- round(exp(coef), digits = 2)
  HR <- paste(HR, CI, sep = " ")
  names(HR) <- names(coef)
  write.csv(HR, paste0("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB", "/output/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

rm(list = ls.str(mode = "numeric"))


