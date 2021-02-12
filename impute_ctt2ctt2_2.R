library(tidyverse)
library(QRISK3)
library(matrixStats)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

# based on ctt, the recoded dataset from the original bd
load(file.path(work_data, "ctt.Rdata"))

'%ni%' <- Negate('%in%')

osgb <- read.csv(file.path(work_data, "osgb.csv"), header = T)
load(file.path(work_data, "ctt2_imp.Rdata"))
corti <- read.csv(file.path(work_data, "cortico.csv"), header = F, encoding = "UTF-8")

ctt2 <- ctt

## Ad hoc imputation

# Ethnicity: NA -> white

ctt2$ethn2[is.na(ctt2$ethn2)] <- "White"

# Smoking status: impute to the majority by gender, age groups and education levels

ctt2 <- ctt2 %>% 
  group_by(sex, age.group, education) %>% 
  add_count(smoke) %>% 
  mutate(smoke.imp = if_else(is.na(smoke), smoke[which.max(n)], smoke)) %>% 
  select(-n) %>% ungroup()

# for QRISK calculation

# Current cigarettes per day: impute to the majority level, i.e. 10-20

ctt2$smoke.qrisk <- ifelse(ctt2$smoke.imp=="Never", 1, 
                          ifelse(ctt2$smoke.imp=="Previous", 2, 
                                 ifelse(ctt2$smoke.imp=="Current", 4, NA)))

# after tabulation the three level, moderate level has the most people: 15438; heavy 13512 and light 7213

# count cigarette per day
ctt2$smoke.qrisk[ctt2$smoke.qrisk==4 & ((bd$f.3456.0.0<10 & bd$f.3456.0.0>0) | bd$f.3456.0.0==-10)] <- 3 # less than 10. bd$f.3456.0.0==-10 means less than 1 

ctt2$smoke.qrisk[ctt2$smoke.qrisk==4 & bd$f.3456.0.0>=20] <- 5 # 20+

# Townsend score: regress on IMD score, years and countries. 
# Manually use Ordnance Survey (OSGB) rounded co-ordinates to check post codes 
# linked to ward Townsend score 2001 for missing for both Townsend and IMD. 
# Just use the authority level average score for the three missing for all these variables

ctt2$townsend <- bd$f.189.0.0

# regression first
# generate the combination of IMD scores from the three countries
ctt2$IMD <- ifelse(!is.na(bd$f.26410.0.0), bd$f.26410.0.0, 
                  ifelse(!is.na(bd$f.26427.0.0), bd$f.26427.0.0, 
                         ifelse(!is.na(bd$f.26426.0.0),bd$f.26426.0.0, NA)))

fit <- lm(townsend~IMD+IMD.source, ctt2)
ctt2$townsend.pred <- NA
ctt2$townsend.pred[!is.na(ctt2$IMD)] <- predict(fit, newdata=ctt2[!is.na(ctt2$IMD),])

ctt2$townsend.imp <- ctt2$townsend
ctt2$townsend.imp[is.na(ctt2$townsend.imp)] <- ctt2$townsend.pred[is.na(ctt2$townsend.imp)]

# manually lookup townsend score
osgb <- osgb %>% rename(id = f.eid) %>% select(id, townsend.2001)
ctt2 <- left_join(ctt2, osgb)
ctt2$townsend.imp[is.na(ctt2$townsend.imp)] <- ctt2$townsend.2001[is.na(ctt2$townsend.imp)]

# finally, check the assessment centres of the remaining three missing values
# f.54.0.0 = 11008 Bury, Townsend 0.0594
# f.54.0.0 = 11011 Bristol, Townsend 1.1034
# another one is missing for all values, will be deleted later
ctt2$townsend.imp[is.na(ctt2$townsend.imp) & bd$f.54.0.0==11008] <- 0.0594
ctt2$townsend.imp[is.na(ctt2$townsend.imp) & bd$f.54.0.0==11011] <- 1.1034

# Qrisk ethnicity preparation

white <- c("White", "British", "Irish", "Any other white background", "Prefer not to answer", "Do not know")
indian <- c("Indian")
pak <- c("Pakistani")
bang <- c("Bangladeshi")
chi <- c("Chinese")
oasia <- c("Asian or Asian British", "Any other Asian background")
bc <- c("Caribbean", "Black or Black British", "Any other Black background")
ba <- c("African")
oeth <- c("Mixed", "Other ethnic group", "White and Black Caribbean",
          "White and Black African", "White and Asian", "Any other mixed background")

text1 <- "ethnicity"
text2 <- "f.21000.0.0"
ctt2[[text1]] <- NA
ctt2[[text1]][bd[[text2]]%in% white] <- 1
ctt2[[text1]][bd[[text2]]%in% indian] <- 2
ctt2[[text1]][bd[[text2]]%in% pak] <- 3
ctt2[[text1]][bd[[text2]]%in% bang] <- 4
ctt2[[text1]][bd[[text2]]%in% oasia] <- 5
ctt2[[text1]][bd[[text2]]%in% bc] <- 6
ctt2[[text1]][bd[[text2]]%in% ba] <- 7
ctt2[[text1]][bd[[text2]]%in% chi] <- 8
ctt2[[text1]][bd[[text2]]%in% oeth] <- 9
# assume missing as white
ctt2$ethnicity[is.na(ctt2$ethnicity)] <- 1

# generate another ethnic vars
ctt2$ethn3 <- ifelse(ctt2$ethnicity ==1, "White", 
                ifelse(ctt2$ethnicity %in% c(6,7), "Black", 
                   ifelse(ctt2$ethnicity %in% c(2,3,4), "South Asian", "Mixed or others")))
ctt2$ethn3 <- as.factor(ctt2$ethn3)
ctt2$ethn3 <- relevel(ctt2$ethn3, ref = "White")

################################################################################
################################################################################
################################################################################
# # mute this chunk if we have already had imputed file
# 
# # multiple imputation
# 
# # generate variables for MI
# ctt2$weight <- bd$f.21002.0.0
# ctt2$height <- sqrt(ctt2$weight/ctt2$BMI)*100
# ctt2$cholesterol <- bd$f.30690.0.0
# ctt2$triglycerides <- bd$f.30870.0.0
# 
# # systolic BP read 1 & 2
# ctt2$SBP1 <- bd$f.4080.0.0
# ctt2$SBP1[is.na(ctt2$SBP1)] <- bd$f.93.0.0[is.na(ctt2$SBP1)]
# 
# ctt2$SBP2 <- bd$f.4080.0.1
# ctt2$SBP2[is.na(ctt2$SBP2)] <- bd$f.93.0.1[is.na(ctt2$SBP2)]
# 
# # Diastolic BP read 1 & 2
# ctt2$DBP1 <- bd$f.4079.0.0
# ctt2$DBP1[is.na(ctt2$DBP1)] <- bd$f.94.0.0[is.na(ctt2$DBP1)]
# 
# ctt2$DBP2 <- bd$f.4079.0.1
# ctt2$DBP2[is.na(ctt2$DBP2)] <- bd$f.94.0.1[is.na(ctt2$DBP2)]
# 
# # select vars for imputation
# 
# ctt2_imp <- ctt2 %>% select(id, weight, height, LDL, HDL, cholesterol, triglycerides, 
#               creatinine, SBP1, SBP2, DBP1, DBP2, age.recruit, sex, ethn2, smoke.imp, 
#               CVhist, hypertension, BPMed, statin, diabetes.fo.pre)
# 
# 
# ## Deletion 
# # Delete the only one case with nearly all information missing
# # id=3462265
# # Some participants quitted from UKB, so we should delete them from analysis as well.
# 
# droplist <- c(3462265, 1330343, 2376313, 3075648, 3123970, 3385704, 3387728, 
#               3809403, 4236833, 4602053, 5063590, 5158442)
# 
# '%ni%' <- Negate('%in%')
# 
# ctt2_imp <- ctt2_imp[ctt2_imp$id %ni% droplist, ]
# 
# # simplify auxiliary variables
# ctt2_imp$CVhist <- ifelse(ctt2_imp$CVhist==0,0,1) 
# 
# ctt2_imp$sex <- factor(ctt2_imp$sex, ordered = F)
# ctt2_imp$CVhist <- as.factor(ctt2_imp$CVhist)
# ctt2_imp$hypertension <- as.factor(ctt2_imp$hypertension)
# ctt2_imp$BPMed <- as.factor(ctt2_imp$BPMed)
# ctt2_imp$statin <- as.factor(ctt2_imp$statin)
# ctt2_imp$diabetes.fo.pre <- as.factor(ctt2_imp$diabetes.fo.pre)
# ctt2_imp$smoke.imp <- factor(ctt2_imp$smoke.imp, ordered = F)
# 
# save(ctt2_imp, file = "ctt2_imp.RData")
# 
# q("yes")
# 
# ## load it in a separate phase, may make the process safer as it will take days
# 
# load(file.path(work_data, "ctt2_imp.Rdata"))
# 
# library(mice)
# library(miceadds)
# library(tidyverse)
# 
# imp <- mice(ctt2_imp, maxit=0)
# predM = imp$predictorMatrix
# meth = imp$method
# 
# # ID won't predict other var.
# predM[, c("id")]=0
# 
# imp <- mice(ctt2_imp, m=20, maxit = 10, predictorMatrix = predM, method = meth, print = T)
# 
# imp_list <- miceadds::mids2datlist(imp)
# 
# save(imp, file = "imp.Rdata")
# save(imp_list, file = "implist.Rdata")
# 
# # all imputations to horizontal form 
# imp_r <- complete(imp, "broad")
# 
# # average of the imputations
# ctt2_imp$weight.imp <- rowMeans(imp_r[, grepl("weight", names(imp_r))])
# ctt2_imp$height.imp <- rowMeans(imp_r[, grepl("height", names(imp_r))])
# ctt2_imp$LDL.imp <- rowMeans(imp_r[, grepl("LDL", names(imp_r))])
# ctt2_imp$HDL.imp <- rowMeans(imp_r[, grepl("HDL", names(imp_r))])
# ctt2_imp$cholesterol.imp <- rowMeans(imp_r[, grepl("cholesterol", names(imp_r))])
# ctt2_imp$triglycerides.imp <- rowMeans(imp_r[, grepl("triglycerides", names(imp_r))])
# ctt2_imp$creatinine.imp <- rowMeans(imp_r[, grepl("creatinine", names(imp_r))])
# ctt2_imp$SBP1.imp <- rowMeans(imp_r[, grepl("SBP1", names(imp_r))])
# ctt2_imp$SBP2.imp <- rowMeans(imp_r[, grepl("SBP2", names(imp_r))])
# ctt2_imp$DBP1.imp <- rowMeans(imp_r[, grepl("DBP1", names(imp_r))])
# ctt2_imp$DBP2.imp <- rowMeans(imp_r[, grepl("DBP2", names(imp_r))])
# 
# save(ctt2_imp, file = "ctt2_imp.Rdata")
################################################################################
################################################################################
################################################################################

### after imputation

# drop one case without any information and other cases who want to quit
# to match the imputed data
droplist <- c(1102194, 1330343, 2376313, 2610320, 3039868, 3075648, 3123970,  3385704, 3387728, 3462265, 3809403, 4236833, 4602053, 4803604, 4890498, 5063590, 5158442)

ctt2 <- ctt2[ctt2$id %ni% droplist, ]

ctt2_imp <- ctt2_imp %>% select(id, weight.imp, height.imp, LDL.imp, HDL.imp, 
        cholesterol.imp, triglycerides.imp, creatinine.imp, SBP1.imp, SBP2.imp, 
        DBP1.imp, DBP2.imp)

ctt2 <- left_join(ctt2, ctt2_imp)

# combine BP
ctt2$SBP.imp <- (ctt2$SBP1.imp+ctt2$SBP2.imp)/2

ctt2$DBP.imp <- (ctt2$DBP1.imp+ctt2$DBP2.imp)/2

# BMI
ctt2$BMI.imp <- ctt2$weight.imp/(ctt2$height.imp/100)^2

# BMI category
text1 <- "BMI_cat.imp"
text2 <- "BMI.imp"
ctt2[[text1]] <- NA
ctt2[[text1]][ctt2[[text2]]<18.5] <- "<18.5"
ctt2[[text1]][ctt2[[text2]]>=18.5 & ctt2[[text2]]<25] <- "18.5-25"
ctt2[[text1]][ctt2[[text2]]>=25 & ctt2[[text2]]<30] <- "25-30"
ctt2[[text1]][ctt2[[text2]]>=30 & ctt2[[text2]]<35] <- "30-35"
ctt2[[text1]][ctt2[[text2]]>=35 & ctt2[[text2]]<40] <- "35-40"
ctt2[[text1]][ctt2[[text2]]>=40] <- "40+"

rm(ctt2_imp, fit, osgb)


################# calculating QRISK3 based on imputed data ###################

d4q <- data.frame(patid=ctt2$id, gender=ctt2$sex, age=ctt2$age.recruit, 
        weight=ctt2$weight.imp, townsend=ctt2$townsend.imp, height=ctt2$height.imp, 
        smoke=ctt2$smoke.qrisk, ethnicity= ctt2$ethnicity, blood_pressure_treatment=ctt2$BPMed)

# recode
d4q$gender <- ifelse(d4q$gender=="Female", 1, 0)
d4q$age <- as.numeric(d4q$age)

# cholesterol/HDL
d4q$cholesterol_HDL_ratio <- ctt2$cholesterol.imp/ctt2$HDL.imp

d4q$systolic_blood_pressure=(ctt2$SBP1.imp+ctt2$SBP2.imp)/2

# sd of SBP
d4q$std_systolic_blood_pressure <- rowSds(as.matrix(ctt2[, c("SBP1.imp", "SBP2.imp")]), na.rm = T)

# diabetes 1
d4q$diabetes1 <- ctt2$diabetes.T1.fo.pre
d4q$diabetes2 <- ctt2$diabetes.T2.fo.pre

################################################################################
# those below need "bd", without cases dropped
# so deal with it separately
# only apply to "bd" or "ctt"

temp <- data.frame(patid=bd$f.eid)

# atrial fibrillation
# use first occurence data
# ICD: I48
temp$atrial_fibrillation <- ifelse(bd$f.131350.0.0<ctt$recruit.date, 1, 0)
temp$atrial_fibrillation[is.na(temp$atrial_fibrillation)] <- 0

# relative angina or heart attack
# use "Heart disease" in UKB instead of ...
# father, mother, siblings
# no age known, so assume <60
temp$heart_attack_relative <- 0
for (j in 0:9) { # father illness
  text2 <- paste0("f.20107.0.", j)
  temp$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}

for (j in 0:10) { # mother illness
  text2 <- paste0("f.20110.0.", j)
  temp$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}

for (j in 0:11) { # sibling illness
  text2 <- paste0("f.20111.0.", j)
  temp$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}
# only with additional information about father/mother's age and age at death, 
# we can only assure some heart disease happened <60, but we cannot rule out 
# which heart disease that definitely happened >=60 


# ED
# only in nurse interview
# erectile dysfunction/impotence 1518
temp$erectile_disfunction <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  temp$erectile_disfunction[bd[[text2]]== 1518] <- 1
}

# refer to drug only also
# reference: Alice Carter. et al. Education inequalities in statin treatment
ed <- c(1140869100, 1140883010, 1141168936, 1141168944, 1141168946, 1141168948, 
        1141187810, 1141187814, 1141187818, 1141192248, 1141192256, 1141192258, 1141192260)
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  temp$erectile_disfunction[bd[[text2]] %in% ed] <- 1
}

# rheumatoid arthritis
# use first occurence data
# ICD: M05, M06
temp$rheumatoid_arthritis <- ifelse(bd$f.131848.0.0<ctt$recruit.date | bd$f.131850.0.0<ctt$recruit.date, 1, 0)
temp$rheumatoid_arthritis[is.na(temp$rheumatoid_arthritis)] <- 0

# migraine
# use first occurrence data
# G43: bd$f.131052.0.0

temp$migraine <- ifelse(bd$f.131052.0.0<ctt$recruit.date, 1, 0)
temp$migraine[is.na(temp$migraine)] <- 0

# systemic lupus erythematosis
# M32: bd$f.131894.0.0
temp$systemic_lupus_erythematosis <- ifelse(bd$f.131894.0.0<ctt$recruit.date, 1, 0)
temp$systemic_lupus_erythematosis[is.na(temp$systemic_lupus_erythematosis)] <- 0

# severe mental illness
# F20: bd$f.130874.0.0
# F23: bd$f.130880.0.0
# F31: bd$f.130892.0.0
# F32: bd$f.130894.0.0
# F33: bd$f.130896.0.0

temp$severe_mental_illness <- ifelse(bd$f.130874.0.0<ctt$recruit.date |
                                       bd$f.130880.0.0<ctt$recruit.date |
                                       bd$f.130892.0.0<ctt$recruit.date |
                                       bd$f.130894.0.0<ctt$recruit.date |
                                       bd$f.130896.0.0<ctt$recruit.date, 1, 0)

temp$severe_mental_illness[is.na(temp$severe_mental_illness)] <- 0
# assumption needed 

# ICD 10 only
# N183 N184 N185 CKD stage 3-5
# ICD 9 has no specified stage of CKD
# baseline interview has no specified stage of CKD
# use HSE record only
# as a result, no satisfying CKD recorded before recruitment
# i guess it is because the classification of CKD is new, so all people with CKD
# stage information is later than recruitment, i.e. 2006-2010
# so we revise the code, use ICD10 N18: chronic renal failure
temp$chronic_kidney_disease <- ifelse(bd$f.132032.0.0<ctt$recruit.date, 1, 0)
temp$chronic_kidney_disease[is.na(temp$chronic_kidney_disease)] <- 0

# steroid
# medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment  
# import the code list 
# the first element is coded incorrectly
# revise it manually 
corti$V1[1] <- 1140853854

corti <- as.numeric(corti$V1)
temp$regular_steroid_tablets <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  temp$regular_steroid_tablets[bd[[text2]] %in% corti] <- 1
}

# atypical psychotics
# medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment 
ap <- c(1140867420, 1140867432, 1140867444, 1140927956, 1140927970, 1140928916, 
        1141152848, 1141152860, 1141153490, 1141167976, 1141177762, 1141195974, 1141202024)

temp$atypical_antipsy <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  temp$atypical_antipsy[bd[[text2]] %in% ap] <- 1
}
##############################################################################

# join temp back to d4q
d4q <- left_join(d4q, temp)
rm(temp)

# calculate QRISK with the package"QRISK3"
# generate Q risk score
qrisk <- QRISK3_2017(data=d4q, 
                     patid = "patid", 
                     gender = "gender", 
                     age = "age", 
                     atrial_fibrillation = "atrial_fibrillation", 
                     atypical_antipsy = "atypical_antipsy", 
                     regular_steroid_tablets = "regular_steroid_tablets", 
                     erectile_disfunction = "erectile_disfunction", 
                     migraine = "migraine", 
                     rheumatoid_arthritis = "rheumatoid_arthritis", 
                     chronic_kidney_disease = "chronic_kidney_disease", 
                     severe_mental_illness = "severe_mental_illness", 
                     systemic_lupus_erythematosis = "systemic_lupus_erythematosis", 
                     blood_pressure_treatment = "blood_pressure_treatment", 
                     diabetes1 = "diabetes1", 
                     diabetes2 = "diabetes2", 
                     weight = "weight", 
                     height = "height", 
                     ethiniciy = "ethnicity", 
                     heart_attack_relative = "heart_attack_relative", 
                     cholesterol_HDL_ratio = "cholesterol_HDL_ratio", 
                     systolic_blood_pressure = "systolic_blood_pressure", 
                     std_systolic_blood_pressure = "std_systolic_blood_pressure", 
                     smoke = "smoke", 
                     townsend = "townsend")

# join the score back to ctt2
qrisk <- qrisk %>% select("patid", "QRISK3_2017") %>% rename(id=patid)
ctt2 <- left_join(ctt2, qrisk)

ctt2$qrisk_cat <- ifelse(ctt2$QRISK3_2017<5, "<5", 
                    ifelse(ctt2$QRISK3_2017>=5 & ctt2$QRISK3_2017<10, "[5,10)", 
                      ifelse(ctt2$QRISK3_2017>=10 & ctt2$QRISK3_2017<15, "[10,15)",
                        ifelse(ctt2$QRISK3_2017>=15 & ctt2$QRISK3_2017<20, "[15,20)", ">=20"))))

# add severe mental illness
ctt2 <- merge(ctt2, d4q[, c("patid", "severe_mental_illness")], by.x = "id", by.y = "patid"  )

# save(d4q, file = paste0(wd, "/processed_data/data4Qrisk.Rdata"))
rm(d4q, qrisk)

# generate excercise indicator
# >=150 min/w moderate or 75 min/week vigorous PA
# 22037 walk min/w
# 22038 moderate PA min/w
# 22039 vigorous PA min/w
# bd %>% mutate(exercise = if_else(bd$f.22037.0.0>=150, 1, 
#                                   if_else(bd$f.22038.0.0>=150, 1,
#                                     if_else(bd$f.22039.0.0>=75, 1, 0)))) %>% 
#   count(exercise)
# use the IPAQ group
# 22032
temp <- bd %>% select(f.eid, f.22032.0.0) %>% 
  rename(id = f.eid, PA_group = f.22032.0.0)
ctt2 <- left_join(ctt2, temp)

ctt2$PA_group[is.na(ctt2$PA_group)] <- "missing"
ctt2$PA_group <- as.factor(ctt2$PA_group)
ctt2$PA_group <- relevel(ctt2$PA_group, ref = "missing")

# neuroticism score
temp <- bd %>% select(f.eid, f.20127.0.0) %>% 
  rename(id = f.eid, neuroticism = f.20127.0.0)
ctt2 <- left_join(ctt2, temp)


# according to the cut off point given by townsen_calc.R based on 2001 census output area in England and Wales
# -2.945689, -1.588635, 0.2767596, 2.896976 
ctt2$townsend.Q5 <- ifelse(ctt2$townsend.imp<=-2.945689, 1, 
                     ifelse(ctt2$townsend.imp>-2.945689 & ctt2$townsend.imp<=-1.588635, 2, 
                      ifelse(ctt2$townsend.imp>-1.588635 & ctt2$townsend.imp<=0.2767596, 3,
                       ifelse(ctt2$townsend.imp>0.2767596 & ctt2$townsend.imp<=2.896976, 4,5)))) 
ctt2$townsend.Q5 <- as.factor(ctt2$townsend.Q5)
ctt2$townsend.Q5 <- relevel(ctt2$townsend.Q5, ref = "3")

save(ctt2, file = file.path(work_data, "ctt2.Rdata"))
