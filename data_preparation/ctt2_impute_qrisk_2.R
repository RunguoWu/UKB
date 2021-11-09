# new multiple imputation
# mainly for adding HbA1c
rm(list = ls())

library(mice)
library(miceadds)
library(tidyverse)
library(parallel)
library(QRISK3)
library(matrixStats)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

ctt <- readRDS(file.path(work_data, "ctt2020.rds") ) 

ctt2 <- ctt


# multiple imputation -----------------------------------------------------

droplist <- c(1102194, 1330343, 2376313, 2610320, 3039868, 3075648, 3123970,  
              3385704, 3387728, 3462265, 3809403, 4236833, 4602053, 4803604, 
              4890498, 5063590, 5158442)

ctt2 <- ctt2[!ctt2$id %in% droplist, ]

# select vars for imputation

ctt2_imp <- ctt2 %>% select(id, weight, height, LDL, HDL, cholesterol, 
              triglycerides, creatinine_ln, SBP1, SBP2, DBP1, DBP2, hba1c_ln, 
              age.recruit, sex, ethn2, smoke.imp, CVhist, hypertension, BPMed_vi, 
              statin, diabetes.fo.pre.duration5)

  
# simplify auxiliary variables
ctt2_imp$CVhist <- ifelse(ctt2_imp$CVhist==0,0,1)
ctt2_imp$sex <- factor(ctt2_imp$sex, ordered = F)
ctt2_imp$CVhist <- as.factor(ctt2_imp$CVhist)
ctt2_imp$hypertension <- as.factor(ctt2_imp$hypertension)
ctt2_imp$BPMed_vi <- as.factor(ctt2_imp$BPMed_vi)
ctt2_imp$statin <- as.factor(ctt2_imp$statin)
ctt2_imp$smoke.imp <- factor(ctt2_imp$smoke.imp, ordered = F)


# pre-define the matrix ----
imp <- mice(ctt2_imp, maxit=0)
predM = imp$predictorMatrix
meth = imp$method

# ID won't predict other var.
predM[, c("id")]=0

# start impute ----
ptm <- proc.time()

imp <- parlmice(data = ctt2_imp, m=20, maxit = 20, n.core = 5, n.imp.core = 4,
                cluster.seed = 1234,
                predictorMatrix = predM, method = meth)

print(proc.time() - ptm)
print(Sys.time())

imp_list <- miceadds::mids2datlist(imp)

saveRDS(imp, file.path(work_data, "imp.rds" ))
saveRDS(imp_list, file.path(work_data, "imp_list.rds" ))


# after imputation data manipulation --------------------------------------


# all imputations to horizontal form
imp_r <- complete(imp, "broad")

# average of the imputations
ctt2_imp$weight.imp <- rowMeans(imp_r[, grepl("weight", names(imp_r))])
ctt2_imp$height.imp <- rowMeans(imp_r[, grepl("height", names(imp_r))])
ctt2_imp$LDL.imp <- rowMeans(imp_r[, grepl("LDL", names(imp_r))])
ctt2_imp$HDL.imp <- rowMeans(imp_r[, grepl("HDL", names(imp_r))])
ctt2_imp$cholesterol.imp <- rowMeans(imp_r[, grepl("cholesterol", names(imp_r))])
ctt2_imp$triglycerides.imp <- rowMeans(imp_r[, grepl("triglycerides", names(imp_r))])
ctt2_imp$creatinine_ln.imp <- rowMeans(imp_r[, grepl("creatinine_ln", names(imp_r))])
ctt2_imp$SBP1.imp <- rowMeans(imp_r[, grepl("SBP1", names(imp_r))])
ctt2_imp$SBP2.imp <- rowMeans(imp_r[, grepl("SBP2", names(imp_r))])
ctt2_imp$DBP1.imp <- rowMeans(imp_r[, grepl("DBP1", names(imp_r))])
ctt2_imp$DBP2.imp <- rowMeans(imp_r[, grepl("DBP2", names(imp_r))])
ctt2_imp$hba1c_ln.imp <- rowMeans(imp_r[, grepl("hba1c_ln", names(imp_r))])

saveRDS(ctt2_imp, file = file.path(work_data, "ctt2_imp.rds"))

ctt2_imp <- ctt2_imp %>% 
  select(id, weight.imp, height.imp, LDL.imp, HDL.imp, 
         cholesterol.imp, triglycerides.imp, creatinine_ln.imp, SBP1.imp, SBP2.imp, 
         DBP1.imp, DBP2.imp, hba1c_ln.imp)

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

ctt2$BMI_cat.imp <- factor(ctt2$BMI_cat.imp, 
            levels = c("18.5-25", "<18.5", "25-30", "30-35", "35-40", "40+"))


ctt2$creatinine.imp <- exp(ctt2$creatinine_ln.imp)

ctt2$hba1c.imp <- exp(ctt2$hba1c_ln.imp)


# after impute hba1c
# recode people with hba1c>=48 at baseline as diabetes
# diagnosis date as one day before entry
ctt2$diabetes.fo.pre[ctt2$hba1c.imp>=48] <- 1

ctt2$diabetes.fo.pre.date[ctt2$diabetes.fo.pre==1 & is.na(ctt2$diabetes.fo.pre.date)] <- 
  as.Date(ctt2$recruit.date[ctt2$diabetes.fo.pre==1 & is.na(ctt2$diabetes.fo.pre.date)]) - 1

# the post-entry diabetes must be 0 if baseline hba1c.imp>=48
ctt2$diabetes.fo.post[ctt2$hba1c.imp>=48] <- 0

ctt2$diabetes.fo.post.date[ctt2$diabetes.fo.post==0] <- NA

# assume people without diabetes diagnosis and hba1c.imp>=48 as T2DM at beaseline
ctt2$diabetes.T2.fo.pre[ctt2$hba1c.imp>=48] <- 1

ctt2$diabetes.T2.fo.pre.date[ctt2$diabetes.T2.fo.pre==1 & is.na(ctt2$diabetes.T2.fo.pre.date)] <- 
  as.Date(ctt2$recruit.date[ctt2$diabetes.T2.fo.pre==1 & is.na(ctt2$diabetes.T2.fo.pre.date)]) - 1

# so recode post-entry T2DM as 0 for these above
ctt2$diabetes.T2.fo.post[ctt2$hba1c.imp>=48] <- 0

ctt2$diabetes.T2.fo.post.date[ctt2$diabetes.T2.fo.post==0] <- NA

# T1 isn't affected. we treat it as a separate disease as T2


# calculate QRISK ---------------------------------------------------------

d4q <- data.frame(patid=ctt2$id, 
                  gender=ctt2$sex, 
                  age=ctt2$age.recruit, 
                  weight=ctt2$weight.imp, 
                  townsend=ctt2$townsend.imp, 
                  height=ctt2$height.imp, 
                  smoke=ctt2$smoke.qrisk, 
                  ethnicity= ctt2$ethnicity, 
                  blood_pressure_treatment=ctt2$BPMed_vi, 
                  systolic_blood_pressure=ctt2$SBP.imp, 
                  atrial_fibrillation=ctt2$atrial_fibrillation,
                  heart_attack_relative=ctt2$heart_attack_relative,
                  erectile_disfunction=ctt2$erectile_disfunction,
                  rheumatoid_arthritis=ctt2$rheumatoid_arthritis,
                  migraine=ctt2$migraine,
                  systemic_lupus_erythematosis=ctt2$systemic_lupus_erythematosis,
                  severe_mental_illness=ctt2$severe_mental_illness,
                  chronic_kidney_disease=ctt2$chronic_kidney_disease,
                  regular_steroid_tablets=ctt2$regular_steroid_tablets,
                  atypical_antipsy=ctt2$atypical_antipsy, 
                  diabetes1=ctt2$diabetes.T1.fo.pre,
                  diabetes2=ctt2$diabetes.T2.fo.pre
                  )

# recode
d4q$gender <- ifelse(d4q$gender=="Female", 1, 0)
d4q$age <- as.numeric(d4q$age)

# cholesterol/HDL
d4q$cholesterol_HDL_ratio <- ctt2$cholesterol.imp/ctt2$HDL.imp

# sd of SBP
d4q$std_systolic_blood_pressure <- rowSds(as.matrix(ctt2[, c("SBP1.imp", "SBP2.imp")]), na.rm = T)


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

ctt2$qrisk_cat <- factor(ctt2$qrisk_cat, levels = c("<5", "[5,10)", "[10,15)", "[15,20)", ">=20"))



# other stuff -------------------------------------------------------------


# according to the cut off point given by townsen_calc.R based on 2001 census output area in England and Wales
# -2.945689, -1.588635, 0.2767596, 2.896976 
ctt2$townsend.Q5 <- ifelse(ctt2$townsend.imp<=-2.945689, 1, 
                           ifelse(ctt2$townsend.imp>-2.945689 & ctt2$townsend.imp<=-1.588635, 2, 
                                  ifelse(ctt2$townsend.imp>-1.588635 & ctt2$townsend.imp<=0.2767596, 3,
                                         ifelse(ctt2$townsend.imp>0.2767596 & ctt2$townsend.imp<=2.896976, 4,5)))) 
ctt2$townsend.Q5 <- as.factor(ctt2$townsend.Q5)
ctt2$townsend.Q5 <- relevel(ctt2$townsend.Q5, ref = "3")

# add age 10 year group to ctt2
ctt2$age.group10 <- ifelse(ctt2$age.recruit<50, "40-50", 
                           ifelse(ctt2$age.recruit %in% 50:59, "50-59", 
                                  ifelse(ctt2$age.recruit >=60, "60-70",NA)))



# incorporate new diet vars from bd2 ---------------------------------------
bd2 <- readRDS(file.path(work_data, "bd2.rds"))

bd2 <- bd2 %>% select("f.eid", "f.1309.0.0", "f.1319.0.0", "f.1289.0.0", "f.1299.0.0",
              "f.1329.0.0", "f.1339.0.0", "f.1349.0.0", "f.1369.0.0", "f.1379.0.0",
              "f.1389.0.0", "f.1438.0.0", "f.1448.0.0", "f.1458.0.0", "f.1468.0.0")

# recode 7 components of diets

# fruit and veg ----
# 1309 fresh fruit intake
# 1319 dried fruit intake
# 1289 cooked veg intake
# 1299 salad/raw veg intake

# recode fruit and veg variables such as all <0 as NA except -10 as 0.5
codetab <- cbind(c("fruit1", "fruit2", "veg1", "veg2"),
                 c(1309, 1319, 1299, 1289))
for (i in 1:nrow(codetab)) {
  bd2[[codetab[i, 1]]] <- bd2[[paste0("f.", codetab[i, 2], ".0.0")]]
  
  # code -10 less than one portion as 0.5
  bd2[[codetab[i, 1]]][bd2[[codetab[i, 1]]]==-10] <- 0.5
  
  # code other < 0 dont know/not to say as NA
  bd2[[codetab[i, 1]]][bd2[[codetab[i, 1]]]<0] <- NA
  
  bd2[[codetab[i, 1]]] <- as.numeric(bd2[[codetab[i, 1]]])
}

# fruit score >=3 servings/d as 1
bd2$fruit <- bd2$fruit1 + bd2$fruit2
bd2$fruit_score <- NA
bd2$fruit_score[bd2$fruit>=3 | bd2$fruit1>=3 | bd2$fruit2>=3] <- 1
bd2$fruit_score[bd2$fruit<3] <- 0

# veg score >=3 servings/d as 1
bd2$veg <- bd2$veg1 + bd2$veg2
bd2$veg_score <- NA
bd2$veg_score[bd2$veg>=3 | bd2$veg1>=3 | bd2$veg2>=3] <- 1
bd2$veg_score[bd2$veg<3] <- 0

# fish and meat ----
# 1329/1339 oily/non-oily fish intake
# 1349 processed meat intake
# 1369/1379/1389 beef/lamb or mutton/pork intake

# recode fish and meat variables such as all <0 as NA
codetab2 <- cbind(c("fish1", "fish2", "pmeat", "beef", "lamb", "pork"),
                 c(1329, 1339, 1349, 1369, 1379, 1389))
for (i in 1:nrow(codetab2)) {
  
  bd2[[codetab2[i, 1]]] <- bd2[[paste0("f.", codetab2[i, 2], ".0.0")]]
  
  # code other < 0 dont know/not to say as NA
  bd2[[codetab2[i, 1]]][bd2[[codetab2[i, 1]]]<0] <- NA
  
  bd2[[codetab2[i, 1]]] <- as.numeric(bd2[[codetab2[i, 1]]])
}

# fish score >=2 servings/w as 1
bd2$fish_score <- NA
#0 never; 1 = less than once a week; 2 = once a week; 3 = 2-4 times a week
# assume less than once a week as 0.5 a week
bd2$fish_score[bd2$fish1>=3 | bd2$fish2>=3 | bd2$fish1 + bd2$fish2==4] <- 1
bd2$fish_score[bd2$fish1 + bd2$fish2<4 & bd2$fish1<3 & bd2$fish2<3] <- 0

# processed meat <=1 serving/w as 1
bd2$pmeat_score <- ifelse(bd2$pmeat<=2, 1, 0)

# read meat <=1.5 servings/w as 1
bd2$rmeat_score <- NA
bd2$rmeat_score[bd2$beef>=3 | bd2$lamb>=3 | bd2$pork>=3 | 
                bd2$beef + bd2$lamb + bd2$pork >3] <- 0
bd2$rmeat_score[bd2$beef + bd2$lamb + bd2$pork <=3 & 
                  bd2$beef<3 & bd2$lamb<3 & bd2$pork<3] <- 1


# bread and cereal ----
# 1438 bread intake
# 1448 bread type
# 1458 cereal intake
# 1468 cereal type

# recode bread and cereal variables such as all <0 as NA except -10 as 0.5
codetab3 <- cbind(c("bread", "cereal"), 
                  c(1438, 1458))
for (i in 1:nrow(codetab3)) {
  
  bd2[[codetab3[i, 1]]] <- bd2[[paste0("f.", codetab3[i, 2], ".0.0")]]

  # code -10 less than one portion as 0.5
  bd2[[codetab3[i, 1]]][bd2[[codetab3[i, 1]]]==-10] <- 0.5
  
  # code other < 0 dont know/not to say as NA
  bd2[[codetab3[i, 1]]][bd2[[codetab3[i, 1]]]<0] <- NA
  
  bd2[[codetab3[i, 1]]] <- as.numeric(bd2[[codetab3[i, 1]]])
}

bd2$bread_type <- bd2$f.1448.0.0
bd2$cereal_type <- bd2$f.1468.0.0
# code other < 0 dont know/not to say as NA
bd2$bread_type[bd2$bread_type < 0] <- NA
bd2$cereal_type[bd2$cereal_type < 0] <- NA

bd2$bread_type <- as.numeric(bd2$bread_type)
bd2$cereal_type <- as.numeric(bd2$cereal_type)

# bread type
# 1: white
# 2: brown
# 3: wholemeal or wholegrain bread
# 4: other type
# NA dont know/not to answer

# cereal type 
# 1: bran cereal
# 2: biscuit cereal
# 3: oat cereal
# 4: muesli
# 5: other (cornflakes, frosties)
# NA dont know/not to answer

# count bread type 3 and cereal type 1,2,3 as whole grains
# count bread type 1 and cereal type 5 as refined grains

# whole grains >=3 serving/day as 1
# the period of bread and cereal is week
bd2$wholeGrain_score <- NA
bd2$wholeGrain_score[(bd2$bread_type==3 & bd2$bread/7>=3) | 
                     (bd2$cereal_type %in% 1:3 & bd2$cereal/7>=3) |
                     (bd2$bread_type==3 & bd2$cereal_type %in% 1:3 &
                      (bd2$bread/7)+(bd2$cereal/7)>=3)] <- 1

# must be 0 situation, regardless of type
bd2$wholeGrain_score[(bd2$bread/7)+(bd2$cereal/7) < 3] <- 0

# for both no wholegrain bread and cereal eaters
bd2$wholeGrain_score[bd2$bread_type!=3 & bd2$cereal_type %in% 4:5 ] <- 0
# for single wholegrain eater
bd2$wholeGrain_score[bd2$bread_type==3 & bd2$cereal_type %in% 4:5 & 
                       bd2$bread/7 < 3] <- 0
bd2$wholeGrain_score[bd2$bread_type!=3 & bd2$cereal_type %in% 1:3 & 
                       bd2$cereal/7 < 3] <- 0
# for one type missing
bd2$wholeGrain_score[is.na(bd2$bread_type) & bd2$cereal_type %in% 4:5 & 
                       bd2$bread/7 < 3] <- 0
bd2$wholeGrain_score[bd2$bread_type!=3 & is.na(bd2$cereal_type) & 
                       bd2$cereal/7 < 3] <- 0


# refined grains <= 1.5 servings/day as 1
bd2$refinedGrain_score <- NA
bd2$refinedGrain_score[(bd2$bread_type==1 & bd2$bread/7>1.5) | 
                       (bd2$cereal_type ==5 & bd2$cereal/7>1.5) |
                       (bd2$bread_type==1 & bd2$cereal_type ==5 &
                        (bd2$bread/7)+(bd2$cereal/7)>1.5)] <- 0

# must be 1 situation, regardless of type
bd2$refinedGrain_score[(bd2$bread/7)+(bd2$cereal/7)<=1.5] <- 1
# for both no refined eaters
bd2$refinedGrain_score[bd2$bread_type!=1 & bd2$cereal_type !=5] <- 1
# for single refined eater
bd2$refinedGrain_score[bd2$bread_type==1 & bd2$cereal_type !=5 & 
                       bd2$bread/7 <= 1.5] <- 1
bd2$refinedGrain_score[bd2$bread_type!=1 & bd2$cereal_type ==5 & 
                       bd2$cereal/7 <= 1.5] <- 1

# for one type missing
bd2$refinedGrain_score[is.na(bd2$bread_type) & bd2$cereal_type !=5 & 
                       bd2$bread/7 <= 1.5] <- 1
bd2$refinedGrain_score[bd2$bread_type!=1 & is.na(bd2$cereal_type) & 
                       bd2$cereal/7 <= 1.5] <- 1


# extract useful information
diet_score <- bd2 %>% 
  select(fruit_score, veg_score, fish_score, pmeat_score,
         rmeat_score, wholeGrain_score, refinedGrain_score) %>% 
  # calculate sum of score, ignoring NA, but record number of NA
  mutate(diet_tot = rowSums(., na.rm = T), miss_tot = rowSums(is.na(.)))

diet_score$id <- bd2$f.eid
diet_score$healthyDiet <- "uncertain"
diet_score$healthyDiet[diet_score$diet_tot>=4] <- "healthy"
diet_score$healthyDiet[diet_score$diet_tot + diet_score$miss_tot <4] <- "unhealthy" 

# merge with ctt2 ----
ctt2 <- merge(ctt2, diet_score, by.x = "id", by.y = "id", all.x = TRUE)

ctt2$healthyDiet <- "uncertain"
ctt2$healthyDiet[ctt2$diet_tot>=4] <- "healthy"
ctt2$healthyDiet[ctt2$diet_tot + ctt2$miss_tot <4] <- "unhealthy" 


saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"))

