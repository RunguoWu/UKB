library(tidyverse)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

load(file.path(work_data, "bd.Rdata"))

# read inpatient bunk data 
# diagnosis table
hesin_diag <- read.delim(file.path(raw_data, "hesin_diag.txt"))
# use hesin to incorporate dates
hesin <- read.delim(file.path(raw_data, "hesin.txt"))
# operation table
hesin_oper <- read.delim(file.path(raw_data, "hesin_oper.txt"))
# death tables
death <- read.delim(file.path(raw_data, "death.txt"))
death_cause <- read.delim(file.path(raw_data, "death_cause.txt"))

# calculate townsend
osgb <- read.csv(file.path(work_data, "osgb.csv"), header = T)
# for qrisk
corti <- read.csv(file.path(work_data, "cortico.csv"), header = F, encoding = "UTF-8")

# use the main dataset to create a new dataset for recoding variables
# keep the original dataset - "bd" intact
# start with the basic information
ctt <- data.frame(id = bd$f.eid, 
                  birth.year = bd$f.34.0.0, 
                  sex = bd$f.31.0.0,
                  age.recruit = bd$f.21022.0.0, 
                  recruit.date = as.Date(bd$f.53.0.0))
# age group
ctt$age.group <- ifelse(ctt$age.recruit<45, "<45", 
                   ifelse(ctt$age.recruit %in% 45:49, "45-49", 
                     ifelse(ctt$age.recruit %in% 50:54, "50-54",
                       ifelse(ctt$age.recruit %in% 55:59, "55-59",
                         ifelse(ctt$age.recruit %in% 60:64, "60-64",
                           ifelse(ctt$age.recruit>=65, "65+",NA))))))

# bulk data ---------------------------------------------------------------


# the variables such as MI.hes, stroke.hes and CRV.hes derived from bulk HES data 
# are basically the same as those from main dataset
# However, bulk data are updated more frequently than main dataset
# I'd recommend to use those from bulk HES data

# code empty cases of icd 10 and icd 9 as NA
hesin_diag$diag_icd10[hesin_diag$diag_icd10==""] <- NA
hesin_diag$diag_icd9[hesin_diag$diag_icd9==""] <- NA

# transform hesin_diag from long to wide by arr_index
# only keep three id vars. and icd10 and icd9
# converting data from long to wide 
hesin_diag <- hesin_diag %>% 
  select("eid", "ins_index", "arr_index", "diag_icd10", "diag_icd9")%>% 
  group_by(eid, ins_index) %>% 
  gather("diag_icd10", "diag_icd9", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code)

# remove columns with all NAs
hesin_diag <- hesin_diag[, colSums(is.na(hesin_diag))<nrow(hesin_diag)]

# generate diagnosis date, using episode start date as the primary proxy
# admission date as the secondary proxy.
hesin$date <- hesin$epistart
hesin$date[is.na(hesin$date)] <- hesin[is.na(hesin$date), "admidate"]
# convert integer date to date format
# hesin <- transform(hesin, date=as.Date(as.character(date), "%Y%m%d"))
# RW: 2021-03-23: the 2020 updated hesin use different date format
hesin <- transform(hesin, date=as.Date(as.character(date), "%d/%m/%Y"))

hesin <- hesin %>% select("eid", "ins_index", "date", "dsource")

# merge hesin into hesin_diag
hesin_diag <- merge(hesin_diag, hesin, by=c("eid", "ins_index"), all.x = T)

# refer to ctt for recruitment date
hesin_diag <- merge(hesin_diag, ctt[, c("id","recruit.date")], by.x = "eid", 
                    by.y = "id", all.x = T)

# only keep instance where date >= recruitment date
hesin_diag <- hesin_diag[hesin_diag$date>=hesin_diag$recruit.date & !is.na(hesin_diag$date),]

# # censor at 29 Feb 2020 for England and Scotland, and 28 Feb for Wales. RW: 2021-07-12
# # Note, have not implemented this yet
# hesin_diag <- hesin_diag %>% filter((dsource != "PEDW" & date<="2020-02-29") | date<="2018-02-28")

# again, remove columns with all NAs
hesin_diag <- hesin_diag[, colSums(is.na(hesin_diag))<nrow(hesin_diag)]
# all icd-9 codes are effectively removed because they all happened before recruitment

### code HES event ----
# MI first
# follow-up MI, without old MI
hesin_diag$MI.hes <- 0
# MI codes include: all I21X, all I22X, all I23X, I241in ICD-10
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$MI.hes[grepl("^I21|^I22|^I23|I241", hesin_diag[[text2]])] <- 1 # exclude I252, old MI
}

# only keep those with MI=1
# generate a new instance index, because the old one has been discrete due to filtering
# convert from long to wide
hesin_MI <- hesin_diag %>% filter(MI.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_MI$earliest <- apply(hesin_MI[, -1], 1, function(x) min(x, na.rm = T))
hesin_MI <- hesin_MI[, c(1, 57, 2:56)]

# output
hesin_MI <- hesin_MI %>% select("eid", "earliest") %>% mutate(MI.hes=1) %>% 
  rename(MI.hes.date = earliest)
hesin_MI$MI.hes.date <- as.Date(hesin_MI$MI.hes.date)

ctt <- merge(ctt, hesin_MI, by.x = "id", by.y = "eid", all.x = T)
ctt$MI.hes[is.na(ctt$MI.hes)] <- 0
# table(ctt$MI.hes, ctt$MI.inpatient.post)

# stroke
hesin_diag$stroke.hes <- 0
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$stroke.hes[grepl("^I60|^I61|^I62|^I63|I64", hesin_diag[[text2]])] <- 1
}

hesin_stroke <- hesin_diag %>% filter(stroke.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_stroke$stroke.hes.date <- apply(hesin_stroke[, -1], 1, function(x) min(x, na.rm = T))
hesin_stroke$stroke.hes.date <- as.Date(hesin_stroke$stroke.hes.date)
# output
hesin_stroke <- hesin_stroke %>% select("eid", "stroke.hes.date") %>% mutate(stroke.hes=1)

ctt <- merge(ctt, hesin_stroke, by.x = "id", by.y = "eid", all.x = T)
ctt$stroke.hes[is.na(ctt$stroke.hes)] <- 0
# table(ctt$stroke.hes, ctt$stroke.inpatient.post)

# cancer
hesin_diag$cancer.hes <- 0
for (i in 0:19) { 
  text2 <- paste0("diag_icd10_",i)
  hesin_diag$cancer.hes[grepl("^C", hesin_diag[[text2]]) &
                          !grepl("^C44", hesin_diag[[text2]])] <- 1
}

hesin_cancer <- hesin_diag %>% filter(cancer.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_cancer$cancer.hes.date <- apply(hesin_cancer[, -1], 1, function(x) min(x, na.rm = T))
hesin_cancer$cancer.hes.date <- as.Date(hesin_cancer$cancer.hes.date)
# output
hesin_cancer <- hesin_cancer %>% select("eid", "cancer.hes.date") %>% mutate(cancer.hes=1)

ctt <- merge(ctt, hesin_cancer, by.x = "id", by.y = "eid", all.x = T)
ctt$cancer.hes[is.na(ctt$cancer.hes)] <- 0


### apply the same method to operation codes

# code empty cases of oper4 as NA
# oper3 has been coded as NA 
hesin_oper$oper4[hesin_oper$oper4==""] <- NA

hesin_oper <- hesin_oper %>% 
  select("eid", "ins_index", "arr_index", "oper4", "oper3")%>% 
  group_by(eid, ins_index) %>% 
  gather("oper4", "oper3", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% spread(icd_array, code)

# remove columns with all NAs
hesin_oper <- hesin_oper[, colSums(is.na(hesin_oper))<nrow(hesin_oper)]

# merge hesin 
hesin_oper <- merge(hesin_oper, hesin, by=c("eid", "ins_index"), all.x = T)

# refer to ctt for recruitment date
cruit_date <- ctt[, c("id","recruit.date")] 
hesin_oper <- merge(hesin_oper, cruit_date, by.x = "eid", by.y = "id", all.x = T)
rm(cruit_date)

# only keep instance where date >= recruitment date
hesin_oper <- hesin_oper[hesin_oper$date>=hesin_oper$recruit.date,]

# # censor at 29 Feb 2020 for England and Scotland, and 28 Feb for Wales. RW: 2021-07-12
# # Note, have not implemented this yet
# hesin_oper <- hesin_oper %>% filter((dsource != "PEDW" & date<="2020-02-29") | date<="2018-02-28")


# again, remove columns with all NAs
hesin_oper <- hesin_oper[, colSums(is.na(hesin_oper))<nrow(hesin_oper)]
# all opcs 3 codes are effectively removed because they all happened before recruitment

# scan the code row by row to check if there is our interested ones

# CRV
hesin_oper$CRV.hes <- 0
for (i in 0:23) { 
  text2 <- paste0("oper4_",i)
  hesin_oper$CRV.hes[grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", 
                           hesin_oper[[text2]])] <- 1 
}

hesin_CRV <- hesin_oper %>% filter(CRV.hes==1) %>% 
  select("eid", "ins_index", "date") %>% group_by(eid) %>% 
  mutate(ins = 1:n(), ins_index=NULL) %>% spread(ins, date)

# finally select the earliest date among the all
hesin_CRV$CRV.hes.date <- apply(hesin_CRV[, -1], 1, function(x) min(x, na.rm = T))
hesin_CRV$CRV.hes.date <- as.Date(hesin_CRV$CRV.hes.date)
# output
hesin_CRV <- hesin_CRV %>% select("eid", "CRV.hes.date") %>% mutate(CRV.hes=1)

ctt <- merge(ctt, hesin_CRV, by.x = "id", by.y = "eid", all.x = T)
ctt$CRV.hes[is.na(ctt$CRV.hes)] <- 0

# remove hesin data
rm(hesin_diag, hesin, hesin_oper, hesin_MI, hesin_stroke, hesin_CRV, hesin_cancer)


# death based on bulk data --------------------------------------------------

# direcely use bulk data
# bulk deaths include all main dataset deaths and have a few extra

# convert integer date to date format
death <- transform(death, date=as.Date(as.character(date_of_death), "%d/%m/%Y"))
death <- death %>% select(eid, ins_index, date)

death_cause_prim <- death_cause %>% filter(level==1)
death_cause_sec <- death_cause %>% filter(level==2)

# transform from long to wide
death_cause_prim <- death_cause_prim %>% 
  select("eid","ins_index","arr_index","cause_icd10") %>% 
  group_by(eid, ins_index) %>% 
  gather("cause_icd10", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code) 

death_cause_sec <- death_cause_sec %>% 
  select("eid","ins_index","arr_index","cause_icd10") %>% 
  group_by(eid, ins_index) %>% 
  gather("cause_icd10", key = icd, value = code) %>% 
  unite(icd_array, icd, arr_index) %>% 
  spread(icd_array, code) 

# merge
death <- merge(death, death_cause_prim, by=c("eid", "ins_index"), all.x = T)
colnames(death)[4] <- "ICD10_prim"
death <- merge(death, death_cause_sec, by=c("eid", "ins_index"), all.x = T)
colnames(death) <- sub("cause_icd10", "ICD10_sec", colnames(death))

# # censor death at 29 Feb 2020 due to covid. RW: 2021-07-12
# # Note, have not implemented this yet
# death <- death %>% filter(date <= "2020-02-29")


# vascular death
# all ICD-10 I category: circulatory
# all ICD-10 R category: unclassfied elsewhere
death$VD <- 0
death$VD[grepl("^I", death$ICD10_prim) | 
           grepl("^R", death$ICD10_prim)] <- 1
# on the condition of CVD as a secondary cause 
# Y832: anastomosis, bypass or graft
# Y835: amputation of limb
# W19: unspecified fall
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$VD[grepl("^I", death[[text1]]) & (death$ICD10_prim %in% 
                                            c("Y832", "Y835") | grepl("^W19", death$ICD10_prim))] <- 1
}

death_VD <- death %>% filter(VD==1) %>% 
  select(eid, date) %>% 
  distinct(eid, .keep_all = T) %>% 
  mutate(death.vascular=1) %>% 
  rename(death.vascular.date=date)

# merge to ctt
ctt <- merge(ctt, death_VD, by.x = "id", by.y = "eid", all.x = T)
ctt$death.vascular[is.na(ctt$death.vascular)] <- 0

# all death
death_all <- death %>% select(eid, date) %>% 
  distinct(eid, .keep_all = T) %>% 
  mutate(death.allcause=1) %>% 
  rename(death.date=date)

ctt <- merge(ctt, death_all, by.x = "id", by.y = "eid", all.x = T)
ctt$death.allcause[is.na(ctt$death.allcause)] <- 0

# non-vascular death
ctt$death.nonvascular <- ifelse(ctt$death.allcause==1 & 
                                  ctt$death.vascular==0, 1, 0)

ctt$death.nonvascular.date <- ctt$death.date
ctt$death.nonvascular.date[ctt$death.vascular==1] <- NA



########################
##### main dataset #####
########################


# MI ----------------------------------------------------------------------

# baseline MI use UKB algorithm 

# combination of verbal interview and inpatient record before recruitment
ctt$MI.baseline <- ifelse(!is.na(bd$f.42001.0.0) & 
                            bd$f.42000.0.0<ctt$recruit.date, 1, 0)

#### incident MI
# first
# MI death
# death with a primary or secondary reason as MI
# bd$f.40001 has two instance although
# cases in the second instance is actually included in the first instance

# # primary reason as MI
# ctt$MI.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   ctt$MI.death[grepl("^I21|^I22|^I23|I241", bd[[text2]])] <- 1
# }
# 
# # secondary reason as MI
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     ctt$MI.death[grepl("^I21|^I22|^I23|I241", bd[[text2]])] <- 1
#   }
# }
# 
# ctt$MI.death.date <- as.Date(bd$f.40000.0.0)
# ctt$MI.death.date[ctt$MI.death==0] <- NA

# use bulk death data instead
death$MI.death <- 0
# primary reason
death$MI.death[grepl("^I21|^I22|^I23|I241", death$ICD10_prim)] <- 1
# secondary reason
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$MI.death[grepl("^I21|^I22|^I23|I241", death[[text1]])] <- 1
}

MI_death <- death %>% filter(MI.death==1) %>% 
  select(eid, date, MI.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(MI.death.date=date)

# merge to ctt
ctt <- merge(ctt, MI_death, by.x = "id", by.y = "eid", all.x = T)
ctt$MI.death[is.na(ctt$MI.death)] <- 0


# second
# first occurrence date, including primary care
# check I21,22,23 only, since not all I24 belongs to MI
# I21: 131298
# I22: 131300
# I23: 131302
# because it lack some code, just used for complement only

temp <- data.frame(I21=bd$f.131298.0.0, I22=bd$f.131300.0.0, I23=bd$f.131302.0.0)
# generate MI ever
ctt$MI.fo <- 0
ctt$MI.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for MI, we should choose the first one
# before this, check "1901-01-01", "1902-02-02", "1903-03-03" and "2037-07-07". recode them as NA
# summary(temp$I21)
# summary(temp$I22)
# summary(temp$I23)
# answer is NO

earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row names as names
# join using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

ctt$MI.fo.date <- temp$date
ctt$MI.fo.date <- as.Date(ctt$MI.fo.date)

# generate pre-recruitment and post-recruitment
ctt$MI.fo.pre <- ifelse(is.na(ctt$MI.fo.date), 0,
                        ifelse(ctt$MI.fo.date<ctt$recruit.date, 1, 0))
ctt$MI.fo.pre.date <- ctt$MI.fo.date
ctt$MI.fo.pre.date[ctt$MI.fo.pre==0] <- NA

ctt$MI.fo.post <- ifelse(is.na(ctt$MI.fo.date), 0,
                         ifelse(ctt$MI.fo.date>=ctt$recruit.date, 1, 0))
ctt$MI.fo.post.date <- ctt$MI.fo.date
ctt$MI.fo.post.date[ctt$MI.fo.post==0] <- NA


# third
# combine MI death, inpatient MI record and first occurrence post
# date choose whichever the earliest
# data from HES inpatient records contain all inpatient records in the main dataset 
# we'd use HES records instead of that from main dataset inpatient records.  

#ctt$MI.all.post <- ctt$MI.inpatient.post
ctt$MI.all.post <- ctt$MI.hes
ctt$MI.all.post[ctt$MI.death==1] <- 1
ctt$MI.all.post[ctt$MI.fo.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(ctt[ctt$MI.hes==1 | ctt$MI.death==1 | ctt$MI.fo.post==1, 
    c("MI.hes.date", "MI.death.date", "MI.fo.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(ctt$MI.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 
# check if the inpatient.date is always the earliest
# sum(temp$`ctt$MI.inpatient.date`!=temp$date, na.rm = T)

temp$date <- as.Date(temp$date)

# whichever the earliest
ctt$MI.all.date <- temp$date

# finally, also complement MI baseline
ctt$MI.baseline[ctt$MI.fo.pre==1] <- 1



# stroke ------------------------------------------------------------------

# first inpatient record
# We finally use bulk HES data, but don't mute the follows, because baseline stroke 
# refers to the information.
# baseline MI basically follows the UKB algorithm so do not need inpatient data

ctt$stroke.inpatient <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  ctt$stroke.inpatient[grepl("^I60|^I61|^I62|^I63|I64 ", bd[[text2]])] <- 1
}
# ICD-9 codes all happened years before recruitment.
# they are included here only to supplement identifying baseline stroke
# so they are not considered in identifying dates
for (i in 0:46) { # f.41271 has 47 arrays
  text2 <- paste0("f.41271.0.",i)
  ctt$stroke.inpatient[bd[[text2]]=="4309" |
                         bd[[text2]]=="4319" |
                         bd[[text2]]=="4349" |
                         bd[[text2]]=="4369" |
                         bd[[text2]]=="4320" |
                         bd[[text2]]=="4321"] <- 1

}
# diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to stroke diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(ctt$stroke.inpatient)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("stroke.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]]) &
                  !is.na(bd[[text2]])] <-
    bd[[text3]][grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]]) &
                  !is.na(bd[[text2]])]
}

# if there are multiple date for stroke ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:212){
  text1 <- paste0("stroke.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < ctt$recruit.date] <- NA
}
temp$`ctt$stroke.inpatient` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) <
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T)
# give the value to ctt
ctt$stroke.inpatient.date <- temp$date
ctt$stroke.inpatient.date <- as.Date(ctt$stroke.inpatient.date)

# stroke.inpatient.date is only post-recruitment

# dates linked to ICD-9 codes are no later than 1996-03-30
# don't work on ICD-9 dates

ctt$stroke.inpatient.post <- ctt$stroke.inpatient
ctt$stroke.inpatient.post[is.na(ctt$stroke.inpatient.date)] <- 0

# baseline stroke
# base on f.20002, self-reported non-cancer illness
text1 <- "stroke.baseline"
ctt[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  ctt[[text1]][bd[[text2]]== 1081 |
                 bd[[text2]]== 1583 |
                 bd[[text2]]== 1086 |
                 bd[[text2]]== 1491 |
                 bd[[text2]]== 1083] <- 1
}

# and inpatient record before recruitment
# is.na(ctt$stroke.inpatient.date) means the date precede recruit date
ctt$stroke.baseline[ctt$stroke.inpatient==1 & is.na(ctt$stroke.inpatient.date)] <- 1

# second
# stroke death
# death with a primary or secondary reason as stroke
# bd$f.40001 has two instance although, 
# cases in the second instance is actually included in the first instance
# primary reason as stroke
# ctt$stroke.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   ctt$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]])] <- 1
# }
# 
# # secondary reason as stroke
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     ctt$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", bd[[text2]])] <- 1
#   }
# }
# 
# ctt$stroke.death.date <- bd$f.40000.0.0
# ctt$stroke.death.date[ctt$stroke.death==0] <- NA

# use bulk death data instead
death$stroke.death <- 0
# primary reason
death$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", death$ICD10_prim)] <- 1
# secondary reason
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$stroke.death[grepl("^I60|^I61|^I62|^I63|I64", death[[text1]])] <- 1
}

stroke_death <- death %>% filter(stroke.death==1) %>% 
  select(eid, date, stroke.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(stroke.death.date=date)

# merge to ctt
ctt <- merge(ctt, stroke_death, by.x = "id", by.y = "eid", all.x = T)
ctt$stroke.death[is.na(ctt$stroke.death)] <- 0


# third
# first occurrence date, including primary care
# check I21,22,23 only, since not all I24 belongs to MI
# I60: 131360
# I61: 131362
# I62: 131364
# I63: 131366
# I64: 131368

temp <- data.frame(I60=bd$f.131360.0.0, I61=bd$f.131362.0.0, I62=bd$f.131364.0.0, 
                  I63=bd$f.131366.0.0, I64=bd$f.131368.0.0)
# generate MI ever
ctt$stroke.fo <- 0
ctt$stroke.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for stroke, we should choose the first one
# before this, check "1901-01-01", "1902-02-02", "1903-03-03" and "2037-07-07". recode them as NA
# summary(temp$I60)
# summary(temp$I61)
# summary(temp$I62)
# summary(temp$I63)
# summary(temp$I64)
# answer is NO

earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row names as names
# join using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

ctt$stroke.fo.date <- temp$date
ctt$stroke.fo.date <- as.Date(ctt$stroke.fo.date)

# generate pre-recruitment and post-recruitment
ctt$stroke.fo.pre <- ifelse(is.na(ctt$stroke.fo.date), 0,
                            ifelse(ctt$stroke.fo.date<ctt$recruit.date, 1, 0))
ctt$stroke.fo.pre.date <- ctt$stroke.fo.date
ctt$stroke.fo.pre.date[ctt$stroke.fo.pre==0] <- NA

ctt$stroke.fo.post <- ifelse(is.na(ctt$stroke.fo.date), 0,
                             ifelse(ctt$stroke.fo.date>=ctt$recruit.date, 1, 0))
ctt$stroke.fo.post.date <- ctt$stroke.fo.date
ctt$stroke.fo.post.date[ctt$stroke.fo.post==0] <- NA

# forth
# combine stroke death, inpatient stroke record and first occurrence post
# date choose whichever the earliest

# ctt$stroke.all.post <- ctt$stroke.inpatient.post
ctt$stroke.all.post <- ctt$stroke.hes
ctt$stroke.all.post[ctt$stroke.death==1] <- 1
ctt$stroke.all.post[ctt$stroke.fo.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(ctt[ctt$stroke.hes==1 | ctt$stroke.death==1 | ctt$stroke.fo.post==1, 
   c("stroke.hes.date", "stroke.death.date", "stroke.fo.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(ctt$stroke.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

# whichever the earliest
ctt$stroke.all.date <- temp$date

# finally also complement stroke baseline
ctt$stroke.baseline[ctt$stroke.fo.pre==1] <- 1


# CRV ---------------------------------------------------------------------

# only use hospital diagnoses
# no baseline CRV

ctt$CRV.inpatient <- 0
# MI codes include: all K49X, K501, all K75X, all K76X, all K40X
# all K41X, all K42X, all K43X, all K44X, all K45X, all K46X in OPCS-4
for (i in 0:116) { # f.41272 has 117 arrays
  text2 <- paste0("f.41272.0.",i)
  ctt$CRV.inpatient[grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]])] <- 1
}

# CRV date
temp <- as.data.frame(ctt$CRV.inpatient)
# OPCS-4 only
# dates linked to OPCS-3 codes are no later than 1989, so don't use them
for (i in 0:116) { # f.41282 has 117 arrays, corresponding to 41272
  text1 <- paste0("CRV.inpatient.date.", i)
  text2 <- paste0("f.41272.0.",i)
  text3 <- paste0("f.41282.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^K49|^K75|^K76|^K40|^K41|^K42|^K43|^K44|^K45|^K46|K501|K504", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

for (i in 0:116){
  text1 <- paste0("CRV.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < ctt$recruit.date] <- NA
}
temp$`ctt$CRV.inpatient` <- NULL
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$CRV.inpatient.date <- temp$date
ctt$CRV.inpatient.date <- as.Date(ctt$CRV.inpatient.date)

# post-recruitment 
ctt$CRV.inpatient.post <- ctt$CRV.inpatient
ctt$CRV.inpatient.post[is.na(ctt$CRV.inpatient.date)] <- 0

# CRV.hes include all CRV.inpatient.post
ctt$CRV.all.post <- ctt$CRV.hes

# so just pick up the earliest date among the two
earliest <- apply(ctt[ctt$CRV.inpatient.post==1 | ctt$CRV.hes, 
      c("CRV.inpatient.date", "CRV.hes.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(ctt$CRV.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

ctt$CRV.all.date <- temp$date

# check 
# identical(ctt$CRV.all.date, ctt$CRV.hes.date)
# the two are the same, so we actually only need to use CRV.hes data


# cancer ------------------------------------------------------------------

# first
# hospital diagnoses

ctt$cancer.inpatient <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:212) { # f.41270 has 213 arrays
  text2 <- paste0("f.41270.0.",i)
  # all malignant neoplasms but C44 family of non-melanoma skin cancers
  ctt$cancer.inpatient[grepl("^C", bd[[text2]]) &
                         !grepl("^C44", bd[[text2]])] <- 1
}
# 140X-208X in ICD-9
# ICD-9 codes all happened years before recruitment. 
# they are included here only to supplement identifying baseline cancer 
# so they are not considered in identifying dates
for (i in 0:46) { # f.41271 has 47 arrays
  text2 <- paste0("f.41271.0.",i)
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  ctt$cancer.inpatient[(grepl("^14", bd[[text2]]) | 
                          grepl("^15", bd[[text2]]) |
                          grepl("^16", bd[[text2]]) |
                          grepl("^17", bd[[text2]]) |
                          grepl("^18", bd[[text2]]) |
                          grepl("^19", bd[[text2]]) |
                          grepl("^20", bd[[text2]])) &
                         !grepl("^173", bd[[text2]])] <- 1
}
# cancer diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to cancer diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(ctt$cancer.inpatient)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("cancer.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:212){
  text1 <- paste0("cancer.inpatient.date.", i)
  temp[[text1]][temp[[text1]] < ctt$recruit.date] <- NA
}
temp$`ctt$cancer.inpatient` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$cancer.inpatient.post.date <- temp$date
ctt$cancer.inpatient.post.date <- as.Date(ctt$cancer.inpatient.post.date)

ctt$cancer.inpatient.post <- ctt$cancer.inpatient
ctt$cancer.inpatient.post[is.na(ctt$cancer.inpatient.post.date)] <- 0


# inpatient records date before recruitment

# ICD-10
temp <- as.data.frame(ctt$cancer.inpatient)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("cancer.inpatient.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

for (i in 0:46) { # f.41271 has 47 arrays
  text1 <- paste0("cancer.inpatient.date.", i+213)
  text2 <- paste0("f.41271.0.",i)
  text3 <- paste0("f.41281.0.",i)
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

# if there are multiple date for cancer ICD, use the earliest one before recruitment date
# get rid of date after recruitment date and first column
for (i in 0:259){
  text1 <- paste0("cancer.inpatient.date.", i)
  temp[[text1]][temp[[text1]] >= ctt$recruit.date] <- NA
}
temp$`ctt$cancer.inpatient` <- NULL

# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$cancer.inpatient.pre.date <- temp$date
ctt$cancer.inpatient.pre.date <- as.Date(ctt$cancer.inpatient.pre.date)

# there is a date, there is a baseline cancer
ctt$cancer.inpatient.pre <- ctt$cancer.inpatient
ctt$cancer.inpatient.pre[is.na(ctt$cancer.inpatient.pre.date)] <- 0

# second
# cancer registry

ctt$cancer.registry <- 0
# MI codes include: all CXXX in ICD-10
for (i in 0:16) { # f.40006: type of cancer ICD-10 has 17 instances 
  text2 <- paste0("f.40006.",i,".0")
  # all malignant neoplasms but C44 family of non-melanoma skin cancers
  ctt$cancer.registry[grepl("^C", bd[[text2]]) &
                        !grepl("^C44", bd[[text2]])] <- 1
}
# 140X-208X in ICD-9
# ICD-9 codes all happened years before recruitment. 
# they are included here only to supplement identifying baseline cancer 
# so they are not considered in identifying dates
for (i in 0:14) { # f.40013: type of cancer ICD-9 has 15 instances
  text2 <- paste0("f.40013.",i, ".0")
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  ctt$cancer.registry[(grepl("^14", bd[[text2]]) | 
                         grepl("^15", bd[[text2]]) |
                         grepl("^16", bd[[text2]]) |
                         grepl("^17", bd[[text2]]) |
                         grepl("^18", bd[[text2]]) |
                         grepl("^19", bd[[text2]]) |
                         grepl("^20", bd[[text2]])) &
                        !grepl("^173", bd[[text2]])] <- 1
}
# cancer diagnosis date
# ICD-10
# create a new date frame to collect all dates linked to cancer diagnoses
# one individual may have more than one diagnoses and linked dates
temp <- as.data.frame(ctt$cancer.registry)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:16) { # f.40005: date of cancer diagnosis has 213 instances
  text1 <- paste0("cancer.registry.date.", i)
  text2 <- paste0("f.40006.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date before recruitment date and first column
for (i in 0:16){
  text1 <- paste0("cancer.registry.date.", i)
  temp[[text1]][temp[[text1]] < ctt$recruit.date] <- NA
}
temp$`ctt$cancer.registry` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$cancer.registry.post.date <- temp$date
ctt$cancer.registry.post.date <- as.Date(ctt$cancer.registry.post.date)

ctt$cancer.registry.post <- ctt$cancer.registry
ctt$cancer.registry.post[is.na(ctt$cancer.registry.post.date)] <- 0


# cancer register date before recruitment

# ICD-10

temp <- as.data.frame(ctt$cancer.registry)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:16) { # f.40005: date of cancer diagnosis has 17 instances
  text1 <- paste0("cancer.registry.date.", i)
  text2 <- paste0("f.40006.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^C", bd[[text2]]) &
                  !grepl("^C44", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# ICD-9
for (i in 0:14) { # f.40013 has 15 arrays
  text1 <- paste0("cancer.registry.date.", i+17)
  text2 <- paste0("f.40013.",i,".0")
  text3 <- paste0("f.40005.",i,".0")
  # all malignant neoplasms but 173 family of non-melanoma skin cancers
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^14", bd[[text2]]) | 
                   grepl("^15", bd[[text2]]) |
                   grepl("^16", bd[[text2]]) |
                   grepl("^17", bd[[text2]]) |
                   grepl("^18", bd[[text2]]) |
                   grepl("^19", bd[[text2]]) |
                   grepl("^20", bd[[text2]])) &
                  !grepl("^173", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}
# if there are multiple date for cancer ICD, use the earliest one after recruitment date
# get rid of date after recruitment date and first column
for (i in 0:31){
  text1 <- paste0("cancer.registry.date.", i)
  temp[[text1]][temp[[text1]] >= ctt$recruit.date] <- NA
}
temp$`ctt$cancer.registry` <- NULL

# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$cancer.registry.pre.date <- temp$date
ctt$cancer.registry.pre.date <- as.Date(ctt$cancer.registry.pre.date)

# cancer registry before recruitment
ctt$cancer.registry.pre <- ctt$cancer.registry
ctt$cancer.registry.pre[is.na(ctt$cancer.registry.pre.date)] <- 0


# third
# baseline cancer
# base on f.20001, self-reported cancer
# exclude non-melanoma skin cancers

text1 <- "cancer.baseline"
ctt[[text1]] <- 0
for (j in 0:5) {
  text2 <- paste0("f.20001.0.", j)
  ctt[[text1]][!is.na(bd[[text2]]) & 
                 bd[[text2]] != 1060 &
                 bd[[text2]] != 1061 &
                 bd[[text2]] != 1062 &
                 bd[[text2]] != 1073] <- 1
}

# cancer baseline date 
temp <- as.data.frame(ctt$cancer.baseline)
for (j in 0:5) {
  text1 <- paste0("cancer.baseline.date.", j)
  text2 <- paste0("f.20001.0.", j)
  text3 <- paste0("f.20006.0.", j)
  temp[[text1]] <- NA
  temp[[text1]][!is.na(bd[[text2]]) & 
                  bd[[text2]] != 1060 &
                  bd[[text2]] != 1061 &
                  bd[[text2]] != 1062 &
                  bd[[text2]] != 1073] <- 
    bd[[text3]][!is.na(bd[[text2]]) & 
                  bd[[text2]] != 1060 &
                  bd[[text2]] != 1061 &
                  bd[[text2]] != 1062 &
                  bd[[text2]] != 1073]
}
# select the earliest
temp$`ctt$cancer.baseline` <- NULL
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# select the earliest date
# the following command only apply to rows with at least one non-missing values
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row
# names as names
# join to temp using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 
# give the value to ctt
ctt$cancer.baseline.date <- temp$date
# recode date=-1 NA
ctt$cancer.baseline.date[ctt$cancer.baseline.date==-1] <- NA

# RW 2020-12-17
# because baseline date is only year, assume it is the mid day of the year
# the computer interestingly assume the date as today, though year is different
# so we minus the difference between the current date XXXX-XX-XX and the XXXX-01-01
# then add 182 to the mid of the year
dif <- as.numeric(Sys.Date() - as.Date(paste0(substr(as.character(Sys.Date()), 1, 4), "-01-01")))
ctt$cancer.baseline.date <- as.Date(as.character(ctt$cancer.baseline.date), format = "%Y")-dif+182

# a problem is that some baseline cancer dates exceed the recruitment date (n =374) 
# we have to assume them 1 day before the recruitment date
ctt$cancer.baseline.date[!is.na(ctt$cancer.baseline.date) & ctt$cancer.baseline.date > ctt$recruit.date] <- 
  ctt$recruit.date[!is.na(ctt$cancer.baseline.date) & ctt$cancer.baseline.date > ctt$recruit.date]-1

# and inpatient record before recruitment
# is.na(ctt$cancer.inpatient.post.date) means the date precede recruit date
ctt$cancer.baseline.all <- ctt$cancer.baseline
# inpatient
ctt$cancer.baseline.all[ctt$cancer.inpatient.pre==1] <- 1
# cancer register
ctt$cancer.baseline.all[ctt$cancer.registry.pre==1] <- 1

# fourth
# cancer death
# death with a primary or secondary reason as cancer
# bd$f.40001 has two instance although, 
# cases in the second instance is actually included in the first instance
# primary reason as cancer
# ctt$cancer.death <- 0 # RW 2020-10-12
# for (i in 0:1) {
#   text2 <- paste0("f.40001.",i, ".", 0)
#   ctt$cancer.death[grepl("^C", bd[[text2]]) &
#                      !grepl("^C44", bd[[text2]])] <- 1
# }
# # secondary reason as cancer
# for (i in 0:1) {
#   for (j in 0: 13) {
#     text2 <- paste0("f.40002.",i, ".", j)
#     ctt$cancer.death[grepl("^C", bd[[text2]]) &
#                        !grepl("^C44", bd[[text2]])] <- 1
#     
#   }
# }
# 
# ctt$cancer.death.date <- as.Date(bd$f.40000.0.0)
# ctt$cancer.death.date[ctt$cancer.death==0] <- NA

# use bulk death data

death$cancer.death <- 0
# primary reason as cancer
death$cancer.death[grepl("^C", death$ICD10_prim) & 
           !grepl("^C44", death$ICD10_prim)] <- 1
# secondary reason as cancer
for (i in 1:14) {
  text1 <- paste0("ICD10_sec_",i)
  death$cancer.death[grepl("^C", death[[text1]]) & !grepl("^C44", death[[text1]])] <- 1
}

cancer_death <- death %>% filter(cancer.death==1) %>% 
  select(eid, date, cancer.death) %>% 
  distinct(eid, .keep_all = T) %>% 
  rename(cancer.death.date=date)

# merge to ctt
ctt <- merge(ctt, cancer_death, by.x = "id", by.y = "eid", all.x = T)
ctt$cancer.death[is.na(ctt$cancer.death)] <- 0


# fifth
# combine cancer death, inpatient cancer record and first occurrence post
# date choose whichever the earliest
# as MI and stroke, cancer.hes replaces cancer.inpatient.post

# ctt$cancer.all.post <- ctt$cancer.inpatient.post
ctt$cancer.all.post <- ctt$cancer.hes
ctt$cancer.all.post[ctt$cancer.death==1] <- 1
ctt$cancer.all.post[ctt$cancer.registry.post==1] <- 1

# pick up the earliest date among the three
earliest <- apply(ctt[ctt$cancer.hes==1 | ctt$cancer.death==1 | 
          ctt$cancer.registry.post==1, c("cancer.hes.date", "cancer.death.date", 
                "cancer.registry.post.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(ctt$cancer.hes.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

temp$date <- as.Date(temp$date)

# whichever the earliest
ctt$cancer.all.date <- temp$date

# baseline cancer diagnosis date
# interview, inpatient, register
# missing values exist
# pick up the earliest date among the three
earliest <- apply(ctt[!is.na(ctt$cancer.inpatient.pre.date) | 
        !is.na(ctt$cancer.baseline.date) | !is.na(ctt$cancer.registry.pre.date),
        c("cancer.inpatient.pre.date", "cancer.baseline.date", 
          "cancer.registry.pre.date")], 1, function(x) min(x, na.rm = T))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))

temp <- as.data.frame(ctt$cancer.inpatient.post.date)
temp$key <- as.numeric(rownames(temp))
temp <- merge(temp, earliest, all = T) 

# whichever the earliest
ctt$cancer.baseline.all.date <- as.Date(temp$date)
# there are a few missing dates 
# because there are 77 missing dates for baseline interview cancer data.
# assume them to be the median of the cancer.baseline.date - "2002-07-01"
ctt$cancer.baseline.all.date[is.na(ctt$cancer.baseline.all.date) & ctt$cancer.baseline.all==1] <- 
  as.Date(as.numeric(summary(ctt$cancer.baseline.date)["Median"]), origin= "1970-01-01")

# incident cancer that has no baseline cancer

ctt$cancer.incident.only <- ctt$cancer.all.post
ctt$cancer.incident.only[ctt$cancer.baseline.all==1] <- 0
ctt$cancer.incident.only.date <- ctt$cancer.all.date
ctt$cancer.incident.only.date[ctt$cancer.baseline.all==1] <- NA


# diabetes ----------------------------------------------------------------

###
# incorporate diabetes medication info from primary care data 
# there are four variables that record prescriptions:
# read 2 code, BNF code, dmd code and drug name
# All have considerable missing values
# The best combination is read 2 + drug name, with missingness for only five ids, for whom other code variables are also missing

# insulin prescription
# read2
# f1... short-acting insulin preparations
# f2... medium/long-acting insulins
# fw... short with intermediate-acting insulins
gp_scr <- read.delim(file.path(raw_data, "gp_scripts.txt"))
read2 <- c("^f1", "^f2", "^fw")
rd2ptn <- paste(read2, collapse = "|")

# drug name
drug_name <- c("insulin","hypurin","neusulin","quicksol","velosulin","actrapid",
               "humulin","novopen","penject","humaject","pur-in","autopen","BD pen",
               "BD ultra pen","insuman","diapen","exubera","humalog","novorapid",
               "apidra","rapitard","penmix","lentard","neulente","tempulin",
               "monotard","semitard","ultratard","insulatard","monophane","neuphane",
               "initard","mixtard","protaphane","actraphane","protaphane","isophane",
               "lantus","toujeo","abasaglar","levemir","tresiba","xultophy", "novomix")

drugptn <- paste(drug_name, collapse = "|")

# filter to leave those with insulin prescription

temp <- gp_scr %>% filter(grepl(rd2ptn, read_2) | grepl(drugptn, drug_name, ignore.case = T))

# transform to date format
temp <- transform(temp, insulin_date=as.Date(as.character(issue_date), "%d/%m/%Y"))

# keep the earliest insulin use records for each id

temp <- temp %>% group_by(eid) %>% filter(insulin_date==min(insulin_date)) %>% select(eid, insulin_date) %>% distinct(eid, .keep_all = T)

# merge to ctt
ctt <- merge(ctt, temp, by.x = "id", by.y = "eid", all.x = T)

# non metformin anti-diabetic drugs
# the same method as above
read2 <- c("^f3", "^f5", "^f6", "^f8", "^ft")
# drug name: import from mannually generated list
drug_name <- read.csv(file.path(work_data, "non-metformin.csv"), header = F)
drug_name <- drug_name$V1

rd2ptn <- paste(read2, collapse = "|")
drugptn <- paste(drug_name, collapse = "|")

temp <- gp_scr %>% filter(grepl(rd2ptn, read_2) | grepl(drugptn, drug_name, ignore.case = T))
temp <- transform(temp, nonmetformin_date=as.Date(as.character(issue_date), "%d/%m/%Y"))
temp <- temp %>% group_by(eid) %>% filter(nonmetformin_date==min(nonmetformin_date)) %>% select(eid, nonmetformin_date) %>% distinct(eid, .keep_all = T)

ctt <- merge(ctt, temp, by.x = "id", by.y = "eid", all.x = T)

###


# use first occurrence data only
# all inpatient and death diabetes are included in first occurrence
# a very few baseline diabetes are not included 
# check by hand, those missing in first occurrence but available in baseline interview is because the date for the interview data is -1, representing "uncertain or unknown"
# for them, we cannot calculate duration, so will also abandon them

# use first occurrence data 

# five date.field for diabetes
# E10: bd$f.130706.0.0
# E11: bd$f.130708.0.0
# E12: bd$f.130710.0.0 # exclude
# E13: bd$f.130712.0.0
# E14: bd$f.130714.0.0

temp <- data.frame(E10=bd$f.130706.0.0, E11=bd$f.130708.0.0, 
                   E13=bd$f.130712.0.0, E14=bd$f.130714.0.0)

ctt$diabetes.fo <- 0
ctt$diabetes.fo[rowSums(!is.na(temp)) > 0] <- 1
# given one person can have more than one dates for diabetes, we should choose the first one
# before this, bd$f.130706.0.0 and bd$f.130708.0.0 have 1 and 3 "1902-02-02" respectively, 
# which means code has event date matching date of birth. recode them as NA
temp$E10[temp$E10=="1902-02-02"] <- NA
temp$E11[temp$E11=="1902-02-02"] <- NA
# after recoding, sum(rowSums(!is.na(temp)) > 0) has unchanged number
# it means no diabetes record will lose date
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row names as names
# join to diabetes.date using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

ctt$diabetes.fo.date <- temp$date
ctt$diabetes.fo.date <- as.Date(ctt$diabetes.fo.date)

# so far we have diabetes record and linked date as the first diagnosis
ctt$diabetes.fo.pre <- ifelse(is.na(ctt$diabetes.fo.date), 0,
                              ifelse(ctt$diabetes.fo.date<ctt$recruit.date, 1, 0))
ctt$diabetes.fo.pre.date <- ctt$diabetes.fo.date
ctt$diabetes.fo.pre.date[ctt$diabetes.fo.pre==0] <- NA

ctt$diabetes.fo.post <- ifelse(is.na(ctt$diabetes.fo.date), 0,
                               ifelse(ctt$diabetes.fo.date>=ctt$recruit.date, 1, 0))
ctt$diabetes.fo.post.date <- ctt$diabetes.fo.date
ctt$diabetes.fo.post.date[ctt$diabetes.fo.post==0] <- NA

# # duration of baseline diabetes until recruitment
# ctt$diabetes.fo.pre.duration <- ctt$recruit.date-ctt$diabetes.fo.pre.date
# ctt$diabetes.fo.pre.duration <- as.numeric(ctt$diabetes.fo.pre.duration)

######## add medication-identified diabetes
# RW 2021-01-07

# pre

# keep the earlier from the above two prescription date
ctt <- transform(ctt, medidate_pre = pmin(insulin_date, nonmetformin_date, na.rm = T))
# T1, T2 could use this later as well.
# add medication-identified diabetes
ctt$diabetes.fo.pre[ctt$diabetes.fo.pre==0 & 
                      !is.na(ctt$medidate_pre) & ctt$medidate_pre<ctt$recruit.date] <- 1                  

# use the one above as the date of diagnosis  
ctt$diabetes.fo.pre.date[is.na(ctt$diabetes.fo.pre.date) & ctt$diabetes.fo.pre==1] <- 
  ctt$medidate_pre[is.na(ctt$diabetes.fo.pre.date) & ctt$diabetes.fo.pre==1]

# post
ctt$diabetes.fo.post[ctt$diabetes.fo.post==0 &
    ((!is.na(ctt$insulin_date) & ctt$insulin_date>=ctt$recruit.date) | 
    (!is.na(ctt$nonmetformin_date) & ctt$nonmetformin_date>=ctt$recruit.date))] <- 1  

ctt$diabetes.fo.post[ctt$diabetes.fo.pre==1] <- 0

# keep insulin issue date after recruitment 
ctt$medidate1 <- ctt$insulin_date
ctt$medidate1[!is.na(ctt$insulin_date) & ctt$insulin_date < ctt$recruit.date] <- NA
# keep nonmetformin diabetic drug issue date after recruitment
ctt$medidate2 <- ctt$nonmetformin_date
ctt$medidate2[!is.na(ctt$nonmetformin_date) & ctt$nonmetformin_date < ctt$recruit.date] <- NA
# keep the earlier from the above two
ctt <- transform(ctt, medidate = pmin(medidate1, medidate2, na.rm = T))
# T1, T2 could use this later as well.

# use the one above as the date of diagnosis  
ctt$diabetes.fo.post.date[is.na(ctt$diabetes.fo.post.date) & ctt$diabetes.fo.post==1] <- 
  ctt$medidate[is.na(ctt$diabetes.fo.post.date) & ctt$diabetes.fo.post==1]

ctt$diabetes.fo.post.date[ctt$diabetes.fo.post==0] <- NA

###################################################
# split T1 T2 diabetes
# leave E13 and E14 unspecified diabetes at the end  

# T1 first
# E10: bd$f.130706.0.0
ctt$diabetes.T1.fo.date <- bd$f.130706.0.0
ctt$diabetes.T1.fo.date[ctt$diabetes.T1.fo.date=="1902-02-02"] <- NA

ctt$diabetes.T1.fo.date <- as.Date(ctt$diabetes.T1.fo.date)
ctt$diabetes.T1.fo <- ifelse(is.na(ctt$diabetes.T1.fo.date), 0 , 1)

# then T2

# E11: bd$f.130708.0.0
ctt$diabetes.T2.fo.date <- bd$f.130708.0.0
ctt$diabetes.T2.fo.date[ctt$diabetes.T2.fo.date=="1902-02-02"] <- NA

ctt$diabetes.T2.fo.date <- as.Date(ctt$diabetes.T2.fo.date)
ctt$diabetes.T2.fo <- ifelse(is.na(ctt$diabetes.T2.fo.date), 0 , 1)

# other/unspecified
# E13: bd$f.130712.0.0
# E14: bd$f.130714.0.0

temp <- data.frame(E13=bd$f.130712.0.0, E14=bd$f.130714.0.0)

ctt$diabetes.ukn.fo <- 0
ctt$diabetes.ukn.fo[rowSums(!is.na(temp)) > 0] <- 1
earliest <- apply(temp[rowSums(is.na(temp)) < 
                         ncol(temp), ], 1, function(x) min(x, na.rm = T))
# earliest is a array with fewer elements than nrow(ctt), but keeps the same row names as names
# join to diabetes.date using rownames as the key
temp$key <- as.numeric(rownames(temp))
earliest <- data.frame(date=earliest, key=as.numeric(names(earliest)))
temp <- merge(temp, earliest, all = T) 

ctt$diabetes.ukn.fo.date <- temp$date
ctt$diabetes.ukn.fo.date <- as.Date(ctt$diabetes.ukn.fo.date)

# for those with both T1/T2 and unknown code
# replace unknown with T1/T2, but use the first diagnosis date
ctt$diabetes.T1.fo.date[ctt$diabetes.T1.fo==1 & ctt$diabetes.ukn.fo==1 & 
              ctt$diabetes.ukn.fo.date < ctt$diabetes.T1.fo.date] <-
  ctt$diabetes.ukn.fo.date[ctt$diabetes.T1.fo==1 & ctt$diabetes.ukn.fo==1 & 
                             ctt$diabetes.ukn.fo.date < ctt$diabetes.T1.fo.date]

ctt$diabetes.T2.fo.date[ctt$diabetes.T2.fo==1 & ctt$diabetes.ukn.fo==1 & 
                          ctt$diabetes.ukn.fo.date < ctt$diabetes.T2.fo.date] <-
  ctt$diabetes.ukn.fo.date[ctt$diabetes.T2.fo==1 & ctt$diabetes.ukn.fo==1 & 
                             ctt$diabetes.ukn.fo.date < ctt$diabetes.T2.fo.date]

ctt$diabetes.ukn.fo[ctt$diabetes.T1.fo==1 | ctt$diabetes.T2.fo==1] <- 0
ctt$diabetes.ukn.fo.date[ctt$diabetes.ukn.fo==0] <- NA


# finally divide T1, T2 to pre and post recruitment
# pre
# T1
ctt$diabetes.T1.fo.pre <- ifelse(is.na(ctt$diabetes.T1.fo.date), 0,
                                 ifelse(ctt$diabetes.T1.fo.date<ctt$recruit.date, 1, 0))
ctt$diabetes.T1.fo.pre.date <- ctt$diabetes.T1.fo.date
ctt$diabetes.T1.fo.pre.date[ctt$diabetes.T1.fo.pre==0] <- NA
# # duration of baseline diabetes until recruitment
# ctt$diabetes.T1.fo.pre.duration <- ctt$recruit.date-ctt$diabetes.T1.fo.pre.date
# ctt$diabetes.T1.fo.pre.duration <- as.numeric(ctt$diabetes.T1.fo.pre.duration)

# T2
ctt$diabetes.T2.fo.pre <- ifelse(is.na(ctt$diabetes.T2.fo.date), 0,
                                 ifelse(ctt$diabetes.T2.fo.date<ctt$recruit.date, 1, 0))
ctt$diabetes.T2.fo.pre.date <- ctt$diabetes.T2.fo.date
ctt$diabetes.T2.fo.pre.date[ctt$diabetes.T2.fo.pre==0] <- NA
# # duration of baseline diabetes until recruitment
# ctt$diabetes.T2.fo.pre.duration <- ctt$recruit.date-ctt$diabetes.T2.fo.pre.date
# ctt$diabetes.T2.fo.pre.duration <- as.numeric(ctt$diabetes.T2.fo.pre.duration)

# unkown
ctt$diabetes.ukn.fo.pre <- ifelse(is.na(ctt$diabetes.ukn.fo.date), 0,
                                 ifelse(ctt$diabetes.ukn.fo.date<ctt$recruit.date, 1, 0))
ctt$diabetes.ukn.fo.pre.date <- ctt$diabetes.ukn.fo.date
ctt$diabetes.ukn.fo.pre.date[ctt$diabetes.ukn.fo.pre==0] <- NA

# apply the algorithm for pre-cruit dm 2021-01-07
# add unknown type diabetes have insulin within 12 month post-diagnosis to T1
ctt$diabetes.T1.fo.pre[ctt$diabetes.ukn.fo.pre==1 & !is.na(ctt$insulin_date) & 
                         ctt$insulin_date - ctt$diabetes.ukn.fo.pre.date <= 365 & 
                         ctt$insulin_date - ctt$diabetes.ukn.fo.pre.date >= 0] <- 1
ctt$diabetes.T1.fo.pre.date[is.na(ctt$diabetes.T1.fo.pre.date) & 
                              ctt$diabetes.T1.fo.pre==1] <- 
  ctt$diabetes.ukn.fo.pre.date[is.na(ctt$diabetes.T1.fo.pre.date) & 
                                 ctt$diabetes.T1.fo.pre==1]

# others T2
ctt$diabetes.T2.fo.pre[ctt$diabetes.ukn.fo.pre==1 & ctt$diabetes.T1.fo.pre==0] <- 1
ctt$diabetes.T2.fo.pre.date[is.na(ctt$diabetes.T2.fo.pre.date) & 
                              ctt$diabetes.T2.fo.pre==1] <- 
  ctt$diabetes.ukn.fo.pre.date[is.na(ctt$diabetes.T2.fo.pre.date) & 
                                 ctt$diabetes.T2.fo.pre==1]

# add people without any diabetes record (survey  + GP + HES + death)
# but have insulin or nonmetformin prescription after recruitment
# we have already add medication-identified dm for diabetes.fo.pre, so
# T1
# we don't know the diagnosis date for those with medication information only,
# so use insulin + prescription date<=20 yr

# the as.Date function assume the date and month of today for those with only year 
# so we minus the difference between the current date XXXX-XX-XX and the XXXX-01-01
# then add 182 to the mid of the year
dif <- as.numeric(Sys.Date() - as.Date(paste0(substr(as.character(Sys.Date()), 1, 4), "-01-01")))
ctt$birth.date.assumed <- as.Date(as.character(ctt$birth.year), format = "%Y") - dif +182

ctt$diabetes.T1.fo.pre[ctt$diabetes.T1.fo.pre==0 & !is.na(ctt$insulin_date) & 
                         ctt$insulin_date<ctt$recruit.date &  
                         ctt$insulin_date-ctt$birth.date.assumed<=365.25*20] <- 1 

ctt$diabetes.T1.fo.pre.date[is.na(ctt$diabetes.T1.fo.pre.date) & ctt$diabetes.T1.fo.pre==1] <- 
  ctt$medidate_pre[is.na(ctt$diabetes.T1.fo.pre.date) & ctt$diabetes.T1.fo.pre==1]
# above, not necessarily insulin date, but whatever date earilier

# T2
ctt$diabetes.T2.fo.pre[ctt$diabetes.T2.fo.pre==0 & ctt$diabetes.T1.fo.pre==0 & 
                         ctt$diabetes.fo.pre==1] <- 1

ctt$diabetes.T2.fo.pre.date[is.na(ctt$diabetes.T2.fo.pre.date) & ctt$diabetes.T2.fo.pre==1] <- 
  ctt$medidate_pre[is.na(ctt$diabetes.T2.fo.pre.date) & ctt$diabetes.T2.fo.pre==1]


# post
ctt$diabetes.T1.fo.post <- ifelse(ctt$diabetes.fo.pre==1 | is.na(ctt$diabetes.T1.fo.date), 0,  
                                  ifelse(ctt$diabetes.T1.fo.date>=ctt$recruit.date, 1, 0))
ctt$diabetes.T1.fo.post.date <- ctt$diabetes.T1.fo.date
ctt$diabetes.T1.fo.post.date[ctt$diabetes.T1.fo.post==0] <- NA

ctt$diabetes.T2.fo.post <- ifelse(ctt$diabetes.fo.pre==1 | is.na(ctt$diabetes.T2.fo.date), 0,
                                  ifelse(ctt$diabetes.T2.fo.date>=ctt$recruit.date, 1, 0))
ctt$diabetes.T2.fo.post.date <- ctt$diabetes.T2.fo.date
ctt$diabetes.T2.fo.post.date[ctt$diabetes.T2.fo.post==0] <- NA

ctt$diabetes.ukn.fo.post <- ifelse(ctt$diabetes.fo.pre==1 | is.na(ctt$diabetes.ukn.fo.date), 0,
                                   ifelse(ctt$diabetes.ukn.fo.date>=ctt$recruit.date, 1, 0))
ctt$diabetes.ukn.fo.post.date <- ctt$diabetes.ukn.fo.date
ctt$diabetes.ukn.fo.post.date[ctt$diabetes.ukn.fo.post==0] <- NA

# apply the algorithm 2020-12-01
# add unknown type diabetes have insulin within 12 month post-diagnosis to T1
ctt$diabetes.T1.fo.post[ctt$diabetes.ukn.fo.post==1 & !is.na(ctt$insulin_date) & 
      ctt$insulin_date - ctt$diabetes.ukn.fo.post.date <= 365 & 
      ctt$insulin_date - ctt$diabetes.ukn.fo.post.date >= 0] <- 1
ctt$diabetes.T1.fo.post.date[is.na(ctt$diabetes.T1.fo.post.date) & 
                               ctt$diabetes.T1.fo.post==1] <- 
  ctt$diabetes.ukn.fo.post.date[is.na(ctt$diabetes.T1.fo.post.date) & 
                                  ctt$diabetes.T1.fo.post==1]

# others T2
ctt$diabetes.T2.fo.post[ctt$diabetes.ukn.fo.post==1 & ctt$diabetes.T1.fo.post==0] <- 1
ctt$diabetes.T2.fo.post.date[is.na(ctt$diabetes.T2.fo.post.date) & 
                               ctt$diabetes.T2.fo.post==1] <- 
  ctt$diabetes.ukn.fo.post.date[is.na(ctt$diabetes.T2.fo.post.date) & 
                                  ctt$diabetes.T2.fo.post==1]

# add people without any diabetes record (survey  + GP + HES + death)
# but have insulin or nonmetformin prescription after recruitment
# T1
ctt$diabetes.T1.fo.post[ctt$diabetes.fo==0 & !is.na(ctt$insulin_date) &   
      ctt$insulin_date - ctt$recruit.date <= 365 & 
      ctt$insulin_date - ctt$recruit.date >= 0] <- 1

ctt$diabetes.T1.fo.post.date[is.na(ctt$diabetes.T1.fo.post.date) & ctt$diabetes.T1.fo.post==1] <- 
  ctt$medidate[is.na(ctt$diabetes.T1.fo.post.date) & ctt$diabetes.T1.fo.post==1]

# others T2
ctt$diabetes.T2.fo.post[ctt$diabetes.fo==0 & 
          ((!is.na(ctt$insulin_date) & ctt$insulin_date>ctt$recruit.date) | 
    (!is.na(ctt$nonmetformin_date) & ctt$nonmetformin_date>ctt$recruit.date)) &
      ctt$diabetes.T1.fo.post==0] <- 1

ctt$diabetes.T2.fo.post.date[is.na(ctt$diabetes.T2.fo.post.date) & ctt$diabetes.T2.fo.post==1] <- 
  ctt$medidate[is.na(ctt$diabetes.T2.fo.post.date) & ctt$diabetes.T2.fo.post==1]

ctt$medidate <- NULL
ctt$medidate1 <- NULL
ctt$medidate2 <- NULL
ctt$medidate_pre <- NULL

# update T1 and T2 post, because medication information may add new pre dm as post dm exist
ctt$diabetes.T1.fo.post[ctt$diabetes.T1.fo.pre==1] <- 0
ctt$diabetes.T1.fo.post.date[ctt$diabetes.T1.fo.post==0] <- NA

ctt$diabetes.T2.fo.post[ctt$diabetes.T2.fo.pre==1] <- 0
ctt$diabetes.T2.fo.post.date[ctt$diabetes.T2.fo.post==0] <- NA

# update the combination of post and pre
ctt$diabetes.T1.fo <- ctt$diabetes.T1.fo.pre
ctt$diabetes.T1.fo[ctt$diabetes.T1.fo.post==1] <- 1

ctt$diabetes.T1.fo.date <- ctt$diabetes.T1.fo.pre.date
ctt$diabetes.T1.fo.date[is.na(ctt$diabetes.T1.fo.date) & ctt$diabetes.T1.fo==1] <- 
  ctt$diabetes.T1.fo.post.date[is.na(ctt$diabetes.T1.fo.date) & ctt$diabetes.T1.fo==1]

ctt$diabetes.T2.fo <- ctt$diabetes.T2.fo.pre
ctt$diabetes.T2.fo[ctt$diabetes.T2.fo.post==1] <- 1

ctt$diabetes.T2.fo.date <- ctt$diabetes.T2.fo.pre.date
ctt$diabetes.T2.fo.date[is.na(ctt$diabetes.T2.fo.date) & ctt$diabetes.T2.fo==1] <- 
  ctt$diabetes.T2.fo.post.date[is.na(ctt$diabetes.T2.fo.date) & ctt$diabetes.T2.fo==1]



##### check death date earlier than events----

# found one case with stroke.all.date=2014-07-29; death.vascular.date=2014-07-03
# change stroke.all.date=2014-07-03
# rowname=27773
ctt[ctt$death.vascular.date<ctt$stroke.all.date & !is.na(ctt$death.vascular.date) & 
      !is.na(ctt$stroke.all.date), "stroke.all.date"] <- 
  ctt[ctt$death.vascular.date<ctt$stroke.all.date & !is.na(ctt$death.vascular.date) & 
      !is.na(ctt$stroke.all.date), "death.vascular.date"]

# diabetes.fo.post.date=2017-01-08
# death.date=2016-12.31
# rowname=371444
ctt[ctt$death.date<ctt$diabetes.fo.post.date & !is.na(ctt$death.date) & 
      !is.na(ctt$diabetes.fo.post.date), c("diabetes.fo.post.date")] <- 
  ctt[ctt$death.date<ctt$diabetes.fo.post.date & !is.na(ctt$death.date) & 
        !is.na(ctt$diabetes.fo.post.date), c("death.date")]

# death.date = 2010-08-28
# cancer.all.date = 2010-08-31
# case name = 81073
# recode cancer.all.date = death.date
ctt[ctt$death.date<ctt$cancer.incident.only.date & !is.na(ctt$death.date) & 
      !is.na(ctt$cancer.incident.only.date), "cancer.incident.only.date"] <- 
  ctt[ctt$death.date<ctt$cancer.incident.only.date & !is.na(ctt$death.date) & 
      !is.na(ctt$cancer.incident.only.date), "death.date"]



#############################################
##### other baseline disease algorithms #####
#############################################


# hypertension ------------------------------------------------------------

# verbal interview data
# code 1065: hypertension
# code 1072: essential hypertension
text1 <- "hpt_sr"
ctt[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  ctt[[text1]][bd[[text2]]== 1065 |
                 bd[[text2]]== 1072] <- 1
}

# ICD10: 
# first occurrence data
# I10 essential (primary) hypertension: f.131286.0.0
# I11 hypertensive heart disease: f.131288.0.0
# I12 hypertensive renal disease: f.131290.0.0
# I13 hypertensive heart and renal disease: f.131292.0.0
# I15 secondary hypertension: f.131294.0.0
ctt$hpt_fo <- 0
for (i in c("f.131286.0.0","f.131288.0.0","f.131290.0.0","f.131292.0.0", "f.131294.0.0")) {
  ctt$hpt_fo[bd[[i]]< ctt$recruit.date] <- 1
}

ctt$hypertension <- ifelse(ctt$hpt_fo==1 | ctt$hpt_sr==1, 1, 0)

ctt$hpt_fo <- NULL
ctt$hpt_sr <- NULL

# bp medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment  
bpmed <- read.csv(file.path(work_data, "HBPmed.csv"), header = F, encoding = "UTF-8")

# the first element is coded incorrectly
# revise it manually 
bpmed$V1[1] <- 1140860332

bpmed <- as.numeric(bpmed$V1)

ctt$BPMed_vi <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt$BPMed_vi[bd[[text2]] %in% bpmed] <- 1
}

# hypertension treatment
text1 <- "HBPtx"
text2 <- "hypertension"
text3 <- "BPMed_vi" 
ctt[[text1]] <- NA
ctt[[text1]][ctt[[text2]]==1 & ctt[[text3]]==1] <- 1
ctt[[text1]][ctt[[text2]]==0 | ctt[[text3]]==0] <- 0  




# statin medication -------------------------------------------------------

# 1141146234 atorvastatin
# 1141146138 lipitor 10mg tablet
# 1141192736 ezetimibe
# 1141192740 ezetrol 10mg tablet
# 1140888594 fluvastatin
# 1140864592 lescol 20mg capsule
# 1140888648 pravastatin
# 1140861970 lipostat 10mg tablet
# 1141192410 rosuvastatin
# 1141192414 crestor 10mg tablet
# 1140861958 simvastatin
# 1140881748 zocor 10mg tablet
# 1141200040 zocor heart-pro 10mg tablet
# 1141188146 simvador 10mg tablet
# 1140910632 eptastatin
# 1140910654 velastatin

sti <- c(1141146234, 1141146138, 1141192736, 1141192740, 1140888594, 1140864592, 
         1140888648, 1140861970, 1141192410, 1141192414, 1140861958, 1140881748, 
         1141200040, 1141188146, 1140910632, 1140910654)

# verbal interview
text1 <- "statin_vi"
ctt[[text1]] <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt[[text1]][bd[[text2]] %in% sti] <- 1
}

# # touch screen # we may not use touch screen data
# text1 <- "statin_ts"
# # gather all medication from multiple arrays
# text2 <- "f.6153.0.0"
# text3 <- "f.6153.0.1"
# text4 <- "f.6153.0.2"
# text5 <- "f.6153.0.3"
# text6 <- "f.6177.0.0"
# text7 <- "f.6177.0.1"
# text8 <- "f.6177.0.2"
# 
# ctt[[text1]] <- 0
# ctt[[text1]][bd[[text2]]== "Cholesterol lowering medication" |
#                bd[[text3]]== "Cholesterol lowering medication" |
#                bd[[text4]]== "Cholesterol lowering medication" |
#                bd[[text5]]== "Cholesterol lowering medication" |
#                bd[[text6]]== "Cholesterol lowering medication" |
#                bd[[text7]]== "Cholesterol lowering medication" |
#                bd[[text8]]== "Cholesterol lowering medication" ] <- 1
# # if the first array is missing, the others also are missing
# ctt[[text1]][is.na(bd[[text2]]) & is.na(bd[[text6]])] <- NA
# # only the first array has records of "Prefer not to answer", "Do not know"
# ctt[[text1]][bd[[text2]] %in% c("Prefer not to answer", "Do not know") | bd[[text6]] %in% c("Prefer not to answer", "Do not know")] <- NA

# combine statin
ctt$statin <- ctt$statin_vi
# ctt$statin[ctt$statin_ts==1] <- 1

rm(sti)
ctt$statin_ts <- NULL
ctt$statin_vi <- NULL

# by type
atorva <- c(1141146234, 1141146138)
fluva <- c(1140888594, 1140864592)
prava <- c(1140888648, 1140861970, 1140910632)
rosuva <- c(1141192410, 1141192414)
simva <- c(1140861958, 1140881748, 1141200040, 1141188146, 1140910654)
ezeti <- c(1141192736, 1141192740)

for (sti in c("atorva", "fluva", "prava", "rosuva", "simva", "ezeti")) {
  ctt[[sti]] <- 0
  for (j in 0:47) { # gather all records from multiple arrays
    text2 <- paste0("f.20003.0.", j)
    ctt[[sti]][bd[[text2]] %in% get(sti)] <- 1
  }
}

ctt$statin_sum <- rowSums(ctt[, c("atorva", "fluva", "prava", "rosuva", "simva", "ezeti")])
# # max statin_sum=3
# sum(ctt$statin_sum==3) # only 3 cases use 3 statins
# ctt[ctt$statin_sum==3, c("atorva", "fluva", "prava", "rosuva", "simva", "ezeti")]

# statin type combination
ctt$statin_type <- "none"
ctt$statin_type[ctt$statin_sum==1 & ctt$atorva==1] <- "atorva"
ctt$statin_type[ctt$statin_sum==1 & ctt$fluva==1] <- "fluva"
ctt$statin_type[ctt$statin_sum==1 & ctt$prava==1] <- "prava"
ctt$statin_type[ctt$statin_sum==1 & ctt$rosuva==1] <- "rosuva"
ctt$statin_type[ctt$statin_sum==1 & ctt$simva==1] <- "simva"
ctt$statin_type[ctt$statin_sum==1 & ctt$ezeti==1] <- "ezeti"

ctt$statin_type[ctt$statin_sum==2 & ctt$atorva==1 & ctt$ezeti==1] <- "atorva+ezeti"
ctt$statin_type[ctt$statin_sum==2 & ctt$fluva==1 & ctt$ezeti==1] <- "fluva+ezeti"
ctt$statin_type[ctt$statin_sum==2 & ctt$prava==1 & ctt$ezeti==1] <- "prava+ezeti"
ctt$statin_type[ctt$statin_sum==2 & ctt$rosuva==1 & ctt$ezeti==1] <- "rosuva+ezeti"
ctt$statin_type[ctt$statin_sum==2 & ctt$simva==1 & ctt$ezeti==1] <- "simva+ezeti"

ctt$statin_type[ctt$statin_sum==2 & ctt$atorva==1 & ctt$rosuva==1] <- "atorva+rosuva"
ctt$statin_type[ctt$statin_sum==2 & ctt$atorva==1 & ctt$simva==1] <- "atorva+simva"

ctt$statin_type[ctt$statin_sum==2 & ctt$simva==1 & ctt$rosuva==1] <- "simva+rosuva"
ctt$statin_type[ctt$statin_sum==2 & ctt$simva==1 & ctt$fluva==1] <- "simva+fluva"
ctt$statin_type[ctt$statin_sum==2 & ctt$simva==1 & ctt$prava==1] <- "simva+prava"

ctt$statin_type[ctt$statin_sum==3 & ctt$atorva==1 & ctt$rosuva==1 & ctt$ezeti==1] <- "atorva+rosuva+ezeti"
ctt$statin_type[ctt$statin_sum==3 & ctt$atorva==1 & ctt$simva==1 & ctt$ezeti==1] <- "atorva+simva+ezeti"

######
# collect dosage information from a few branded statin use with dosage
ctt$atorva10 <- 0
ctt$prava10 <- 0
ctt$rosuva10 <- 0
ctt$simva10 <- 0
ctt$fluva20 <- 0

a10 <- c(1141146138)
p10 <- c(1140861970)
r10 <- c(1141192414)
s10 <- c(1140881748, 1141200040, 1141188146)
f20 <- c(1140864592)
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt[["atorva10"]][bd[[text2]] %in% a10] <- 1
  ctt[["prava10"]][bd[[text2]] %in% p10] <- 1
  ctt[["rosuva10"]][bd[[text2]] %in% r10] <- 1
  ctt[["simva10"]][bd[[text2]] %in% s10] <- 1
  ctt[["fluva20"]][bd[[text2]] %in% f20] <- 1
}

ctt$sum_brd_statin <- ctt$atorva10+ctt$prava10+ctt$rosuva10+ctt$simva10+ctt$fluva20 

# check overlap, only one prava 10 + simva 10
# only keep simva 10
# so generate the combined variable for dosage
ctt$statin_dose_vi <- NA
ctt$statin_dose_vi[ctt$prava10==1] <- "prava10"
ctt$statin_dose_vi[ctt$fluva20==1] <- "fluva20"
ctt$statin_dose_vi[ctt$simva10==1] <- "simva10"
ctt$statin_dose_vi[ctt$atorva10==1] <- "atorva10"
ctt$statin_dose_vi[ctt$rosuva10==1] <- "rosuva10"

# merge to ctt2
ctt2 <- merge(ctt2, ctt[, c("id", "statin_dose_vi")])
saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"), compress = F)

ctt$atorva10 <- NULL
ctt$prava10 <- NULL
ctt$rosuva10 <- NULL
ctt$simva10 <- NULL
ctt$fluva20 <- NULL
ctt$sum_brd_statin <- NULL

######


# other CHD ---------------------------------------------------------------

# angina: verbal interview data, self-report
# code 1074: angina
# code 1076: heart failure/pulmonary odema

text1 <- "othCHD_sr"
ctt[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  ctt[[text1]][bd[[text2]]== 1074 | bd[[text2]]==1076] <- 1
}

# ICD 10
# I01: 131272.0.0
# I02: 131274.0.0
# I05-09: 131276.0.0, 131278.0.0, 131280.0.0, 131282.0.0, 131284.0.0
# I11: 131288.0.0
# I13: 131292.0.0
# I20: 131296.0.0
# I24: 131304.0.0
# I25: 131306.0.0 # not include, because old MI in this category
# I26-28: 131308.0.0, 131310.0.0, 131312.0.0
# I30-52: 131314.0.0, 131316.0.0, 131318.0.0, 131320.0.0, 131322.0.0, 131324.0.0, 
# 131326.0.0, 131328.0.0, 131330.0.0, 131332.0.0, 131334.0.0, 131336.0.0, 
# 131338.0.0, 131340.0.0, 131342.0.0, 131344.0.0, 131346.0.0, 131348.0.0, 
# 131350.0.0, 131352.0.0, 131354.0.0

oCHD <- c("131272.0.0", "131274.0.0", "131276.0.0", "131278.0.0", "131280.0.0", 
          "131282.0.0", "131284.0.0", "131288.0.0", "131292.0.0", "131296.0.0", 
          "131304.0.0", "131308.0.0", "131310.0.0", "131312.0.0", "131314.0.0", 
          "131316.0.0", "131318.0.0", "131320.0.0", "131322.0.0", "131324.0.0", 
          "131326.0.0", "131328.0.0", "131330.0.0", "131332.0.0", "131334.0.0", 
          "131336.0.0", "131338.0.0", "131340.0.0", "131342.0.0", "131344.0.0", 
          "131346.0.0", "131348.0.0", "131350.0.0", "131352.0.0", "131354.0.0")

ctt$othCHD_fo <- 0
for (i in oCHD){
  ctt$othCHD_fo[bd[[paste0("f.", i)]] < ctt$recruit.date] <- 1
}

# by eyeball check, those with I241 all have other CHD code, 
# so we can just keep them. no more action needed
# codes are not presented here

# because fo not include I25 due to I252, here use inpatient data to include
# inpatient record for chronic IHD 
# ICD10 I25 but I252
# ICD9 414 
# use temp to collect all dates with code that match the rule
temp <- as.data.frame(ctt$id)
for (i in 0:212) { # f.41270 has 213 arrays
  text1 <- paste0("CHD.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("^I25", bd[[text2]]) & 
                  bd[[text2]]!="I252" & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("^I25", bd[[text2]]) & 
                  bd[[text2]]!="I252" & 
                  !is.na(bd[[text2]])]
}

# all 414* 
for (i in 0:46) { # f.41271 has 47 arrays
  text1 <- paste0("CHD.date.", i+213)
  text2 <- paste0("f.41271.0.",i)
  text3 <- paste0("f.41281.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("414", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("414", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

# get rid of date after recruitment date
for (i in 0:259){
  text1 <- paste0("CHD.date.", i)
  temp[[text1]][temp[[text1]] >= ctt$recruit.date] <- NA
}
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# at least one column have a date (before recruitment), recode 1, otherwise 0 
temp$`ctt$id` <- NULL
temp$CHD.i25 <- ifelse(rowSums(is.na(temp)) < ncol(temp), 1, 0)

ctt$othCHD <- ifelse(ctt$othCHD_sr==1 | ctt$othCHD_fo==1 | temp$CHD.i25==1, 1, 0)
rm(oCHD, temp)



# peripheral arterial disease ---------------------------------------------
# PAD: verbal interview data
# code 1067: PVD
# code 1087: leg claudication/intermittent claudication
# code 1088: arterial embolism
# code 1492: aortic aneurysm
# code 1591: aortic aneurysm rupture
# code 1592: aortic dissection
# further add codes from operations

text1 <- "PVD.sr"
ctt[[text1]] <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  ctt[[text1]][bd[[text2]]== 1067 |
                 bd[[text2]]== 1087 |
                 bd[[text2]]== 1088 |
                 bd[[text2]]== 1492 |
                 bd[[text2]]== 1591 |
                 bd[[text2]]== 1592] <- 1
}
for (k in 0:31) { # gather all records from multiple arrays
  text3 <- paste0("f.20004.0.", k)
  ctt[[text1]][bd[[text3]]== 1071 |
                 bd[[text3]]== 1102 |
                 bd[[text3]]== 1103 |
                 bd[[text3]]== 1555 |
                 bd[[text3]]== 1104 |
                 bd[[text3]]== 1105 |
                 # bd[[text3]]== 1106 |
                 bd[[text3]]== 1107 |
                 bd[[text3]]== 1108 |
                 bd[[text3]]== 1109 |
                 bd[[text3]]== 1110 |
                 bd[[text3]]== 1440 |
                 bd[[text3]]== 1441 |
                 bd[[text3]]== 1442 |
                 bd[[text3]]== 1443] <- 1
}

# first occurrence
# I71, I72, I73, I74, I77
ctt$PVD_fo <- 0
for (i in c("f.131382.0.0","f.131384.0.0","f.131386.0.0", "f.131388.0.0","f.131390.0.0")) {
  ctt$PVD_fo_date <- bd[[i]]
  ctt$PVD_fo_date[ctt$PVD_fo_date >= ctt$recruit.date] <- NA
  ctt$PVD_fo[!is.na(ctt$PVD_fo_date)] <- 1
}

# PAD OPCS code
# OPCS 4
ctt$PVD.opcs <- 0
for (i in 0:116) { # f.41272 has 117 arrays
  text2 <- paste0("f.41272.0.",i)
  ctt$PVD.opcs[grepl("^L16", bd[[text2]]) |
                 grepl("^L18", bd[[text2]]) |
                 grepl("^L19", bd[[text2]]) |
                 grepl("^L2", bd[[text2]]) |
                 grepl("^L30", bd[[text2]]) |
                 grepl("^L31", bd[[text2]]) |
                 grepl("^L37", bd[[text2]]) |
                 grepl("^L38", bd[[text2]]) |
                 grepl("^L39", bd[[text2]]) |
                 grepl("^L4", bd[[text2]]) |
                 grepl("^L5", bd[[text2]]) |
                 grepl("^L60", bd[[text2]]) |
                 grepl("^L62", bd[[text2]]) |
                 grepl("^L63", bd[[text2]]) |
                 grepl("^L65", bd[[text2]]) |
                 grepl("^L66", bd[[text2]]) |
                 grepl("^L67", bd[[text2]]) |
                 grepl("^L68", bd[[text2]]) |
                 grepl("^L70", bd[[text2]]) |
                 grepl("^L71", bd[[text2]]) |
                 grepl("^L74", bd[[text2]]) |
                 grepl("^L75", bd[[text2]]) |
                 grepl("^L76", bd[[text2]]) |
                 grepl("^L89", bd[[text2]])] <- 1
}
# OPCS 3
for (i in 0:15) {
  text2 <- paste0("f.41273.0.",i)
  ctt$PVD.opcs[grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]])] <- 1
}

temp <- as.data.frame(ctt$PVD.opcs)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:116) {
  text1 <- paste0("PVD.opcs.date.", i)
  text2 <- paste0("f.41272.0.",i)
  text3 <- paste0("f.41282.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][(grepl("^L16", bd[[text2]]) |
                   grepl("^L18", bd[[text2]]) |
                   grepl("^L19", bd[[text2]]) |
                   grepl("^L2", bd[[text2]]) |
                   grepl("^L30", bd[[text2]]) |
                   grepl("^L31", bd[[text2]]) |
                   grepl("^L37", bd[[text2]]) |
                   grepl("^L38", bd[[text2]]) |
                   grepl("^L39", bd[[text2]]) |
                   grepl("^L4", bd[[text2]]) |
                   grepl("^L5", bd[[text2]]) |
                   grepl("^L60", bd[[text2]]) |
                   grepl("^L62", bd[[text2]]) |
                   grepl("^L63", bd[[text2]]) |
                   grepl("^L65", bd[[text2]]) |
                   grepl("^L66", bd[[text2]]) |
                   grepl("^L67", bd[[text2]]) |
                   grepl("^L68", bd[[text2]]) |
                   grepl("^L70", bd[[text2]]) |
                   grepl("^L71", bd[[text2]]) |
                   grepl("^L74", bd[[text2]]) |
                   grepl("^L75", bd[[text2]]) |
                   grepl("^L76", bd[[text2]]) |
                   grepl("^L89", bd[[text2]])) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][(grepl("^L16", bd[[text2]]) |
                   grepl("^L18", bd[[text2]]) |
                   grepl("^L19", bd[[text2]]) |
                   grepl("^L2", bd[[text2]]) |
                   grepl("^L30", bd[[text2]]) |
                   grepl("^L31", bd[[text2]]) |
                   grepl("^L37", bd[[text2]]) |
                   grepl("^L38", bd[[text2]]) |
                   grepl("^L39", bd[[text2]]) |
                   grepl("^L4", bd[[text2]]) |
                   grepl("^L5", bd[[text2]]) |
                   grepl("^L60", bd[[text2]]) |
                   grepl("^L62", bd[[text2]]) |
                   grepl("^L63", bd[[text2]]) |
                   grepl("^L65", bd[[text2]]) |
                   grepl("^L66", bd[[text2]]) |
                   grepl("^L67", bd[[text2]]) |
                   grepl("^L68", bd[[text2]]) |
                   grepl("^L70", bd[[text2]]) |
                   grepl("^L71", bd[[text2]]) |
                   grepl("^L74", bd[[text2]]) |
                   grepl("^L75", bd[[text2]]) |
                   grepl("^L76", bd[[text2]]) |
                   grepl("^L89", bd[[text2]])) & 
                  !is.na(bd[[text2]])]
}
# ICD-9
for (i in 0:15) {
  text1 <- paste0("PVD.opcs.date.", i+117)
  text2 <- paste0("f.41273.0.",i)
  text3 <- paste0("f.41283.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]]) & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][grepl("88", bd[[text2]]) & !grepl("^888", bd[[text2]]) & 
                  !is.na(bd[[text2]])]
}

# if there are multiple date use the earliest one
# get rid of date after recruitment date and first column
for (i in 0:132){
  text1 <- paste0("PVD.opcs.date.", i)
  temp[[text1]][temp[[text1]] >= ctt$recruit.date] <- NA
}
temp$`ctt$PVD.opcs` <- NULL
# code PVD.opcs=0 if all dates are NAs, which means no PVD before recruitment
ctt$PVD.opcs[rowSums(is.na(temp)) == ncol(temp)] <- 0

# ctt$PVD_2 <- ifelse(ctt$PVD_fo==1 | ctt$PVD.opcs==1, 1, 0)

ctt$PVD <- ifelse(ctt$PVD.sr==1 | ctt$PVD_fo==1 | ctt$PVD.opcs==1, 1, 0)

ctt$PVD_fo <- NULL
ctt$PVD.sr <- NULL
ctt$PVD_fo_date <- NULL
ctt$PVD.opcs <- NULL


# CVD history (primary/secondary) -----------------------------------------

ctt$CVhist <- ctt$MI.baseline + ctt$stroke.baseline + ctt$PVD + ctt$othCHD

ctt$CVD <- ifelse(ctt$CVhist==0, "None",
                  ifelse(ctt$CVhist==1 & ctt$MI.baseline==1, "MI only",  
                         ifelse(ctt$CVhist==1 & ctt$stroke.baseline==1, "Stroke only", 
                                ifelse(ctt$CVhist==1 & ctt$PVD==1, "PVD only", 
                                       ifelse(ctt$CVhist==1 & ctt$othCHD==1, "other CHD only", "Two or more")))))
ctt$CVD <- relevel(as.factor(ctt$CVD), ref = "None")



#########################################
##### algorithms of other variables #####
#########################################

# smoke status ####
ctt$smoke <- bd$f.20116.0.0
ctt$smoke[ctt$smoke=="Prefer not to answer"] <- NA
ctt$smoke <- factor(ctt$smoke)

# impute to the majority by gender, age groups and education levels
ctt <- ctt %>% 
  group_by(sex, age.group, education) %>% 
  add_count(smoke) %>% 
  mutate(smoke.imp = if_else(is.na(smoke), smoke[which.max(n)], smoke)) %>% 
  select(-n) %>% ungroup()

# for QRISK calculation

# Current cigarettes per day: impute to the majority level, i.e. 10-20

ctt$smoke.qrisk <- ifelse(ctt$smoke.imp=="Never", 1, 
                           ifelse(ctt$smoke.imp=="Previous", 2, 
                                  ifelse(ctt$smoke.imp=="Current", 4, NA)))

# after tabulation the three level, moderate level has the most people: 15438; heavy 13512 and light 7213

# count cigarette per day
ctt$smoke.qrisk[ctt$smoke.qrisk==4 & ((bd$f.3456.0.0<10 & bd$f.3456.0.0>0) | bd$f.3456.0.0==-10)] <- 3 # less than 10. bd$f.3456.0.0==-10 means less than 1 

ctt$smoke.qrisk[ctt$smoke.qrisk==4 & bd$f.3456.0.0>=20] <- 5 # 20+


# LDL & HDL ####
ctt$LDL <- bd$f.30780.0.0
ctt$HDL <- bd$f.30760.0.0

# blood creatinine ####
ctt$creatinine <- bd$f.30700.0.0
ctt$creatinine_ln <- log(ctt$creatinine)

# blood pressure ####
# 93,94 (manual, alternative), 4080,4079 (automatic, preferred), they rarely overlap
# two readings 
# systolic BP
text1 <- "manuSBP"
# two close readings, take average 
text2 <- "f.93.0.0"
text3 <- "f.93.0.1"
ctt[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                       (bd[[text2]]+bd[[text3]])/2, 
                       ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "autoSBP"
# two close readings, take average 
text2 <- "f.4080.0.0"
text3 <- "f.4080.0.1"
ctt[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                       (bd[[text2]]+bd[[text3]])/2, 
                       ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
# automatic reading is the priority
text1 <- "SBP"
text2 <- "autoSBP"
text3 <- "manuSBP"
ctt[[text1]] <- ifelse(!is.na(ctt[[text2]]), ctt[[text2]], ctt[[text3]])
# diastolic BP
text1 <- "manuDBP"
# two close readings, take average 
text2 <- "f.94.0.0"
text3 <- "f.94.0.1"
ctt[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                       (bd[[text2]]+bd[[text3]])/2, 
                       ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "autoDBP"
# two close readings, take average 
text2 <- "f.4079.0.0"
text3 <- "f.4079.0.1"
ctt[[text1]] <- ifelse(!is.na(bd[[text2]]) & !is.na(bd[[text3]]), 
                       (bd[[text2]]+bd[[text3]])/2, 
                       ifelse(is.na(bd[[text2]]), bd[[text3]], bd[[text2]]))
text1 <- "DBP"
text2 <- "autoDBP"
text3 <- "manuDBP"
ctt[[text1]] <- ifelse(!is.na(ctt[[text2]]), ctt[[text2]], ctt[[text3]])
# delete the intermedium variables
ctt$autoSBP <- NULL
ctt$autoDBP <- NULL
ctt$manuDBP <- NULL
ctt$manuSBP <- NULL

# ethnicity ####
white <- c("White", "British", "Irish", "Any other white background")
black <- c("Black or Black British", "Caribbean", "African", "Any other Black background")
s.asian <- c("Asian or Asian British", "Indian", "Pakistani", "Bangladeshi")
other <- c("Mixed", "Chinese", "Other ethnic group", "White and Black Caribbean",
           "White and Black African", "White and Asian", "Any other mixed background",
           "Any other Asian background")
norecord <- c("Prefer not to answer", "Do not know")

text1 <- "ethnicity"
text2 <- "f.21000.0.0"
ctt[[text1]] <- NA
ctt[[text1]][bd[[text2]]%in% white] <- "White"
ctt[[text1]][bd[[text2]]%in% black] <- "Black"
ctt[[text1]][bd[[text2]]%in% s.asian] <- "South Asian"
ctt[[text1]][bd[[text2]]%in% other] <- "Others"
ctt[[text1]][bd[[text2]]%in% norecord] <- NA

ctt$ethn2 <- ifelse(ctt$ethnicity=="Black", "Black", 
                    ifelse(ctt$ethnicity=="White", "White", "Others"))
ctt$ethn2 <- relevel(as.factor(ctt$ethn2), ref = "White")

# BMI ####
text1 <- "BMI"
text2 <- "f.21001.0.0"
ctt[[text1]] <- bd[[text2]]
# BMI category
text1 <- "BMI_cat"
text2 <- "BMI"
ctt[[text1]] <- NA
ctt[[text1]][ctt[[text2]]<18.5] <- "<18.5"
ctt[[text1]][ctt[[text2]]>=18.5 & ctt[[text2]]<25] <- "18.5-25"
ctt[[text1]][ctt[[text2]]>=25 & ctt[[text2]]<30] <- "25-30"
ctt[[text1]][ctt[[text2]]>=30 & ctt[[text2]]<35] <- "30-35"
ctt[[text1]][ctt[[text2]]>=35 & ctt[[text2]]<40] <- "35-40"
ctt[[text1]][ctt[[text2]]>=40] <- "40+"

# highest education ####
# categories as below in hierarchy
degree <- c("College or University degree")
alevel <- c("A levels/AS levels or equivalent", "NVQ or HND or HNC or equivalent"
            , "Other professional qualifications eg: nursing, teaching")
olevel <- c("O levels/GCSEs or equivalent", "CSEs or equivalent")
none <- c("None of the above")

text1 <- "education"
text2 <- "f.6138.0.0"
text3 <- "f.6138.0.1"
text4 <- "f.6138.0.2"
text5 <- "f.6138.0.3"
text6 <- "f.6138.0.4"
text7 <- "f.6138.0.5"
ctt[[text1]] <- NA
ctt[[text1]][bd[[text2]]%in% none | bd[[text3]]%in% none |
               bd[[text4]]%in% none | bd[[text5]]%in% none |
               bd[[text6]]%in% none | bd[[text7]]%in% none] <- "None"
ctt[[text1]][bd[[text2]]%in% olevel | bd[[text3]]%in% olevel |
               bd[[text4]]%in% olevel | bd[[text5]]%in% olevel |
               bd[[text6]]%in% olevel | bd[[text7]]%in% olevel] <- "O level or CSE"
ctt[[text1]][bd[[text2]]%in% alevel | bd[[text3]]%in% alevel |
               bd[[text4]]%in% alevel | bd[[text5]]%in% alevel |
               bd[[text6]]%in% alevel | bd[[text7]]%in% alevel] <- 
  "A level, vocational diploma or professional qualification"
ctt[[text1]][bd[[text2]]%in% degree | bd[[text3]]%in% degree |
               bd[[text4]]%in% degree | bd[[text5]]%in% degree |
               bd[[text6]]%in% degree | bd[[text7]]%in% degree] <- "College or University degree"


# marital status ####
# use relation of people in household to the participant as indicator
cohabit <- c("Husband, wife or partner")

text1 <- "marital"
text2 <- "f.6141.0.0"
text3 <- "f.6141.0.1"
text4 <- "f.6141.0.2"
text5 <- "f.6141.0.3"
text6 <- "f.6141.0.4"
ctt[[text1]] <- "Not living with spouse, partner or cohabitee"
ctt[[text1]][bd[[text2]]== "Prefer not to answer" | 
               bd[[text3]]== "Prefer not to answer" |
               bd[[text4]]== "Prefer not to answer" | 
               bd[[text5]]== "Prefer not to answer" |
               bd[[text6]]== "Prefer not to answer"] <- NA
ctt[[text1]][bd[[text2]]%in% cohabit | bd[[text3]]%in% cohabit |
               bd[[text4]]%in% cohabit | bd[[text5]]%in% cohabit |
               bd[[text6]]%in% cohabit] <- "married/partnership"


# IMD: England, Wales, and Scotland separately ####
# generate residential country using availability of IMD
ctt$England <- ifelse(is.na(bd$f.26410.0.0), 0, 1)
ctt$Wales <- ifelse(is.na(bd$f.26426.0.0), 0, 1)
ctt$Scotland <- ifelse(is.na(bd$f.26427.0.0), 0, 1)

# according to the imd_baseline.pdf from UKB, assign the IMD sources
# IMD2004: 2004, 2005, 2006
# IMD2007: 2007, 2008, 2009
# IMD2010: 2010
# WIMD2005: 2006, 2007
# WIMD2008: 2008, 2009, 2010
# SIMD2006: 2006, 2007, 2008
# SIMD2009: 2009, 2010

ctt$IMD.source <- NA
ctt$IMD.source[ctt$England==1 & ctt$recruit.date < "2007-01-01"] <- "IMD2004"
ctt$IMD.source[ctt$England==1 & ctt$recruit.date >= "2007-01-01" 
               & ctt$recruit.date< "2010-01-01"] <- "IMD2007"
ctt$IMD.source[ctt$England==1 & ctt$recruit.date >= "2010-01-01"] <- "IMD2010"

ctt$IMD.source[ctt$Scotland==1 & ctt$recruit.date < "2009-01-01"] <- "SIMD2006"
ctt$IMD.source[ctt$Scotland==1 & ctt$recruit.date >= "2009-01-01"] <- "SIMD2009"

ctt$IMD.source[ctt$Wales==1 & ctt$recruit.date < "2008-01-01"] <- "WIMD2005"
ctt$IMD.source[ctt$Wales==1 & ctt$recruit.date >= "2008-01-01"] <- "WIMD2008"

# Quintile cutoffs of IMD of different sources
# IMD2004: 8.35; 13.72; 21.15; 34.20
# IMD2007: 8.32; 13.74; 21.22; 34.42
# IMD2010: 8.49; 13.79; 21.35; 34.17
# SIMD2006: 7.75; 13.56; 21.05; 33.70
# SIMD2009: 7.76; 13.76; 21.02; 33.72
# WIMD2005: 9.96; 14.94; 21.16; 32.70
# WIMD2008: 9.8; 14.8; 21.2; 32.5

# generate IMD quintile
ctt$IMD.Q5 <- NA
ctt$IMD.Q5[ctt$IMD.source=="IMD2004" & bd$f.26410.0.0 <= 8.35] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="IMD2004" & bd$f.26410.0.0 > 8.35 & bd$f.26410.0.0 <= 13.72] <- 2
ctt$IMD.Q5[ctt$IMD.source=="IMD2004" & bd$f.26410.0.0 > 13.72 & bd$f.26410.0.0 <= 21.15] <- 3
ctt$IMD.Q5[ctt$IMD.source=="IMD2004" & bd$f.26410.0.0 > 21.15 & bd$f.26410.0.0 <= 34.20] <- 4
ctt$IMD.Q5[ctt$IMD.source=="IMD2004" & bd$f.26410.0.0 > 34.20] <- 5

ctt$IMD.Q5[ctt$IMD.source=="IMD2007" & bd$f.26410.0.0 <= 8.32] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="IMD2007" & bd$f.26410.0.0 > 8.32 & bd$f.26410.0.0 <= 13.74] <- 2
ctt$IMD.Q5[ctt$IMD.source=="IMD2007" & bd$f.26410.0.0 > 13.74 & bd$f.26410.0.0 <= 21.22] <- 3
ctt$IMD.Q5[ctt$IMD.source=="IMD2007" & bd$f.26410.0.0 > 21.22 & bd$f.26410.0.0 <= 34.42] <- 4
ctt$IMD.Q5[ctt$IMD.source=="IMD2007" & bd$f.26410.0.0 > 34.42] <- 5

ctt$IMD.Q5[ctt$IMD.source=="IMD2010" & bd$f.26410.0.0 <= 8.49] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="IMD2010" & bd$f.26410.0.0 > 8.49 & bd$f.26410.0.0 <= 13.79] <- 2
ctt$IMD.Q5[ctt$IMD.source=="IMD2010" & bd$f.26410.0.0 > 13.79 & bd$f.26410.0.0 <= 21.35] <- 3
ctt$IMD.Q5[ctt$IMD.source=="IMD2010" & bd$f.26410.0.0 > 21.35 & bd$f.26410.0.0 <= 34.17] <- 4
ctt$IMD.Q5[ctt$IMD.source=="IMD2010" & bd$f.26410.0.0 > 34.17] <- 5

ctt$IMD.Q5[ctt$IMD.source=="SIMD2006" & bd$f.26427.0.0 <= 7.75] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 7.75 & bd$f.26427.0.0 <= 13.56] <- 2
ctt$IMD.Q5[ctt$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 13.56 & bd$f.26427.0.0 <= 21.05] <- 3
ctt$IMD.Q5[ctt$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 21.05 & bd$f.26427.0.0 <= 33.70] <- 4
ctt$IMD.Q5[ctt$IMD.source=="SIMD2006" & bd$f.26427.0.0 > 33.70] <- 5

ctt$IMD.Q5[ctt$IMD.source=="SIMD2009" & bd$f.26427.0.0 <= 7.76] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 7.76 & bd$f.26427.0.0 <= 13.76] <- 2
ctt$IMD.Q5[ctt$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 13.76 & bd$f.26427.0.0 <= 21.02] <- 3
ctt$IMD.Q5[ctt$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 21.02 & bd$f.26427.0.0 <= 33.72] <- 4
ctt$IMD.Q5[ctt$IMD.source=="SIMD2009" & bd$f.26427.0.0 > 33.72] <- 5

ctt$IMD.Q5[ctt$IMD.source=="WIMD2005" & bd$f.26426.0.0 <= 9.96] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 9.96 & bd$f.26426.0.0 <= 14.94] <- 2
ctt$IMD.Q5[ctt$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 14.94 & bd$f.26426.0.0 <= 21.16] <- 3
ctt$IMD.Q5[ctt$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 21.16 & bd$f.26426.0.0 <= 32.70] <- 4
ctt$IMD.Q5[ctt$IMD.source=="WIMD2005" & bd$f.26426.0.0 > 32.70] <- 5

ctt$IMD.Q5[ctt$IMD.source=="WIMD2008" & bd$f.26426.0.0 <= 9.8] <- 1 # wealthiest
ctt$IMD.Q5[ctt$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 9.8 & bd$f.26426.0.0 <= 14.8] <- 2
ctt$IMD.Q5[ctt$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 14.8 & bd$f.26426.0.0 <= 21.2] <- 3
ctt$IMD.Q5[ctt$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 21.2 & bd$f.26426.0.0 <= 32.5] <- 4
ctt$IMD.Q5[ctt$IMD.source=="WIMD2008" & bd$f.26426.0.0 > 32.5] <- 5

# Townsend score #### 
# regress on IMD score, years and countries

# Manually use Ordnance Survey (OSGB) rounded co-ordinates to check post codes 
# linked to ward Townsend score 2001 for missing for both Townsend and IMD. 
# Just use the authority level average score for the three missing for all these variables

ctt$townsend <- bd$f.189.0.0

# regression first
# generate the combination of IMD scores from the three countries
ctt$IMD <- ifelse(!is.na(bd$f.26410.0.0), bd$f.26410.0.0, 
                   ifelse(!is.na(bd$f.26427.0.0), bd$f.26427.0.0, 
                          ifelse(!is.na(bd$f.26426.0.0),bd$f.26426.0.0, NA)))

fit <- lm(townsend~IMD+IMD.source, ctt)
ctt$townsend.pred <- NA
ctt$townsend.pred[!is.na(ctt$IMD)] <- predict(fit, newdata=ctt[!is.na(ctt$IMD),])

ctt$townsend.imp <- ctt$townsend
ctt$townsend.imp[is.na(ctt$townsend.imp)] <- ctt$townsend.pred[is.na(ctt$townsend.imp)]

# manually lookup townsend score
osgb <- osgb %>% rename(id = f.eid) %>% select(id, townsend.2001)
ctt <- left_join(ctt, osgb)
ctt$townsend.imp[is.na(ctt$townsend.imp)] <- ctt$townsend.2001[is.na(ctt$townsend.imp)]

# finally, check the assessment centres of the remaining three missing values
# f.54.0.0 = 11008 Bury, Townsend 0.0594
# f.54.0.0 = 11011 Bristol, Townsend 1.1034
# another one is missing for all values, will be deleted later
ctt$townsend.imp[is.na(ctt$townsend.imp) & bd$f.54.0.0==11008] <- 0.0594
ctt$townsend.imp[is.na(ctt$townsend.imp) & bd$f.54.0.0==11011] <- 1.1034



# add variable for multiple imputation use --------------------------------

ctt$weight <- bd$f.21002.0.0
ctt$height <- sqrt(ctt$weight/ctt$BMI)*100
ctt$cholesterol <- bd$f.30690.0.0
ctt$triglycerides <- bd$f.30870.0.0

# systolic BP read 1 & 2
ctt$SBP1 <- bd$f.4080.0.0
ctt$SBP1[is.na(ctt$SBP1)] <- bd$f.93.0.0[is.na(ctt$SBP1)]

ctt$SBP2 <- bd$f.4080.0.1
ctt$SBP2[is.na(ctt$SBP2)] <- bd$f.93.0.1[is.na(ctt$SBP2)]

# Diastolic BP read 1 & 2
ctt$DBP1 <- bd$f.4079.0.0
ctt$DBP1[is.na(ctt$DBP1)] <- bd$f.94.0.0[is.na(ctt$DBP1)]

ctt$DBP2 <- bd$f.4079.0.1
ctt$DBP2[is.na(ctt$DBP2)] <- bd$f.94.0.1[is.na(ctt$DBP2)]

# duration of baseline diabetes
ctt$diabetes.duration = as.numeric(ctt$recruit.date - ctt$diabetes.fo.pre.date)

ctt <- ctt %>% 
  mutate(diabetes.fo.pre.duration = 
           case_when(is.na(diabetes.duration) ~ "no diabetes",
                     diabetes.duration %/% (365.25) == 0 ~ "in 1 year",
                     diabetes.duration %/% (365.25) == 1 ~ "1 year ago",
                     diabetes.duration %/% (365.25) == 2 ~ "2 years ago",
                     diabetes.duration %/% (365.25) == 3 ~ "3 years ago",
                     diabetes.duration %/% (365.25) == 4 ~ "4 years ago",
                     diabetes.duration %/% (365.25) == 5 ~ "5 years ago",
                     diabetes.duration %/% (365.25) == 6 ~ "6 years ago",
                     diabetes.duration %/% (365.25) == 7 ~ "7 years ago",
                     diabetes.duration %/% (365.25) == 8 ~ "8 years ago",
                     diabetes.duration %/% (365.25) == 9 ~ "9 years ago",
                     diabetes.duration %/% (365.25) %in% 10:14 ~ "10-14 years ago",
                     diabetes.duration %/% (365.25) %in% 15:19 ~ "15-19 years ago",
                     diabetes.duration %/% (365.25) >=20 ~ "20 plus years ago"
           ))

ctt$diabetes.fo.pre.duration <- factor(ctt$diabetes.fo.pre.duration, 
             levels = c("no diabetes", "in 1 year", "1 year ago", "2 years ago", "3 years ago", 
                  "4 years ago", "5 years ago", "6 years ago",
                  "7 years ago", "8 years ago", "9 years ago",
                  "10-14 years ago", "15-19 years ago", "20 plus years ago"))


ctt <- ctt %>% 
  mutate(diabetes.fo.pre.duration5 = 
           case_when(is.na(diabetes.duration) ~ "no diabetes",
                     diabetes.duration %/% (365.25) <= 4 ~ "<5 years",
                     diabetes.duration %/% (365.25) %in% 5:9 ~ "5-10 years",
                     diabetes.duration %/% (365.25) %in% 10:14 ~ "10-14 years",
                     diabetes.duration %/% (365.25) %in% 15:19 ~ "15-19 years",
                     diabetes.duration %/% (365.25) >=20 ~ "20 plus years"

))

ctt$diabetes.fo.pre.duration5 <- factor(ctt$diabetes.fo.pre.duration5, 
                       levels = c("no diabetes", "<5 years", "5-10 years", 
                          "10-14 years", "15-19 years", "20 plus years"))

# end


# add gp record tag into ctt####
gp_clin <- read.delim(file.path(raw_data, "gp_clinical.txt"))

# keep unique id
gp <- data.frame(id=unique(gp_clin$eid))
gp$gp_record <- 1

ctt <- left_join(ctt, gp, by = "id")
ctt$gp_record[is.na(ctt$gp_record)] <- 0

# add loss of follow-up date
ctt$lossfudate <- as.Date(bd$f.191.0.0)


# add familial hypercholesterolaemia---- 
temp <- as.data.frame(ctt$id)
# above is an easy way to create a date frame with the same structure with ctt
for (i in 0:212) { # f.41280 has 213 arrays, corresponding to 41270
  text1 <- paste0("fam.hypercho.date.", i)
  text2 <- paste0("f.41270.0.",i)
  text3 <- paste0("f.41280.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][bd[[text2]]== "E780" &
                  !is.na(bd[[text2]])] <-
    bd[[text3]][bd[[text2]]== "E780" &
                  !is.na(bd[[text2]])]
}

for (i in 0:46) { # f.41271 has 47 arrays
  text1 <- paste0("fam.hypercho.date.", i+213)
  text2 <- paste0("f.41271.0.",i)
  text3 <- paste0("f.41281.0.",i)
  temp[[text1]] <- NA
  temp[[text1]] <- as.Date(temp[[text1]])
  temp[[text1]][bd[[text2]]=="27200" & 
                  !is.na(bd[[text2]])] <- 
    bd[[text3]][bd[[text2]]=="27200" & 
                  !is.na(bd[[text2]])]
}

# get rid of date after recruitment date
for (i in 0:259){
  text1 <- paste0("fam.hypercho.date.", i)
  temp[[text1]][temp[[text1]] >= ctt$recruit.date] <- NA
}
# remove columns with all NAs
temp <- temp[, colSums(is.na(temp))<nrow(temp)]
# at least one column have a date (before recruitment), recode 1, otherwise 0 
temp$`ctt$id` <- NULL
temp$fam.hypercho <- ifelse(rowSums(is.na(temp)) < ncol(temp), 1, 0)

# combine to ctt
# maybe it is more correct to just call it hypercholaemia
ctt$hypercho <- temp$fam.hypercho

# CKD ----
# use ICD10 N18: chronic renal failure
ctt$CKD <- ifelse(bd$f.132032.0.0<ctt$recruit.date, 1, 0)
ctt$CKD[is.na(ctt$CKD)] <- 0

# Hba1c
ctt$hba1c <- bd$f.30750.0.0

# according to UKB showcase serum_hb1ac.pdf,
# the machine's analytical range is 15-184 mmol/mol
# 5 case >184 coded as 184
# plus the 5 cases defined as no baseline diabetes

ctt$hba1c[ctt$hba1c>184] <- 184

ctt$hba1c_ln <- log(ctt$hba1c)


# use the IPAQ group
# 22032
ctt$PA_group <- bd$f.22032.0.0

ctt$PA_group[is.na(ctt$PA_group)] <- "missing"
ctt$PA_group <- factor(ctt$PA_group, 
                       levels = c("moderate", "low", "high", "missing"))

# neuroticism score
ctt$neuroticism <- bd$f.20127.0.0



# variables for calculate qrisk -------------------------------------------


# atrial fibrillation
# use first occurence data
# ICD: I48
ctt$atrial_fibrillation <- ifelse(bd$f.131350.0.0<ctt$recruit.date, 1, 0)
ctt$atrial_fibrillation[is.na(ctt$atrial_fibrillation)] <- 0

# relative angina or heart attack
# use "Heart disease" in UKB instead of ...
# father, mother, siblings
# no age known, so assume <60
ctt$heart_attack_relative <- 0
for (j in 0:9) { # father illness
  text2 <- paste0("f.20107.0.", j)
  ctt$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}

for (j in 0:10) { # mother illness
  text2 <- paste0("f.20110.0.", j)
  ctt$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}

for (j in 0:11) { # sibling illness
  text2 <- paste0("f.20111.0.", j)
  ctt$heart_attack_relative[bd[[text2]]== "Heart disease"] <- 1
}
# only with additional information about father/mother's age and age at death, 
# we can only assure some heart disease happened <60, but we cannot rule out 
# which heart disease that definitely happened >=60 


# ED
# only in nurse interview
# erectile dysfunction/impotence 1518
ctt$erectile_disfunction <- 0
for (j in 0:33) { # gather all records from multiple arrays
  text2 <- paste0("f.20002.0.", j)
  ctt$erectile_disfunction[bd[[text2]]== 1518] <- 1
}

# refer to drug only also
# reference: Alice Carter. et al. Education inequalities in statin treatment
ed <- c(1140869100, 1140883010, 1141168936, 1141168944, 1141168946, 1141168948, 
        1141187810, 1141187814, 1141187818, 1141192248, 1141192256, 1141192258, 1141192260)
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt$erectile_disfunction[bd[[text2]] %in% ed] <- 1
}

# rheumatoid arthritis
# use first occurence data
# ICD: M05, M06
ctt$rheumatoid_arthritis <- ifelse(bd$f.131848.0.0<ctt$recruit.date | bd$f.131850.0.0<ctt$recruit.date, 1, 0)
ctt$rheumatoid_arthritis[is.na(ctt$rheumatoid_arthritis)] <- 0

# migraine
# use first occurrence data
# G43: bd$f.131052.0.0

ctt$migraine <- ifelse(bd$f.131052.0.0<ctt$recruit.date, 1, 0)
ctt$migraine[is.na(ctt$migraine)] <- 0

# systemic lupus erythematosis
# M32: bd$f.131894.0.0
ctt$systemic_lupus_erythematosis <- ifelse(bd$f.131894.0.0<ctt$recruit.date, 1, 0)
ctt$systemic_lupus_erythematosis[is.na(ctt$systemic_lupus_erythematosis)] <- 0

# severe mental illness
# F20: bd$f.130874.0.0
# F23: bd$f.130880.0.0
# F31: bd$f.130892.0.0
# F32: bd$f.130894.0.0
# F33: bd$f.130896.0.0

ctt$severe_mental_illness <- ifelse(bd$f.130874.0.0<ctt$recruit.date |
                                       bd$f.130880.0.0<ctt$recruit.date |
                                       bd$f.130892.0.0<ctt$recruit.date |
                                       bd$f.130894.0.0<ctt$recruit.date |
                                       bd$f.130896.0.0<ctt$recruit.date, 1, 0)

ctt$severe_mental_illness[is.na(ctt$severe_mental_illness)] <- 0
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
ctt$chronic_kidney_disease <- ifelse(bd$f.132032.0.0<ctt$recruit.date, 1, 0)
ctt$chronic_kidney_disease[is.na(ctt$chronic_kidney_disease)] <- 0

# steroid
# medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment  
# import the code list 
# the first element is coded incorrectly
# revise it manually 
corti$V1[1] <- 1140853854

corti <- as.numeric(corti$V1)
ctt$regular_steroid_tablets <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt$regular_steroid_tablets[bd[[text2]] %in% corti] <- 1
}

# atypical psychotics
# medication from verbal interview
# reference: Alice Carter. et al. Education inequalities in statin treatment 
ap <- c(1140867420, 1140867432, 1140867444, 1140927956, 1140927970, 1140928916, 
        1141152848, 1141152860, 1141153490, 1141167976, 1141177762, 1141195974, 1141202024)

ctt$atypical_antipsy <- 0
for (j in 0:47) { # gather all records from multiple arrays
  text2 <- paste0("f.20003.0.", j)
  ctt$atypical_antipsy[bd[[text2]] %in% ap] <- 1
}


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
ctt[[text1]] <- NA
ctt[[text1]][bd[[text2]]%in% white] <- 1
ctt[[text1]][bd[[text2]]%in% indian] <- 2
ctt[[text1]][bd[[text2]]%in% pak] <- 3
ctt[[text1]][bd[[text2]]%in% bang] <- 4
ctt[[text1]][bd[[text2]]%in% oasia] <- 5
ctt[[text1]][bd[[text2]]%in% bc] <- 6
ctt[[text1]][bd[[text2]]%in% ba] <- 7
ctt[[text1]][bd[[text2]]%in% chi] <- 8
ctt[[text1]][bd[[text2]]%in% oeth] <- 9
# assume missing as white
ctt$ethnicity[is.na(ctt$ethnicity)] <- 1

# generate another ethnic vars
ctt$ethn3 <- ifelse(ctt$ethnicity ==1, "White", 
                     ifelse(ctt$ethnicity %in% c(6,7), "Black", 
                            ifelse(ctt$ethnicity %in% c(2,3,4), "South Asian", "Mixed or others")))

ctt$ethn3 <- factor(ctt$ethn3, levels=c("White", "Black", "South Asian", "Mixed or others"))
# 

# end ---------------------------------------------------------------------


# for diabetes risk -------------------------------------------------------

# refer to qdiabetes.org

# father, mother, siblings diabetes
ctt$diabetes_relative <- 0
for (j in 0:9) { # father illness
  text2 <- paste0("f.20107.0.", j)
  ctt$diabetes_relative[bd[[text2]]== "Diabetes"] <- 1
}

for (j in 0:10) { # mother illness
  text2 <- paste0("f.20110.0.", j)
  ctt$diabetes_relative[bd[[text2]]== "Diabetes"] <- 1
}

for (j in 0:11) { # sibling illness
  text2 <- paste0("f.20111.0.", j)
  ctt$diabetes_relative[bd[[text2]]== "Diabetes"] <- 1
}

# # copy to ctt2
ctt2 <- merge(ctt2, ctt[, c("id", "diabetes_relative")])
saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"), compress = F)


# End stage renal disease, esrd, baseline ----
# ctt <- readRDS(file.path(work_data, "ctt2020.rds"))
# algorithmically ESRD before recruitment
ctt$esrd <- ifelse(as.Date(bd$f.42026.0.0)<=as.Date(bd$f.53.0.0), 1, 0)
ctt$esrd[is.na(ctt$esrd)] <- 0
ctt$esrd.date <- as.Date(bd$f.42026.0.0)

ctt2 <- merge(ctt2, ctt[, c("id", "esrd", "esrd.date")])
saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"), compress = F)

# save ----
saveRDS(ctt, file = file.path(work_data, "ctt2020.rds"), compress = F)


