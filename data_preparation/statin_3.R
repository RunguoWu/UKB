library(tidyverse)

source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

gp_clin <- read.delim(file.path(raw_data, "gp_clinical.txt"))

gp_scr <- read.delim(file.path(raw_data, "gp_scripts.txt"))

ctt <- readRDS(file.path(work_data, "ctt2020.rds"))

ctt2 <- readRDS(file.path(work_data, "ctt2.rds"))

# statin in GP record code
# exclude ciprofibrate 
# exclude cerivastatin/baycol/lipobay

# read in prepared statin read2:drug name table
statin_list <- read.csv(file.path(work_data, "statin_read2_name.csv"))
sti_read2 <- paste(statin_list$read2, collapse = "|")
sti_name <- paste(statin_list[statin_list$match!="", ]$match, collapse = "|")

# extract lines with use of statin record, either read2 or drug name 
gp_statin <- gp_scr %>% filter(grepl(sti_read2, read_2) | 
                               grepl(sti_name, drug_name, ignore.case=T)) 


temp <- gp_statin %>% select(eid, issue_date, read_2, drug_name)
# append the individuals' entry date 
temp <- merge(temp, ctt[, c("id", "recruit.date")], by.x="eid", by.y="id")
temp <- transform(temp, statin_date=as.Date(as.character(issue_date), "%d/%m/%Y"))


# deal with drug name -----------------------------------------------------

# extract all numbers in drug_name strings into a list for each record
temp$dose <- regmatches(temp$drug_name, gregexpr("[[:digit:]]+", temp$drug_name))

# as a result, a list can include up to 3 numbers,  
# dose3 only includes "0" value, not a valid value
# only dose 1 and 2 are valid
temp$dose1 <- sapply(temp$dose, "[", 1)
temp$dose2 <- sapply(temp$dose, "[", 2)
temp$dose3 <- sapply(temp$dose, "[", 3)

# assume the true dosage is dose1 first
# this suits most cases
temp$statin_dose <- as.numeric(temp$dose1)

# tab dose1 values and check corresponding drug name
# handle inappropriate values
# 4 cases to correct

# 1. those corresponding to non-statin/eze drug name
temp$statin_dose[temp$statin_dose %in% c(2010, 2020, 4750) ] <- NA

# 2. change to dose2, when the two first numbers are irrelevant
temp$statin_dose[temp$statin_dose %in% c(24, 42) ] <- 
  as.numeric(temp$dose2[temp$statin_dose %in% c(24, 42)])

# 3. change to dose2, when eze+statin and eze dose goes first
eze_first <- c("Inegy 10mg/40mg tablets (Merck Sharp & Dohme Ltd)",
           "Inegy 10mg/20mg tablets (Merck Sharp & Dohme Ltd)",
           "INEGY tabs 10mg + 40mg", 
           "EZETIMIBE + SIMVASTATIN tabs 10mg + 40mg",
           "Ezetimibe And Simvastatin Tablets 10 mg + 20 mg",
           "Inegy  Tablets  10 mg/40 mg",
           "EZETIMIBE + SIMVASTATIN tabs 10mg + 20mg",
           "INEGY tabs 10mg + 20mg",
           "Inegy 10mg/80mg tablets (Merck Sharp & Dohme Ltd)",
           "Ezetimibe And Simvastatin  Tablets  10 mg + 20 mg",
           "Ezetimibe And Simvastatin  Tablets  10 mg + 40 mg",
           "Inegy  Tablets  10 mg/80 mg",
           "Inegy 10mg/40mg tablets (MSD-SP Ltd)")

temp$statin_dose[temp$drug_name %in% eze_first] <- 
  as.numeric(temp$dose2[temp$drug_name %in% eze_first])

# 4. recode those with ezetimibe only statin dose to NA
# corrected after define statin and ezetimibe names

# statin name
temp$statin <- NA

simva <- "SIMVASTATIN|ZOCOR|SIMVADOR|RANZOLONT|INEGY"
prava <- "PRAVASTATIN|LIPOSTAT"
fluva <- "FLUVASTATIN|LESCOL|LUVINSTA|STEFLUVIN|DORISIN|NANDOVAR"
atorva <- "ATORVASTATIN|LIPITOR"
rosuva <- "CRESTOR|ROSUVASTATIN"

temp$statin[grepl(simva, temp$drug_name, ignore.case=T)] <- "simva"
temp$statin[grepl(prava, temp$drug_name, ignore.case=T)] <- "prava"
temp$statin[grepl(fluva, temp$drug_name, ignore.case=T)] <- "fluva"
temp$statin[grepl(atorva, temp$drug_name, ignore.case=T)] <- "atorva"
temp$statin[grepl(rosuva, temp$drug_name, ignore.case=T)] <- "rosuva"

# ezetimibe name
temp$ezetimibe <- NA

ezeti <- "INEGY|EZETIMIBE|EZETROL"
  
temp$ezetimibe[grepl(ezeti, temp$drug_name, ignore.case=T)] <- "ezeti"


# deal with READ 2 --------------------------------------------------------

# name
temp$statin[grepl("^bxd", temp$read_2)] <- "simva"
temp$statin[grepl("^bxe", temp$read_2)] <- "prava"
temp$statin[grepl("^bxg", temp$read_2)] <- "fluva"
temp$statin[grepl("^bxi", temp$read_2)] <- "atorva"
temp$statin[grepl("^bxk", temp$read_2)] <- "rosuva"

temp$ezetimibe[grepl("^bxl|bxdH.|bxdI.|bxdJ.|bxdw.|bxdx.|bxdy.", temp$read_2)] <- "ezeti"

# dose

mg10 <- c(paste0("bxd", c("1","3","7","A", "D", "G"),"."),
         paste0("bxe", c("1", "3", "5"), "."), 
         paste0("bxi", c("1", "4", "8", "9"), "."), 
         "bxk1.", "bxkx.")

mg20 <- c(paste0("bxd", c("2","4","8","B", "E", "H", "u", "y"),"."),
          paste0("bxe", c("2", "4", "6"), "."),
          "bxg1.","bxg3.",
          paste0("bxi", c("2", "5", "A", "B"), "."), 
          "bxk2.", "bxky.")

mg40 <- c(paste0("bxd", c("5","6","C", "F", "I", "v", "x"),"."),
          paste0("bxe", c("7", "8"), "."),
          "bxg2.", "bxg4.",
          paste0("bxi", c("3", "6"), "."), 
          "bxk3.", "bxkz.")


mg80 <- c(paste0("bxd", c("9","J","K", "w", "z"),"."),
          paste0("bxg", c("5","6", "7", "8", "9", "z"), "."),
          paste0("bxi", c("7", "z"), "."))

mg60 <- c("bxiy.")

mg5 <- c("bxk4.", "bxkw.")

temp$statin_dose[temp$read_2 %in% mg5] <- 5
temp$statin_dose[temp$read_2 %in% mg10] <- 10
temp$statin_dose[temp$read_2 %in% mg20] <- 20
temp$statin_dose[temp$read_2 %in% mg40] <- 40
temp$statin_dose[temp$read_2 %in% mg60] <- 60
temp$statin_dose[temp$read_2 %in% mg80] <- 80


# 4. finally, recode those with ezetimibe only statin dose to NA
temp$statin_dose[is.na(temp$statin) & temp$ezetimibe=="ezeti"] <- NA



# keep the relevant records -----------------------------------------------

# # only keep within 14 months before recruitment and 6 months after recruitment
# temp2 <-  temp %>% filter(statin_date <= recruit.date + 183 &
#                           statin_date >= recruit.date - 427) %>% 
#   select(eid,statin_date, statin, ezetimibe, statin_dose) %>% 
#   rename(issue_date=statin_date)

# the above is a test, now use 12 months before recruitment ONLY
temp2 <-  temp %>% filter(statin_date <= recruit.date &
                            statin_date >= recruit.date - 365.25) %>% 
  select(eid,statin_date, statin, ezetimibe, statin_dose) %>% 
  rename(issue_date=statin_date)

# for ones with multiple records only keep the latest record(s)
# if there are multiple records in the latest issue date
# remove the identical regimen records first

temp2 <- temp2 %>% 
  group_by(eid) %>% 
  filter(issue_date==max(issue_date)) %>% 
  distinct(eid, statin, statin_dose, ezetimibe, .keep_all = T )

# obtain the ids still with duplicated record
id_dup <- temp2[duplicated(temp2$eid),]$eid

# combine stain, ezetimibe names and dosage to one column
temp2$statin[is.na(temp2$statin)] <- ""
temp2$ezetimibe[is.na(temp2$ezetimibe)] <- ""
temp2$statin_dose[is.na(temp2$statin_dose)] <- ""

# deal with duplicated ids
temp_dup <- temp2 %>% filter(eid %in% id_dup) %>% 
  mutate(statin_regimen=paste0(statin," ", statin_dose)) %>% # combine statin and dosage first
  mutate(statin_regimen=if_else(statin_regimen==" ", "", statin_regimen)) %>% 
  group_by(eid) %>% 
  arrange(desc(statin_regimen), .by_group=T) %>% # let ezeti placed always after statin
  mutate(comb_regimen = paste(statin_regimen, ezetimibe, collapse = "+")) %>% # combine all stain and ezeti by id
  mutate(comb_regimen = gsub(" ", "", comb_regimen, fixed = T)) %>%  # remove all whitespace 
  distinct(eid, .keep_all = T) %>% # keep distinct ids 
  select(eid, issue_date, comb_regimen) %>% 
  rename(statin_issue_date=issue_date, 
         statin_regimen=comb_regimen)

# deal with no duplicated ids
temp_nodup <- temp2 %>% filter(!eid %in% id_dup) %>% 
  mutate(comb_regimen=paste0(statin,statin_dose), 
         comb_regimen=if_else(ezetimibe=="ezeti" & statin=="","ezeti", comb_regimen), 
         comb_regimen=if_else(ezetimibe=="ezeti" & statin!="", paste0(comb_regimen,"+ezeti"), comb_regimen)) %>% 
  select(eid, issue_date, comb_regimen) %>% 
  rename(statin_issue_date=issue_date, 
         statin_regimen=comb_regimen)

temp3 <- rbind(temp_dup, temp_nodup)


# merge to ctt
ctt <- merge(ctt, temp3, by.x = "id", by.y = "eid", all.x = T)
saveRDS(ctt, file = file.path(work_data, "ctt2020.rds"), compress = F)

# merge to ctt2
ctt2 <- merge(ctt2, temp3, by.x = "id", by.y = "eid", all.x = T)
saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"))



# gp registration data ----------------------------------------------------

gp_reg <- read.delim(file.path(raw_data, "gp_registrations.txt"))

gp_reg <- transform(gp_reg, reg.date=as.Date(as.character(reg_date), "%d/%m/%Y"))

gp_reg <- transform(gp_reg, remove.date=as.Date(as.character(deduct_date), "%d/%m/%Y"))

# code wrong dates to NA
gp_reg$remove.date[gp_reg$remove.date == "1901-01-01"] <- NA
gp_reg$remove.date[gp_reg$remove.date == "1902-02-02"] <- NA
gp_reg$remove.date[gp_reg$remove.date == "1903-03-03"] <- NA
gp_reg$remove.date[gp_reg$remove.date == "2037-07-07"] <- NA                   

gp_reg$reg.date[gp_reg$reg.date == "1901-01-01"] <- NA
gp_reg$reg.date[gp_reg$reg.date == "1902-02-02"] <- NA
gp_reg$reg.date[gp_reg$reg.date == "1903-03-03"] <- NA
gp_reg$reg.date[gp_reg$reg.date == "2037-07-07"] <- NA


# with/out gp registration record
gp <- data.frame(id=unique(gp_reg$eid))
gp$gp_record_reg <- 1

ctt <- left_join(ctt, gp, by = "id")
ctt$gp_record_reg[is.na(ctt$gp_record_reg)] <- 0

ctt2 <- left_join(ctt2, gp, by = "id")
ctt2$gp_record_reg[is.na(ctt2$gp_record_reg)] <- 0


# # remove those who were removed before 6 months after recruitment
# gp_reg2 <- merge(gp_reg, ctt[, c("id","recruit.date")], by.x = "eid", by.y = "id", all.x = T )
# gp_reg2 <- gp_reg2 %>% filter(is.na(remove.date) | remove.date >= recruit.date + 183)

# above is a test, use before recruitment
gp_reg2 <- merge(gp_reg, ctt[, c("id","recruit.date")], by.x = "eid", by.y = "id", all.x = T )
gp_reg2 <- gp_reg2 %>% filter(is.na(remove.date) | remove.date >= recruit.date)

gp <- data.frame(id=unique(gp_reg2$eid))
gp$gp_record_reg2 <- 1

ctt <- left_join(ctt, gp, by = "id")
ctt$gp_record_reg2[is.na(ctt$gp_record_reg2)] <- 0

ctt2 <- left_join(ctt2, gp, by = "id")
ctt2$gp_record_reg2[is.na(ctt2$gp_record_reg2)] <- 0


# with/out gp prescription records
gp <- data.frame(id=unique(gp_scr$eid))
gp$gp_record_scr <- 1

ctt <- left_join(ctt, gp, by = "id")
ctt$gp_record_scr[is.na(ctt$gp_record_scr)] <- 0

ctt2 <- left_join(ctt2, gp, by = "id")
ctt2$gp_record_scr[is.na(ctt2$gp_record_scr)] <- 0

saveRDS(ctt, file = file.path(work_data, "ctt2020.rds"), compress = F)
saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"))



# statin effect - based on regimen + type ---------------------------------

# first, use gp prescription data----
# export statin regimen categories 
# statin_regimen is the primary data source, because it has dosage info
tab <- table(ctt2$statin_regimen, useNA = "ifany")
tab
# expo the table to generate the dictionary
write.csv(tab, file = file.path(output, "statin_regimen.csv"))

# recode conditions of two statin together
# only keep the most effective one
# one case simva 40 = atorva 10 in %LDL reduction, keep atorva 
# read in a dictionary 
statin_dic <- read.csv(file.path(work_data, "statin_regimen_effect_dic.csv"))

# merge into ctt2
# use temp, because by.x = "statin_regimen", by.y = "regimen" will reorder the rows
# we'd better keep the id order all along
temp <- ctt2[, c("id", "statin_regimen")]

temp <- merge(temp, statin_dic[, c("regimen", "revised_regimen")], 
              by.x = "statin_regimen", by.y = "regimen", all.x = T)

ctt2$statin_regimen <- NULL

ctt2 <- merge(ctt2, temp, by = "id") # by id


# read in treatment effect file for %LDL reduction
ldl_red <- readRDS(file.path(vali_data, "tx.rds"))
ctt2$LDL_red <- ldl_red$statin_ldl_red[ctt2$revised_regimen]



# second, use nurse interview data----
# for those without gp statin prescription
# impute nurse interview result
ctt2$revised_regimen[is.na(ctt2$statin_regimen)] <- 
                  ctt2$statin_type[is.na(ctt2$statin_regimen)]


# recode conditions of two statin together
# rosuva > atorva > simva > prava > fluva
# simva > prava > fluva
# simva > fluva

ctt2$revised_regimen[ctt2$revised_regimen=="atorva+simva"] <- "atorva"
ctt2$revised_regimen[ctt2$revised_regimen=="simva+rosuva"] <- "rosuva"
ctt2$revised_regimen[ctt2$revised_regimen=="atorva+rosuva"] <- "rosuva"
ctt2$revised_regimen[ctt2$revised_regimen=="atorva+simva+ezeti"] <- "atorva+ezeti"
ctt2$revised_regimen[ctt2$revised_regimen=="simva+prava"] <- "simva"
ctt2$revised_regimen[ctt2$revised_regimen=="atorva+rosuva+ezeti"] <- "rosuva+ezeti"
ctt2$revised_regimen[ctt2$revised_regimen=="simva+fluva"] <- "simva"


# for no dosage information
# assume average reduction according to weighted frequency
# by pp/sp
# calculate them in excel
ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="atorva"] <- 0.429 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="atorva"] <- 0.45 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="atorva+ezeti"] <- 0.563 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="atorva+ezeti"] <- 0.573 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="fluva"] <- 0.268 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="fluva"] <- 0.264 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="fluva+ezeti"] <- 0.386 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="fluva+ezeti"] <- 0.375

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="prava"] <- 0.251 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="prava"] <- 0.261 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="prava+ezeti"] <- 0.389 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="prava+ezeti"] <- 0.394

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="rosuva"] <- 0.435 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="rosuva"] <- 0.443 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="rosuva+ezeti"] <- 0.613 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="rosuva+ezeti"] <- 0.612

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="simva"] <- 0.350 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="simva"] <- 0.355 

ctt2$LDL_red[ctt2$CVhist==0 & ctt2$revised_regimen=="simva+ezeti"] <- 0.556 
ctt2$LDL_red[ctt2$CVhist>0 & ctt2$revised_regimen=="simva+ezeti"] <- 0.563

ctt2$LDL_red[ctt2$revised_regimen=="ezeti"] <- 0.185

# replace the values for those with verbal interview statin only
# but also have dosage information
# with exact %LDL reduction
ctt2$LDL_red[ctt2$revised_regimen=="atorva" & ctt2$statin_dose_vi=="atorva10"] <- 0.37
ctt2$LDL_red[ctt2$revised_regimen=="atorva+ezeti" & ctt2$statin_dose_vi=="atorva10"] <- 0.53

ctt2$LDL_red[ctt2$revised_regimen=="fluva" & ctt2$statin_dose_vi=="fluva20"] <- 0.21
ctt2$LDL_red[ctt2$revised_regimen=="fluva+ezeti" & ctt2$statin_dose_vi=="fluva20"] <- 0.32

ctt2$LDL_red[ctt2$revised_regimen=="prava" & ctt2$statin_dose_vi=="prava10"] <- 0.20
ctt2$LDL_red[ctt2$revised_regimen=="prava+ezeti" & ctt2$statin_dose_vi=="prava10"] <- 0.34

ctt2$LDL_red[ctt2$revised_regimen=="rosuva" & ctt2$statin_dose_vi=="rosuva10"] <- 0.43
ctt2$LDL_red[ctt2$revised_regimen=="rosuva+ezeti" & ctt2$statin_dose_vi=="rosuva10"] <- 0.597

ctt2$LDL_red[ctt2$revised_regimen=="simva" & ctt2$statin_dose_vi=="simva10"] <- 0.27
ctt2$LDL_red[ctt2$revised_regimen=="simva+ezeti" & ctt2$statin_dose_vi=="simva10"] <- 0.471


# finally, recode LDL reduction missing =  zero
ctt2$LDL_red[is.na(ctt2$LDL_red)] <- 0

# calculate theoretical LDL concentration without statin effect
ctt2$LDL_nostatin <- ctt2$LDL.imp/(1-ctt2$LDL_red)

# add LDL categories
ctt2 <- ctt2 %>% 
  mutate(ldl_cat = case_when(LDL_nostatin < 3.4 ~ "<3.4", 
                             LDL_nostatin >= 3.4 & LDL_nostatin < 4.1~ "3.4-4.1",
                             LDL_nostatin >= 4.1 ~ ">=4.1"))

ctt2$ldl_cat <- factor(ctt2$ldl_cat, levels = c("<3.4", "3.4-4.1", ">=4.1"))

saveRDS(ctt2, file = file.path(work_data, "ctt2.rds"))









