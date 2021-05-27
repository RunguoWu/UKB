f_sum_lys_qaly <- function(rst_i){
  
  tot_d <- rst_i %>% 
    select(id_new, cycle, vd, nvd, qol, mi, stroke, dm, cancer_icd) %>% 
    mutate(d = vd + nvd) %>%  
    select(id_new, cycle, d, qol, mi, stroke, dm, cancer_icd, vd, nvd)
  
  tot_d1 <- tot_d %>% group_by(id_new) %>% filter(cycle != max(cycle))
  
  tot_d2 <- tot_d1 %>%  
    # to adjust the last cycle death rate, 
    # to make sure 100% cumulative death rate
    group_by(id_new) %>% 
    summarise(cycle = max(cycle), d = sum(d)) %>% 
    mutate(cycle = cycle + 1, d = 1 - d)
  
  tot_d2_qol <- tot_d %>% group_by(id_new) %>% 
    filter(cycle == max(cycle)) %>% 
    select(id_new, qol, mi, stroke, dm, cancer_icd, vd, nvd)
  
  tot_d2 <- left_join(tot_d2, tot_d2_qol)
  
  tot_lys <- bind_rows(tot_d1, tot_d2) %>% 
    arrange(id_new, cycle) %>% 
    mutate(ly = d * (cycle-0.5)) %>% 
    # d indicates death rate in the particular cycle, 
    # in other words, survival probability up to this cycle
    # for one individual, sum(d)=1
    # so the ly equals to sum of product of each possible year (cycle) and its probability (d) 
    group_by(id_new) %>% 
    summarise(ly = sum(ly), qaly = sum(qol), mi = sum(mi), stroke = sum(stroke), 
              dm = sum(dm), cancer = sum(cancer_icd), vd = sum (vd), nvd = sum(nvd))
  output <- tot_lys 
  return(output)
}



# # secondary prevention
# ukb_sim_sec <- 
#   readRDS(file.path(vali_output, "ukb_vali_sim_summarised_sec.rds"))
# 
# ly_sec <- f_sum_lys_qaly(ukb_sim_sec)
# 
# sec_id_male <- as.data.frame(ukb_b_sec) %>% select(id_new, male)
# sec_id_age <- as.data.frame(ukb_t_sec) %>% select(id_new,CurrAge_cent)
# 
# sec_id <- merge(sec_id_age, sec_id_male)
# 
# sec_id$age3 <- ifelse(sec_id$CurrAge_cent < -1, "<50", 
#                    ifelse(sec_id$CurrAge_cent >= -1 & sec_id$CurrAge_cent <= -0.1, "50-59", "60+")) 
# 
# ly_sec <- merge(ly_sec, sec_id)
# 
# save(ly_sec, file=file.path(vali_output, "ly_sec.Rdata"))
# 
# # primary prevention
# 
# ukb_sim_prim <- 
#   readRDS(file.path(vali_output, "ukb_vali_sim_summarised_prim.rds"))
# 
# ly_prim <- f_sum_lys_qaly(ukb_sim_prim)
# 
# prim_id_male <- as.data.frame(ukb_b_prim) %>% select(id_new, male)
# 
# prim_id_age <- as.data.frame(ukb_t_prim) %>% select(id_new,CurrAge_cent)
# 
# prim_id <- merge(prim_id_age, prim_id_male)
# 
# prim_id$age3 <- ifelse(prim_id$CurrAge_cent < -1, "<50", 
#                       ifelse(prim_id$CurrAge_cent >= -1 & prim_id$CurrAge_cent <= -0.1, "50-59", "60+")) 
# 
# ly_prim <- merge(ly_prim, prim_id)
# 
# ly_prim <- merge(ly_prim, ukb_base_prim[,c("id_new","RG5")])
# 
# save(ly_prim, file=file.path(vali_output, "ly_prim.Rdata"))
# 
# # explore
# # ly_prim$RG5 <- factor(ly_prim$RG5, levels = c("<5", "[5,10)", "[10,15)", "[15,20)", ">=20"))
# 
# ly_prim %>% filter(male==1) %>% group_by(age3, RG5) %>% summarise(mean(ly), mean(qaly))
# 
# ly_prim %>% filter(male==0) %>% group_by(age3, RG5) %>% summarise(mean(ly), mean(qaly))
# 
# ly_sec %>% group_by(age3, male) %>% summarise(mean(ly), mean(qaly))
# 
# 
# 
# 
# 
