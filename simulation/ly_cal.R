f_sum_lys_qaly <- function(rst_i){
  
  tot_d <- rst_i %>% 
    select(id_new, cycle, vd, nvd, qol) %>% 
    mutate(d = vd + nvd) %>%  
    select(id_new, cycle, d, qol)
  
  tot_d1 <- tot_d %>% group_by(id_new) %>% filter(cycle != max(cycle))
  
  tot_d2 <- tot_d1 %>%  
    group_by(id_new) %>% 
    summarise(cycle = max(cycle), d = sum(d)) %>% 
    mutate(cycle = cycle + 1, d = 1 - d)
  
  tot_d2_qol <- tot_d %>% group_by(id_new) %>% 
    filter(cycle == max(cycle)) %>% select(id_new, qol)
  
  tot_d2 <- left_join(tot_d2, tot_d2_qol)
  
  tot_lys <- bind_rows(tot_d1, tot_d2) %>% 
    arrange(id_new, cycle) %>% 
    mutate(ly = d * (cycle-0.5)) %>% 
    group_by(id_new) %>% 
    summarise(ly = sum(ly), qaly = sum(qol))
  output <- tot_lys 
  return(output)
}
