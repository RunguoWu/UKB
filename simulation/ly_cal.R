f_sum_lys_qaly <- function(rst_i, years){
  
  tot_d <- rst_i %>% 
    select(id_new, cycle, vd, nvd, qol, mi, stroke,crv, dm, cancer_icd) %>% 
    mutate(d = vd + nvd) %>%  
    select(id_new, cycle, d, qol, mi, stroke, crv, dm, cancer_icd, vd, nvd)
  
  tot_d1 <- tot_d %>% group_by(id_new) %>% filter(cycle != max(cycle)) %>% 
    select(id_new, cycle, d, qol)
  
  tot_d2 <- tot_d1 %>%  
    # to adjust the last cycle death rate, 
    # to make sure 100% cumulative death rate
    group_by(id_new) %>% 
    summarise(cycle = max(cycle), d = sum(d)) %>% 
    mutate(cycle = cycle + 1, d = 1 - d)
  
  tot_d2_qol <- tot_d %>% 
    group_by(id_new) %>% 
    filter(cycle == max(cycle)) %>% 
    select(id_new, qol)
  
  tot_d2 <- left_join(tot_d2, tot_d2_qol)
  
  tot_lys <- bind_rows(tot_d1, tot_d2) %>% 
    arrange(id_new, cycle) %>% 
    mutate(ly = d * (cycle-0.5)) %>% 
    # d indicates death rate in the particular cycle, 
    # in other words, survival probability up to this cycle
    # for one individual, sum(d)=1
    # so the ly equals to sum of product of each possible year (cycle) and its probability (d) 
    group_by(id_new) %>% 
    summarise(ly = sum(ly), qaly = sum(qol))
  
  tot_e <- tot_d %>% 
    arrange(id_new, cycle) %>% 
    group_by(id_new) %>% 
    filter(cycle<=years) %>% 
    summarise(mi = sum(mi), stroke = sum(stroke), crv = sum(crv),
  dm = sum(dm), cancer = sum(cancer_icd), vd = sum (vd), nvd = sum(nvd))
  
  output <- merge(tot_lys, tot_e) 
  
  return(output)
} 


f_sum <- function(rst_i, years){

  tot_lys <- rst_i %>% 
    group_by(id_new) %>% 
    summarise(ly = sum(ly), qaly = sum(qol))
  
  tot_e <- rst_i %>% 
    group_by(id_new) %>% 
    filter(cycle<=years) %>% 
    summarise(mi = sum(mi), stroke = sum(stroke), crv = sum(crv),
              dm = sum(dm), cancer = sum(cancer_icd), vd = sum (vd), 
              nvd = sum(nvd), mvevd=sum(mvevd))
  
  output <- merge(tot_lys, tot_e) 
  
  return(output)
} 

# summarise number and age, QRISK and LDL
# by var1 var2
sum_sample <- function(dt, var1, var2=NULL) {
  
  if (!is.null(var2)) {
    a1 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(age=mean(age.recruit)) %>% spread(2, age, sep='') # 2 means 2nd column,i.e. var2
    
    a2 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(sd=sd(age.recruit)) %>% spread(2, sd, sep='')
    
    b1 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(risk=mean(QRISK3_2017)) %>% spread(2, risk, sep='')
    
    b2 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(sd=sd(QRISK3_2017)) %>% spread(2, sd, sep='')
    
    c1 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(ldl=mean(LDL_nostatin)) %>% spread(2, ldl, sep='')
    
    c2 <- dt %>% group_by(get(var1), get(var2)) %>% 
      summarise(sd=sd(LDL_nostatin)) %>% spread(2, sd, sep='')
    
    d <- dt %>% group_by(get(var1), get(var2)) %>% count() %>% spread(2,n)
    
  } else {
    
    a1 <- dt %>% group_by(get(var1)) %>% summarise(age=mean(age.recruit)) 
    
    a2 <- dt %>% group_by(get(var1)) %>% summarise(sd=sd(age.recruit)) 
    b1 <- dt %>% group_by(get(var1)) %>% summarise(risk=mean(QRISK3_2017)) 
    
    b2 <- dt %>% group_by(get(var1)) %>% summarise(sd=sd(QRISK3_2017)) 
    
    c1 <- dt %>% group_by(get(var1)) %>% summarise(ldl=mean(LDL_nostatin)) 
    
    c2 <- dt %>% group_by(get(var1)) %>% summarise(sd=sd(LDL_nostatin)) 
    
    d <- dt %>% group_by(get(var1)) %>% count()
  }
  
  a1 <- round(a1[,-1] , digits = 2)
  a2 <- round(a2[,-1] , digits = 2)
  b1 <- round(b1[,-1] , digits = 2)
  b2 <- round(b2[,-1] , digits = 2)
  c1 <- round(c1[,-1] , digits = 2)
  c2 <- round(c2[,-1] , digits = 2)
  
  N <- nrow(d)
  
  e <- c()
  cat <- c()
  for (i in 1:N) {
    add <- rbind(d[i,-1],
                 paste0(a1[i,], " (", a2[i,], ")"),
                 paste0(b1[i,], " (", b2[i,], ")"),
                 paste0(c1[i,], " (", c2[i,], ")"))
    e <- rbind(e, add)
    
    nam <- c(as.character(d[i, 1]), "Age","QRISK", "LDL")
    
    cat <-c(cat, nam)
  }
  
  expo <- cbind(cat,e)
  
  return(expo)
  
}



