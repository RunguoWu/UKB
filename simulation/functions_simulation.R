
### Survival equations -------------------------------------

### cumulative hazards

# TODO: re-write this using a "calibrated_eqns" argument & function 
# so that only need to pass one parameter

# exponential
cumhaz_exp <- function(lam_b, cf_t, x_t, t, shape,
                       scaleby, intercept_calibrated_i) {
  lam_old <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  log_lam <- log(lam_old) * scaleby + intercept_calibrated_i
  return(exp(log_lam) * t)
}

# weibull
cumhaz_wei <- function(lam_b, cf_t, x_t, t, shape,
                       scaleby, intercept_calibrated_i) {
  lam_old <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  log_lam <- log(lam_old) * scaleby + intercept_calibrated_i
  return(exp(log_lam) * (t^shape))
}

# gompertz
cumhaz_gom <- function(lam_b, cf_t, x_t, t, shape,
                       scaleby, intercept_calibrated_i) {
  lam_old <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  log_lam <- log(lam_old) * scaleby + intercept_calibrated_i
  return((exp(log_lam) / shape) * (exp(shape * t) - 1))
}

# cox 
cumhaz_cox <- function(lam_b, cf_t, x_t, t, shape,
                       scaleby, intercept_calibrated_i) {
  # shape is baseline hazard here
  lam_old <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  log_lam <- log(lam_old) * scaleby + intercept_calibrated_i
  return(log_lam * shape[t + 1])
}

### probability of an event

# TODO: this is one of the most time-consuming parts
# TODO: re-write as one common cumhaz function instead of calculating cumhaz_0 & cumhaz_1
# TODO: remove the get() part; refer to a list of functions instead
p_event <- function(dist, lam_b, cf_t, x_t, j, shape = NULL,
                    scaleby = 1, intercept_calibrated_i = 0) {
  cumhaz_fn <- str_c("cumhaz_", dist)
  cumhaz_0 <- get(cumhaz_fn)(lam_b = lam_b, cf_t = cf_t, x_t = x_t, t = (j - 1), shape = shape,
                             scaleby = scaleby, intercept_calibrated_i = intercept_calibrated_i)
  cumhaz_1 <- get(cumhaz_fn)(lam_b = lam_b, cf_t = cf_t, x_t = x_t, t = j, shape = shape,
                             scaleby = scaleby, intercept_calibrated_i = intercept_calibrated_i) 
  return(1 - exp(cumhaz_0 - cumhaz_1))
}


gen_next_cycle <- function(prim_flag, j, v_j, 
                           v_b_int, int_names_b,
                           nonlinage, v_b_int_nonlinage, int_names_b_nonlinage,
                           p_crv, 
                           events_nf_to_predict, events_nf_to_update, events_f, 
                           vars_t_eb, vars_t_e, eventsb,
                           t_since_eb_i, t_since_e_i,
                           alive,
                           dist, lam_b_i, cf_t, shape,
                           scaleby, intercept_calibrated_i,
                           use_cost = FALSE, 
                           side_effects,
                           had_rhabdo, had_myop, statin_high) {
  
  
  ### time-updated characteristics ------------
  
  CurrAge <- v_j["CurrAge_cent"] # RW: don't update age at the beginning
  v_j["cycle"] <- j
  
  # interactions
  
  for (v_int in int_names_b)
    v_j[str_c("CurrAge_cent_int_", v_int)] <- CurrAge * v_b_int[v_int]
  
  if (nonlinage) {
    # TODO: record flexpoint as parameter
    if (CurrAge < 1) {
      v_j["CurrAge_cent_1"] <- CurrAge
      v_j["CurrAge_cent_2"] <- 0
    } else {
      v_j["CurrAge_cent_1"] <- 1
      v_j["CurrAge_cent_2"] <- CurrAge - 1
    }
    for (v_int in int_names_b_nonlinage)
      for (gp in 1:2)
        v_j[str_c("CurrAge_cent_", gp, "_int_", v_int)] <- 
          v_j[str_c("CurrAge_cent_", gp)] * v_b_int_nonlinage[v_int]
  }
  
  # event indicators for those events that have happened within 0-3 years

  # update time-updated events
  for (e in events_nf_to_update){
    t <- t_since_e_i[[e]] + 1
    # check whether increase changes time bracket
    vars_t_e_temp <- vars_t_e[[e]]
    change_bracket <- (t %in% vars_t_e_temp)
    if (change_bracket) {
      t_right_ind <- min(which(vars_t_e_temp >= t))
      # varname of the previous time bracket
      varname_prev <- if (t_right_ind == 1)
        str_c(e, "_", 0, "_", vars_t_e_temp[t_right_ind]) else 
          str_c(e, "_", vars_t_e_temp[t_right_ind - 1], "_", vars_t_e_temp[t_right_ind])
      # varname of the new time bracket
      # use tolower so that Inf is transformed to inf
      varname_new <- tolower(str_c(e, "_", vars_t_e_temp[t_right_ind], "_", vars_t_e_temp[t_right_ind + 1]))
      # update v_j
      v_j[varname_prev] <- 0
      v_j[varname_new] <- 1
    }
    # update time since event
    t_since_e_i[[e]] <- t
  }
  # TODO: could remove e from events_nf_to_update once we are in the last bracket
  
  
  # update side effect parametre values
  or_dm_low <- side_effects$or_dm_low
  or_dm_high <- side_effects$or_dm_high
  p_myop <- side_effects$p_myop
  p_rhabdo <- side_effects$p_rhabdo
  p_nvd_rhabdo <- side_effects$p_nvd_rhabdo


  ### predict mi, stroke, crv & cancer -----------
  
  events_nf_rand <- sample(events_nf_to_predict)
  
  # if want CRV straight after MI
  if (!is.null(p_crv) & ("crv" %in% events_nf_rand) & ("mi" %in% events_nf_rand)) {
    # re-arrange if randomly draw a number >p_crv
    if (runif(1) > p_crv) {
      without_crv <- events_nf_rand[events_nf_rand != "crv"]
      without_crv_length <- length(without_crv)
      mi_n <- which(without_crv == "mi")
      events_nf_rand <- c(without_crv[1:mi_n], "crv")
      if (mi_n < without_crv_length)
        events_nf_rand <- c(events_nf_rand, without_crv[(mi_n + 1):without_crv_length])
    }
  }
  
  #e <- events_nf_rand[2]
  
  for (e in events_nf_rand) {
    
    # TODO check: if !is.na(t_since_e_i[[e]]), don't run p_event 
    
    p_rand_e <- runif(1)
    # p_rand_e <- 0
    p_e <- p_event(dist = dist[[e]], lam_b = lam_b_i[[e]], 
                   cf_t = cf_t[[e]], x_t = v_j, j = j, 
                   shape = shape[[e]],
                   scaleby = scaleby[[e]], 
                   intercept_calibrated_i = intercept_calibrated_i[[e]])
    
    # treatment side effect on incident diabetes - RW 2021-05-04
    if (e == "dm") {
      p_e <- p_e*or_dm_low/(1 - p_e + p_e*or_dm_low)
      p_e <- p_e*or_dm_high/(1 - p_e + p_e*or_dm_high)
    }
    
    if (p_e > p_rand_e) {
      events_nf_to_predict <- setdiff(events_nf_to_predict, e)
      # update time since event
      t_since_e_i[[e]] <- 0
      # update 0_t variable
      varname_0_1 <- str_c(e, "_0_", vars_t_e[[e]][1])
      # read off name of the 0_t variable
      v_j[varname_0_1] <- 1
    }
  }
  
  ### treatment side effects ------------
  # TODO: need to condition on being on Tx; not on having experienced rhabdo
  # no need Tx or had_rhabdo == 0 any more. 
  # we instead control values of p_myop p_rhabdo through the side_effect list outside cycle
  # if had_rhabdo > 0, side_effect list reset to default  
  myop <- (p_myop > runif(1))
  rhabdo <- (p_rhabdo > runif(1))
  
  if (myop & had_myop==0) had_myop <- 1
  
  ### predict fatal events --------------
  
  events_f_rand <- sample(events_f)
  
  #e <- "nvd"
  
  for (e in events_f_rand) {
    
    p_rand_e <- runif(1)
    p_e <- p_event(dist = dist[[e]], lam_b = lam_b_i[[e]], 
                   cf_t = cf_t[[e]], x_t = v_j, j = j, 
                   shape = shape[[e]],
                   scaleby = scaleby[[e]], 
                   intercept_calibrated_i = intercept_calibrated_i[[e]])
    # adjust if nvd
    if (rhabdo & had_rhabdo == 0) {
      if (e == "nvd") {
        p_e <- p_e + p_nvd_rhabdo
        had_rhabdo <- 1
      }
    }
    
    if (p_e > p_rand_e) {
      alive <- 0
      v_j[e] <- 1
      # update d_without_e indicator for events that have not happened
      for (e_pre_death in events_nf_to_predict)
        v_j[str_c("d_without_", e_pre_death)] <- 1
      break
    }
  }
  
  ### cost -----------------------------
  
  if (use_cost) {
    # OLS so just beta*X
    # this consists of baseline & time-updated parts
    cf_t_cost <- cf_t[["cost"]]
    cost <- lam_b_i[["cost"]] + v_j[names(cf_t_cost)] %*% cf_t_cost
    v_j["cost"] <- cost
  }
  
  #####
  # qol
  # RW 2021-03-17
  cf_t_qol <- cf_t[["qol"]]
  qol <- lam_b_i[["qol"]] + v_j[names(cf_t_qol)] %*% cf_t_qol
  # side effects
  if (myop) qol <- qol - 0.017*30/365.25 # 0.017 decrement over 30 days
  if (rhabdo) qol <- qol - 0.5*qol*7.5/365.25 - 0.2*qol*30/365.25 
  # 50% decrement for 7.5d followed by 20% decrement for 30d
  
  if (alive==0) qol <- qol*0.5
  v_j["qol"] <- qol
  #####
  
  
  
  # update time-varying baseline events: cancer, diabetes ---- 	
  
  # RW: as baseline event have been initially defined, we should update it after using it in prediction 
  # update time-updated baseline variables
  # TODO: perhaps re-write some bits as a function
  # as repeated later with time-updated events
  
  dic <- c("dm"=0, "cancer_bsl"=1) 	
  # baseline diabetes start as 0_5, while baseline cancer start at 1_2 	
  for (eb in eventsb){
    # does the timer need to be updated?
    # ie, was there an event?
    if (!is.na(t_since_eb_i[[eb]])) {
      t <- t_since_eb_i[[eb]] + 1
      # check whether increase changes time bracket
      vars_t_eb_temp <- vars_t_eb[[eb]]
      change_bracket <- (t %in% vars_t_eb_temp)
      if (change_bracket) {
        t_right_ind <- min(which(vars_t_eb_temp >= t))
        # varname of the previous time bracket
        varname_prev <- if (t_right_ind == 1)
          str_c(eb, "_", dic[eb], "_", vars_t_eb_temp[t_right_ind]) else # RW: 0 ----> dic[eb]
            str_c(eb, "_", vars_t_eb_temp[t_right_ind - 1], "_", vars_t_eb_temp[t_right_ind])
        # varname of the new time bracket
        # use tolower so that Inf is transformed to inf
        varname_new <- tolower(str_c(eb, "_", vars_t_eb_temp[t_right_ind], "_", vars_t_eb_temp[t_right_ind + 1]))
        # update v_j
        v_j[varname_prev] <- 0
        v_j[varname_new] <- 1
      }
      # update time since event
      t_since_eb_i[[eb]] <- t
    }
  }
  
  # RW: update age at the end of cycle
  v_j["CurrAge_cent"] <- CurrAge + 1 / 10
  
  # update study indicators & treatment flag at the end of year 1
  # from then on these will be the values to carry over
  # TODO: add flag for updating; eg not needed in lifetime simulation
  if (prim_flag == "sec" & j == 1) {
    v_j["ind6_1"] <- 0
    v_j["ind7_1"] <- 0
    v_j["wstatinayr2And"] <- v_j["wstatinayr1"]
    v_j["wstatinayr1"] <- 0
  }
  
  ### collate output & return
  ret_list <- list(v_j = v_j, 
                   events_nf_to_predict = events_nf_to_predict,  
                   alive = alive,
                   t_since_eb_i = t_since_eb_i, t_since_e_i = t_since_e_i,
                   had_rhabdo = had_rhabdo, had_myop = had_myop)
  
  
  return(ret_list)
}

###############################################################################
###############################################################################
###############################################################################

master <- function(.input_data_dir, 
                   cf_filename, 
                   use_cost = FALSE,
                   prim_flag, 
                   mx_b_filename_prefix,
                   mx_t_filename_prefix,
                   events_list,
                   pf_filename,
                   adjust_crv = TRUE, 
                   calibrated_eqns = FALSE,
                   nonlinage = FALSE,
                   tx = FALSE,
                   dist,
                   stop_expr, 
                   sample_pat = FALSE,
                   save_by_pat = FALSE,
                   n_sim,
                   n_cores, 
                   output_dir, 
                   output_filename_prefix,
                   regimen = "none",
                   side_effects_flag = FALSE,
                   side_dm = TRUE,
                   sens_multiplier = 1
){
  
  # set the seed
  set.seed(1234)
  
  
  # load data ---------------------------------------------------------------
  
  cf_all <- readRDS(file.path(.input_data_dir, 
                              str_c(cf_filename, ".rds")))
  mx_b <- readRDS(file.path(.input_data_dir, 
                            str_c(mx_b_filename_prefix, "_", prim_flag, ".rds")))
  mx_t <- readRDS(file.path(.input_data_dir, 
                            str_c(mx_t_filename_prefix, "_", prim_flag, ".rds")))
  
  pf <- readRDS(file.path(.input_data_dir,
                          str_c(pf_filename, ".rds")))
  
  p_crv <- if (adjust_crv == FALSE) NULL else 
                  pf[["p_crv"]][[prim_flag]]
  
  # events 
  events_nf <- events_list$events_nf
  events_f <- events_list$events_f
  events_all <- c(events_nf, events_f)
  events_to_predict <- pf[["events_to_predict"]][[prim_flag]]
    

  # data on time-updated covariates
  vars_t <- pf[["vars_t"]]

  
  # duration since baseline cancer/diabetes
  t_since_eb <- pf[["t_since_eb"]][[prim_flag]]
    

  # treatment effects and side effects
  if (!identical(tx, FALSE)) {
    tx <- readRDS(file.path(.input_data_dir, str_c(tx, ".rds")))
    # when tx is added, baseline LDL is replaced by pre-treated LDL 
    mx_b[,"NEWB_LDL_CL_cent"] <- mx_b[,"LDL_nostatin_cent"]
  }
  
  # prepare inputs ----------------------------------------------------------

  # identify interaction terms
  int_names_b <- sapply(strsplit(grep("_int_", colnames(mx_t), value = TRUE), "_int_"), "[", 2)
  int_names_b_nonlinage <- if (nonlinage)
    sapply(strsplit(grep("1_int_", colnames(mx_t), value = TRUE), "_int_"), "[", 2) else
      NULL
  # RW: 2021-03-29: avoid error when no interaction terms
  if (length(int_names_b_nonlinage)==0) int_names_b_nonlinage <- NULL
  

  # extract correct equations
  if (calibrated_eqns) {
    # TODO: this seems to only extract calibrated equations for secondary?
    # need to check
    cf_temp <- list( # RW
      mi = cf_all[["mi"]][["calibrated"]][[prim_flag]][[dist[["mi"]]]],
      stroke = cf_all[["stroke"]][["calibrated"]][[prim_flag]][[dist[["stroke"]]]],
      crv = cf_all[["crv"]][["calibrated"]][[prim_flag]][[dist[["crv"]]]],
      cancer_icd = cf_all[["cancer"]][["calibrated"]][[prim_flag]][[dist[["cancer_icd"]]]], # RW 
      vd = cf_all[["vd"]][["calibrated"]][[prim_flag]][[dist[["vd"]]]],
      nvd = cf_all[["nvd"]][["calibrated"]][[prim_flag]][[dist[["nvd"]]]], 
      dm = cf_all[["diabetes"]][[prim_flag]][[dist[["dm"]]]]
      ) 
  } else {
    cf_temp <- list( #RW
      mi = cf_all[["mi"]][[prim_flag]][[dist[["mi"]]]],
      stroke = cf_all[["stroke"]][[prim_flag]][[dist[["stroke"]]]],
      crv = cf_all[["crv"]][[prim_flag]][[dist[["crv"]]]],
      cancer_icd = cf_all[["cancer"]][[prim_flag]][[dist[["cancer_icd"]]]], # RW
      vd = cf_all[["vd"]][[prim_flag]][[dist[["vd"]]]],
      nvd = cf_all[["nvd"]][[prim_flag]][[dist[["nvd"]]]],
      dm = cf_all[["diabetes"]][[prim_flag]][[dist[["dm"]]]]
      )
  }
  
  # cost equation
  if (use_cost)
    cf_temp[["cost"]] <- cf_all[["cost"]]
  
  # QoL 
  cf_temp[["qol"]] <- cf_all[["qol"]][[prim_flag]]


  # baseline time-updated variables
  vars_t_eb <- vars_t[["b"]]
  eventsb <- names(vars_t_eb)
  # within-simulation updated variables
  vars_t_e <- vars_t[["e"]]

  
  # lam_b, cf_t and shape for all patients
  lam_b <- list()
  cf_t <- list()
  shape <- list()
  scaleby <- list()
  intercept_calibrated <- list()
  # e <- events_all[1]
  for (e in events_all) {
    cf_b <- cf_temp[[e]][["cf_b"]]
    lam_b[[e]] <- exp(mx_b[, names(cf_b)] %*% cf_b)
    cf_t[[e]] <- cf_temp[[e]][["cf_t"]]
    shape[[e]] <- cf_temp[[e]][["shape"]]
    
    # add calibrated parameters if relevant
    if (calibrated_eqns){
      scaleby[[e]] <- cf_temp[[e]][["scaleby"]]
      cf_intercept <- cf_temp[[e]][["intercept_calibrated"]]
      intercept_calibrated[[e]] <- mx_b[, "Intercept"] * cf_intercept
    } else {
      #TODO: this solution assumes that the calibration action is always performed
      # if calibrated_eqns == FALSE, the linear predictor LP is transformed as LP * 1 + 0
      scaleby[[e]] <- 1
      intercept_calibrated[[e]] <- rep(0, nrow(mx_b))
    }
    # diabetes is solely based on UKB, so no calibration
    scaleby[["dm"]] <- 1
    intercept_calibrated[["dm"]] <- rep(0, nrow(mx_b))
  }
  
  ##########################################
  # TODO
  # prepare tx effect here # RW 2021-05-05
  # treatment regimen
  # get the absolute LDL reduction
  ldl_red <- (mx_b[, "NEWB_LDL_CL_cent"] + 3.6) * tx$statin_ldl_red[regimen]
  mx_b <- cbind(mx_b, ldl_red)
  
  # treat effect
  noeffect <- rep(1, nrow(mx_b))
  tx_effect <- list("mi"=noeffect, "stroke"=noeffect , "crv"=noeffect, "vd"=noeffect)
  if (!identical(tx, FALSE)) {
    for (e in c("mi", "stroke", "crv", "vd")){
      # generate a variable "ldl_red" of LDL reduction in mmol/L
      # divide scaleby, because scale will be multiplied later in cycles
      tx_effect[[e]] <- exp(mx_b[, "ldl_red"] * tx$tx_effect[[e]]/cf_temp[[e]][["scaleby"]])
      lam_b[[e]] <- lam_b[[e]] * tx_effect[[e]]
    }
  }
  # high-intensity statin or not
  statin_high <- tx$statin_ldl_red[regimen]>=0.45
  ##########################################
  # treatment side effects
  side_effect_default <- list(or_dm_low = 1,
                              or_dm_high = 1,
                              p_myop = 0,
                              p_rhabdo = 0,
                              p_nvd_rhabdo = 0)
  
  side_effects <- side_effect_default

  if (regimen=="none") side_effects_flag <- FALSE
  
  if (side_effects_flag) side_effects <- tx$side_effects
  
  if (side_effects_flag & !statin_high) side_effects$or_dm_high <- 1
  
  if (side_effects_flag & !side_dm) {
    side_effects$or_dm_low <- 1
    side_effects$or_dm_high <- 1
  }
  
  side_effects$p_myop <- side_effects$p_myop * sens_multiplier
  side_effects$p_rhabdo <- side_effects$p_rhabdo * sens_multiplier
  ############################################# 
  
  
  # same for cost
  # TODO: no calibration etc incorporated here
  if (use_cost) {
    for (e in c("cost")) {
      cf_b <- cf_temp[[e]][["cf_b"]]
      lam_b[[e]] <- mx_b[, names(cf_b)] %*% cf_b
      cf_t[[e]] <- cf_temp[[e]][["cf_t"]]
    }
  }
  
  # QoL
  cf_b <- cf_temp[["qol"]][["cf_b"]]
  lam_b[["qol"]] <- mx_b[, names(cf_b)] %*% cf_b
  cf_t[["qol"]] <- cf_temp[["qol"]][["cf_t"]]


  
  # template output matrix
  # independent of each simulation
  out_names <- c(colnames(mx_t)) 
  N <- length(out_names)
  # +2 for qol and nsim
  # if age spline used, +2 for two extra age spline
  if (nonlinage) v_output_null <- matrix(nrow = 0, ncol = N + 4) else
    v_output_null <- matrix(nrow = 0, ncol = N + 2)
  

  vars_to_summarise_e <- str_c(names(vars_t_e), "_0_", unlist(lapply(vars_t_e, "[", 1)))
  # death variables
  vars_to_summarise_d <- grep("vd|without", out_names, value = TRUE)
  # combine
  vars_to_summarise <- c(vars_to_summarise_e, vars_to_summarise_d)
  
  # QoL
  vars_to_summarise <- c(vars_to_summarise, "qol")

  
  # vars_to_summarise <- grep("_0_1|dm_0_10|vd|without", out_names, value = TRUE) # RW: "0_1" will unexpectedly include cancer_10_15 etc
  
  # # TODO: replace previous line with the below if cost equation used
  if (use_cost)
    vars_to_summarise <- grep("0_1|vd|without|cost", out_names, value = TRUE)
  

  # loop across patients ----------------------------------------------------

  ids <- if (sample_pat == FALSE)  1:nrow(mx_b) else sort(sample(1:nrow(mx_b), sample_pat))

  # initiate parallel
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("cumhaz_exp", "cumhaz_wei", "cumhaz_gom", "cumhaz_cox", 
                      "p_event", "gen_next_cycle")) 
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(data.table))
  registerDoParallel(cl)
  
  retval <- foreach (i = ids) %dopar% {

    
    ### Monte Carlo simulation -----------------
    
    # initialise output dataframe
    df_output_i <- v_output_null
    
    ### events to model
    events_to_predict_i <- events_to_predict[i, -1]
    events_nf_i <- names(events_to_predict_i)[which(events_to_predict_i == 1)]
    
    events_nf_to_update <- c() 
    
    # initiate baseline values
    # independent for each simulation
    v_b <- mx_b[i, ]
    v_b_int <- v_b[int_names_b]
    v_b_int_nonlinage <- v_b[int_names_b_nonlinage]
    
    # set.seed(1234)
    #n <- 1
    # TODO: save default value of v_j, J & mx_output_null outside & then just reset it
    for (n in 1 : n_sim) {
      
      # baseline value of time-update characteristics
      v_j <- mx_t[i, ]
      
      # stopping value
      J <- max(eval(stop_expr) - 1, 1)
      
      # default output matrix 
      if (nonlinage) {
        mx_output_null <- matrix(nrow = J, ncol = N + 2) # append 2 age spline
        colnames(mx_output_null) <- c(out_names, "CurrAge_cent_1", "CurrAge_cent_2")
      } else {
        mx_output_null <- matrix(nrow = J, ncol = N)
        colnames(mx_output_null) <- out_names
      }
      
      # add QoL
      qol <- cbind("qol"=rep(NA, J))
      mx_output_null <- cbind(mx_output_null, qol)

      # add nsim
      nsim <- cbind("nsim"=rep(NA, J))
      mx_output_null <- cbind(mx_output_null, nsim)
      
      # RW: we must keep the sequence to add splines, qol and nsim, 
      # because in the cycle, splines are added first, and then qol.
      # nsim is added at the end of each cycle, outside the cycle function
      
      # define / reset baseline values
      # TODO: better 
      events_nf_to_predict <- events_nf_i
      alive <- 1
      had_rhabdo <- 0
      had_myop <- 0
      
      # lam_b (to be updated if off tx)
      lam_b_i <- list()
      
      intercept_calibrated_i <- list()
      e_list <- if (use_cost)
        c(events_all, "cost") else
          events_all
      
      # QoL
      e_list <- c(e_list, "qol")

      
      for (e in e_list) {
        lam_b_i[[e]] <- lam_b[[e]][i]
        # intercept_calibrated_i[[e]] <- intercept_calibrated[[e]][i] # RW
      }
      
      for (e in events_all){ # cost and qol do not have calibrated intercept # so does tx effect
        intercept_calibrated_i[[e]] <- intercept_calibrated[[e]][i] # RW 2021-03-17
      }
      
      # add individual tx effect # RW 2021-05-05
      tx_effect_i <- list()
      for (e in c("mi", "stroke", "crv", "vd")) {
        tx_effect_i[[e]] <- tx_effect[[e]][i]
      }
      
      # time since events
      # NA indicates no event
      # baseline events
      t_since_eb_i <- list()
      for (eb in eventsb)
        t_since_eb_i[[eb]] <- as.numeric(t_since_eb[[eb]][which(t_since_eb[[eb]][, "ids"] == i), "t"])
      # simulated events
      t_since_e_i <- list()
      for (e in events_nf_i) 
        t_since_e_i[[e]] <- NA
      
      ### loop across years --------------------
      
      # output for each simulation
      mx_output <- mx_output_null
      
      #j <- 1
      for (j in 1:J) {
        
        alpha <- gen_next_cycle(prim_flag = prim_flag, 
                                j = j, 
                                v_j = v_j, 
                                v_b_int = v_b_int, 
                                int_names_b = int_names_b,
                                nonlinage = nonlinage,
                                v_b_int_nonlinage = v_b_int_nonlinage, 
                                int_names_b_nonlinage = int_names_b_nonlinage,
                                p_crv = p_crv,
                                events_nf_to_predict = events_nf_to_predict, 
                                events_nf_to_update = events_nf_to_update, 
                                events_f = events_f,
                                eventsb = eventsb,
                                vars_t_eb = vars_t_eb, vars_t_e = vars_t_e,
                                t_since_eb_i = t_since_eb_i, t_since_e_i = t_since_e_i,
                                alive = alive,
                                dist = dist, 
                                lam_b_i = lam_b_i, 
                                cf_t = cf_t, 
                                shape = shape,
                                scaleby = scaleby, 
                                intercept_calibrated_i = intercept_calibrated_i,
                                use_cost = use_cost,
                                side_effects = side_effects,
                                had_rhabdo = had_rhabdo,
                                had_myop = had_myop,
                                statin_high = statin_high
        )
        
        # save recursive output
        v_j <- alpha$v_j
        events_nf_to_predict <- alpha$events_nf_to_predict
        # TODO: actually only need to update for events that happened 0-3 years ago
        events_nf_to_update <- setdiff(events_nf_i, events_nf_to_predict)
        alive <- alpha$alive
        had_rhabdo <- alpha$had_rhabdo
        had_myop <- alpha$had_myop
        if (had_rhabdo == 1 | had_myop == 1){
          # whichever turn 1, then turn 2; whichever stay 0, continue to stay 0
          had_rhabdo <- 2 * had_rhabdo
          had_myop <- 2 * had_myop
          # stop statin
          # update baseline lambdas
          for (e in c("mi", "stroke", "crv", "vd")) 
            lam_b_i[[e]] <- lam_b_i[[e]] / tx_effect_i[[e]]
          # turn off side effect to default values
          side_effects <- side_effect_default
        }
        
        # values of the time-updated covariates
        t_since_eb_i <- alpha$t_since_eb_i
        t_since_e_i <- alpha$t_since_e_i
        
        # record in the output matrix
        mx_output[j, ] <- c(v_j, n) # add sim number here
        
        # check whether the patient still alive
        if (alive == 0) {
          break
        }
      }

      ## add to df_output
      df_output_i <- rbind(df_output_i, mx_output)
    }
    
    # cleanup & return
    # TODO: this is still quite time-consuming
    # TODO: perhaps divide by n/sim outside of the cycle?
    if (save_by_pat)
      saveRDS(df_output_i, compress = F, # RW: use saveRDS
        file = file.path(output_dir, 
          str_c(output_filename_prefix, "_pat_", i, "_", prim_flag, ".rds")))
    

    # RW: because incident dm_0_10 will continue to be 1 across cycles once it changes from 0 to 1
    # sum(x)/n_sim will overestimate, because dm_0_10 accumulate over cycles
    # need correct: value in cycle i = value in cycle i - value in cycle (i-1) 
    df_output_i <- df_output_i %>% as.data.frame() %>% group_by(nsim) %>%
      mutate(dm_0_10 = pmax(dm_0_10 - lag(dm_0_10, default = 0), 0)) 
    # use pmax to avoid minus when jumping from dm_0_10 to dm_10_inf

    # summarise
    df_combined <- data.table(df_output_i)[!is.na(id_new),
                                           lapply(.SD, function(x) sum(x) / n_sim), 
                                           by = .(id_new, cycle),
                                           .SDcols = vars_to_summarise]
    
    return(df_combined)
  }
  stopCluster(cl)
  
  # combine and cleanup
  df_output <- do.call("rbind", retval) 
  
  # rename event columns
  # TODO; perhaps do it outside?
  for (e in events_nf) {
    varname_0_1 <- str_c(e, "_0_", vars_t_e[[e]][1])
    colnames(df_output)[which(colnames(df_output) == varname_0_1)] <- e
  }
  
  
  # save the output
  tx_tag <- if (!identical(tx, FALSE)) "tx" else "notx"
  
  regimen_tag <- str_replace_all(regimen, fixed(" "), "") # remove space
  
  side_effct_tag <- if (side_effects_flag) "sideon" else "sideoff"
  
  if (side_effects_flag & !side_dm) side_effct_tag <-"sideondmoff"
  
  length_tag <- if (is.expression(stop_expr)) "lifetime" else paste0(stop_expr,"yr")
  
  saveRDS(df_output, compress = F, # RW: use saveRDS
          file = file.path(output_dir, 
                           str_c(output_filename_prefix,"_",
                                 prim_flag, "_",
                                 tx_tag, "_",
                                 regimen_tag, "_", 
                                 side_effct_tag, "_",
                                 length_tag, "_",
                                 n_sim, "_",
                                 sens_multiplier,".rds")))
  return(df_output)
}
