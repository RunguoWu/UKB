###########################################################
### External validation, With/out re-calibration 
###########################################################
library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)

rm(list = ls())

# batch control -----------------------------------------------------------

# simva, atorva, rosuva, prava, fluva
# 5mg, 10mg, 20mg, 40mg, 80mg
# ex. regimen could be "atorva 40mg", "none", or "atorva 40mg + ezetimibe"

# s1 <- c(tx="tx", regimen="none", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=1)
# s2 <- c(tx="tx", regimen="atorva 20mg", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=1)
# s3 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=1)
# s4 <- c(tx="tx", regimen="atorva 80mg", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=1)
# s5 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=FALSE, 
#         sens_multiplier=1)
# s6 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=FALSE, side_dm=FALSE, 
#         sens_multiplier=1)
# s7 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=3)
# s8 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE, 
#         sens_multiplier=30)

# scenario <- rbind(s1, s2, s3, s4, s5, s6, s7, s8)

s1 <- c(tx="tx", regimen="none", side_effects_flag=TRUE, side_dm=TRUE,
         sens_multiplier=1)

s2 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=FALSE, side_dm=FALSE,
        sens_multiplier=1)

s3 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=FALSE,
        sens_multiplier=1)

s4 <- c(tx="tx", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
        sens_multiplier=1)

scenario <- rbind(s1, s2, s3, s4)

# run loop ----------------------------------------------------------------


for (s in 1:nrow(scenario)) {
  
  # preparation -------------------------------------------------------------
  
  rm(list = setdiff(ls(), c("scenario","s")))

  # load all path names
  source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")
  
  # run function file
  source(file.path(vali, "functions_simulation.R"))
  
  
  # treatment effects control -----------------------------------------------
  
  # filename or FALSE
  # if tx != FALSE, baseline LDL will be replaced by pre-treated LDL in simulation
  tx <- scenario[s, "tx"]
  
  regimen <- scenario[s, "regimen"]
  
  # side effect flag: TRUE or FLASE
  # if regimen is none, side_effects_flag is forced to be FALSE in simulation
  side_effects_flag <- as.logical(scenario[s, "side_effects_flag"])
  
  side_dm <- as.logical(scenario[s, "side_dm"])
  
  sens_multiplier <- as.numeric(scenario[s, "sens_multiplier"])
  
  
  # i/o control -------------------------------------------------------------
  
  # input data directory
  .input_data_dir <- vali_data
  
  # output directory
  output_dir <- vali_output
  
  # output file name
  output_filename_prefix <- "sim" 
  
  # input coefficients file
  cf_filename <- "cf_20210508"
  
  # input baseline and time-varying data name prefixes
  mx_b_filename_prefix <- "ukb_b"
  mx_t_filename_prefix <- "ukb_t"
  
  # all events
  events_list <- list(
    events_nf = c ("mi", "stroke", "crv", "cancer_icd", "dm"),
    events_f = c("vd", "nvd"))
  
  # events_to_predict_filename_prefix <- "events_to_predict"
  # vars_t_filename <- "vars_t"
  # t_since_eb_filename_prefix <- "t_since_eb"
  # p_crv_filename <- "p_crv"
  
  # the 4 above files are combined into one pf file
  pf_filename <- "pf"
  
  # simulation parameter control --------------------------------------------
  
  # after new cf file came
  dist_list <- list( 
    prim = list(
      mi = "gom", 
      stroke = "wei", 
      crv = "wei",
      cancer_icd = "gom", 
      dm = "wei",
      vd = "gom",
      nvd = "wei"),
    sec = list(
      mi = "wei",
      stroke = "gom",
      crv = "gom", 
      cancer_icd = "gom", 
      dm = "wei",
      vd = "gom", 
      nvd = "gom") 
  )
  
  
  # nonlinear age (TRUE/FALSE)
  nonlinage <- TRUE # for QoL only
  
  # use re-calibrated equation (TRUE/FALSE)
  calibrated_eqns <- TRUE
  
  # just check if adjust crv-mi order
  adjust_crv = TRUE
  
  # Length of simulation
  # length = 12 # years
  # stop_expr <- expr(length + 1)
  
  # lifetime # maximum 101
  stop_expr <- expression(floor(101 - (v_j["CurrAge_cent"]*10 + 60)))
  
  # Should saving be done overall (FALSE) or on patient-level (TRUE)
  save_by_pat <- F
  # Should patients be sampled
  # Could be FALSE (for complete sample) or an integer
  sample_pat <- 10000
  
  # number of simulations
  n_sim <- 500
  
  # number of cores
  n_cores <- 24
  
  # populations to be simulated 
  pop <- c("prim")
  
  # prim_flag <- "sec"
  
  # simulation command ------------------------------------------------------
  
  for (prim_flag in pop) { # RW
    
    # distributions
    dist <- dist_list[[prim_flag]]
    
    # output filename parameters
    
    ptm <- proc.time()
    
    # run the master function
    
    use_cost <- FALSE
    
    alpha <- master(
      .input_data_dir = .input_data_dir, 
      cf_filename = cf_filename, 
      prim_flag = prim_flag,
      mx_t_filename_prefix = mx_t_filename_prefix,
      mx_b_filename_prefix = mx_b_filename_prefix,
      events_list = events_list,
      pf_filename = pf_filename,
      adjust_crv = adjust_crv,
      calibrated_eqns = calibrated_eqns,
      nonlinage = nonlinage,
      tx = tx,
      dist = dist,
      stop_expr = stop_expr, 
      sample_pat = sample_pat,
      save_by_pat = save_by_pat,
      n_sim = n_sim, 
      n_cores = n_cores, 
      output_dir = output_dir, 
      output_filename_prefix = output_filename_prefix, 
      regimen = regimen,
      side_effects_flag = side_effects_flag,
      side_dm = side_dm,
      sens_multiplier = sens_multiplier
    )
    
    print(proc.time() - ptm)
  }
}




