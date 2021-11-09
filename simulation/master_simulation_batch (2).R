###########################################################
### External validation, With/out re-calibration 
###########################################################

rm(list = ls())

library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)

# batch control -----------------------------------------------------------

# simva, atorva, rosuva, prava, fluva
# 5mg, 10mg, 20mg, 40mg, 80mg
# ex. regimen could be "atorva 40mg", "none", or "atorva 40mg + ezetimibe"

# s0: the same as validation
s0 <- c(tx=FALSE, regimen="none", side_effects_flag=FALSE, side_dm=FALSE,
        side_mus=FALSE, sens_multiplier=1, pop=c("sec"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=250, id_list=FALSE)

s1 <- c(tx="tx", regimen="none", side_effects_flag=FALSE, side_dm=FALSE,
        side_mus=FALSE, sens_multiplier=1, pop=c("sec"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

s2 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=FALSE, side_dm=FALSE,
        side_mus=FALSE, sens_multiplier=1, pop=c("sec"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

s3 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
         side_mus=FALSE, sens_multiplier=1, pop=c("sec"), sample_nodm = FALSE,
         dm_predict = 1, n_sim=500, id_list=TRUE)

s4 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
        side_mus=TRUE, sens_multiplier=1, pop=c("sec"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)


s5 <- c(tx="tx", regimen="none", side_effects_flag=FALSE, side_dm=FALSE,
        side_mus=FALSE, sens_multiplier=1, pop=c("prim"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

s6 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=FALSE, side_dm=FALSE,
        side_mus=FALSE, sens_multiplier=1, pop=c("prim"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

s7 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
        side_mus=FALSE, sens_multiplier=1, pop=c("prim"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

s8 <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
        side_mus=TRUE, sens_multiplier=1, pop=c("prim"), sample_nodm = FALSE,
        dm_predict = 1, n_sim=500, id_list=TRUE)

# s1a <- c(tx="tx", regimen="none", side_effects_flag=FALSE, side_dm=FALSE,
#         side_mus=FALSE, sens_multiplier=2, pop=c("sec"), sample_nodm = FALSE,
#         dm_predict = 1, n_sim=500, id_list=TRUE)
# 
# s4a <- c(tx="tx88", regimen="atorva 40mg", side_effects_flag=TRUE, side_dm=TRUE,
#         side_mus=FALSE, sens_multiplier=2, pop=c("sec"), sample_nodm = FALSE,
#         dm_predict = 1, n_sim=500, id_list=TRUE)


# scenario <- rbind(s7, s8)


# st <- c(tx=FALSE, regimen="none", side_effects_flag=FALSE, side_dm=FALSE,
#         side_mus=FALSE, sens_multiplier=1, pop=c("prim"), sample_nodm = FALSE,
#         dm_predict = 1, n_sim=100)

# scenario <- rbind(s2, s3, s5, s6, s7)
scenario <- rbind(s7)

# run loop ----------------------------------------------------------------


for (s in 1:nrow(scenario)) {
  
  # preparation -------------------------------------------------------------
  
  rm(list = setdiff(ls(), c("scenario","s")))

  # load all path names
  source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")
  
  # run function file
  source(file.path(simu, "functions_simulation.R"))
  
  
  # treatment effects control -----------------------------------------------
  
  # filename or FALSE
  # if tx != FALSE, baseline LDL will be replaced by pre-treated LDL in simulation
  tx <- scenario[s, "tx"]
  
  if (tx == "FALSE") tx <- as.logical(tx)
  
  regimen <- scenario[s, "regimen"]
  
  # side effect flag: TRUE or FALSE
  # if regimen is none, side_effects_flag is forced to be FALSE in simulation
  side_effects_flag <- as.logical(scenario[s, "side_effects_flag"])
  
  side_dm <- as.logical(scenario[s, "side_dm"])
  
  side_mus <- as.logical(scenario[s, "side_mus"])
  
  sens_multiplier <- as.numeric(scenario[s, "sens_multiplier"])
  
  
  # i/o control -------------------------------------------------------------
  
  # input data directory
  .input_data_dir <- vali_data
  
  # output directory
  output_dir <- simu_output
  
  # output file name
  output_filename_prefix <- "s211106_vdnvdWEI" 
  
  # input coefficients file
  cf_filename <- "cf_20210714_cali1022_final"
  
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
      mi = "wei", 
      stroke = "wei", 
      crv = "wei",
      cancer_icd = "exp", 
      dm = "wei",  # gom previous
      vd = "wei",
      nvd = "wei"),
    sec = list(
      mi = "exp",
      stroke = "exp",
      crv = "gom", 
      cancer_icd = "exp", 
      dm = "wei",
      vd = "wei", 
      nvd = "wei")  
  )
  
  
  # nonlinear age (TRUE/FALSE)
  nonlinage <- TRUE # for QoL only
  
  # use re-calibrated equation (TRUE/FALSE)
  calibrated_eqns <- TRUE
  
  # just check if adjust crv-mi order
  adjust_crv = TRUE
  
  # Length of simulation
  # 12 years
  # length = 12
  # stop_expr <- expr(length + 1)
  
  # lifetime # maximum 101
  stop_expr <- expression(floor(101 - (v_j["CurrAge_cent"]*10 + 60)))
  
  # Should saving be done overall (FALSE) or on patient-level (TRUE)
  save_by_pat <- FALSE
  # Should patients be sampled
  # Could be FALSE (for complete sample) or an integer
  sample_pat <- FALSE
  
  # use pre-selected groups of individuals?
  id_list <- as.logical(scenario[s, "id_list"])
  
  # sample 800 people without baseline dm 
  sample_nodm <- as.logical(scenario[s, "sample_nodm"])  
  
  dm_predict <- as.numeric(scenario[s, "dm_predict"]) 
  
  # number of simulations
  n_sim <- as.numeric(scenario[s, "n_sim"])
  
  # number of cores
  n_cores <- round(detectCores() * 2/3)
  
  # populations to be simulated 
  prim_flag <- scenario[s, "pop"]
  
  # simulation command ------------------------------------------------------
  
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
    side_mus = side_mus,
    sens_multiplier = sens_multiplier,
    sample_nodm = sample_nodm,
    dm_predict = dm_predict,
    id_list = id_list
  )
  
  print(proc.time() - ptm)
  print(Sys.time())
}




