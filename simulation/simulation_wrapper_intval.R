###########################################################
### Internal validation, WITH CRV re-calibration & non-calibrated equations
###########################################################

### preamble ------------------------

rm(list = ls())

library(tidyverse)
library(parallel)
library(doParallel)
library(data.table)

### directories & filenames -------------------

# load all path names
source("Z:/PCTU/HEALTH ECONOMICS/CVD_HE/UKB/path.R")

# run function file
source(file.path(vali, "simulation_functions_v2.R"))

### parameters ---------------------

# Could split into 'common to all simulations' and 'specific to this simulation'
# For the moment copy-paste everything, for more transparency

# Also note that many parameters have default values and therefore there is no need to call those
# I will try to remember to provide an analogous example but with specifying only minimum number of parameters

### input folders and filenames

# use temp directory names not to overwrite global variables
.cf_dir <- vali
.input_data_dir <- vali_data

cf_filename <- "cf_cali_sp_recali" 
mx_b_filename_prefix <- "ukb_b"
mx_t_filename_prefix <- "ukb_t"

use_diab <- TRUE

### events to model

# TODO: this can be pre-saved
# decided against for now, because can change manually if needed
events_list <- list(
        events_nf = c ("mi", "stroke", "crv", "cancer_icd"), # RW: cancer_icd indicate incident cancer 
        events_f = c("vd", "nvd"))

if (use_diab)
        events_list[["events_nf"]] <- c(events_list[["events_nf"]], "dm")

events_to_predict_filename_prefix <- "events_to_predict"

# information on time-updated events

vars_t_filename <- "vars_t"
t_since_eb_filename_prefix <- "t_since_eb"

### calibration

# CRV re-calibration
# If FALSE, no re-calibration
# Otherwise a filename containing re-calibration proportions in a pre-defined format
.p_crv_dir <- vali_data
p_crv_filename <- "p_crv" # RW 

# use re-calibrated equation (TRUE/FALSE)
calibrated_eqns <- TRUE

# nonlinear age (TRUE/FALSE)
nonlinage <- TRUE

### treatment effects

# see examples below for non-zero treatment effects

tx <- FALSE

### distributions

# best AIC
dist_list <- list( # RW
        prim = list(
                mi = "gom", # previously "exp"
                stroke = "wei", # previously "gom"
                crv = "wei",
                cancer_icd = "gom", # previously "exp"
                vd = "exp",
                nvd = "wei"),
        sec = list(
                mi = "wei",
                stroke = "gom",
                crv = "wei",
                cancer_icd = "gom", # previously "exp"
                vd = "gom", # RW 2021-03-18 previously "exp"
                nvd = "wei")
)

# add distribution for diabetes if relevant
if (use_diab) {
        dist_list[["prim"]][["dm"]] <- "wei" # previously "gom"
        dist_list[["sec"]][["dm"]] <- "wei"
}
        

### simulation parameters

# Length of simulation
# stop_expr <- expr(12) # 5 years (including baseline year coded as 0) # RW: 9 years
# lifetime
stop_expr <- expression(floor(101 - (v_j["CurrAge_cent"]*10 + 60)))

# Should saving be done overall (FALSE) or on patient-level (TRUE)
# Note that saving on patient-level will be much more time-consuming
# Could be used to determine speed of convergence or debugging
save_by_pat <- F
# Should patients be sampled
# Could be FALSE (for complete sample) or an integer
sample_pat <- F

# use id list for conference abstract, primary prevention only 
id_list <- FALSE # RW 2021-03-19 

# number of simulations
n_sim <- 200

# number of cores
n_cores <- 24

### output directory

output_dir <- vali_output
# output_filename_prefix <- "ukb_vali_sim"

output_filename_prefix <- "ukb_qol_test" 

### simulation ---------------------

# prim_flag <- "prim"
# prim_flag <- "sec"
pop <- c("sec")

for (prim_flag in pop) { # RW
        
        # distributions
        dist <- dist_list[[prim_flag]]
        
        # output filename parameters
        
        ptm <- proc.time()
        
        # run the master function
        
        use_cost <- FALSE
        
        alpha <- master(
                .cf_dir = .cf_dir,
                .input_data_dir = .input_data_dir, 
                cf_filename = cf_filename, 
                prim_flag = prim_flag,
                mx_t_filename_prefix = mx_t_filename_prefix,
                mx_b_filename_prefix = mx_b_filename_prefix,
                events_list = events_list,
                events_to_predict_filename_prefix = events_to_predict_filename_prefix,
                vars_t_filename = vars_t_filename,
                t_since_eb_filename_prefix = t_since_eb_filename_prefix,
                .p_crv_dir = .p_crv_dir,
                p_crv_filename = p_crv_filename,
                calibrated_eqns = calibrated_eqns,
                nonlinage = nonlinage,
                use_diab = use_diab,
                tx = tx,
                dist = dist,
                stop_expr = stop_expr, 
                sample_pat = sample_pat,
                save_by_pat = save_by_pat,
                n_sim = n_sim, 
                n_cores = n_cores, 
                output_dir = output_dir, 
                output_filename_prefix = output_filename_prefix,
                id_list=id_list # RW: whether use pre-defined id list
        )
        
        print(proc.time() - ptm)
}