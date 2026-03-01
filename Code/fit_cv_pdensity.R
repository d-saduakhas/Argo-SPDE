# In this script we perform CV for potential density based on the fitted bivariate models
# Steps involved:
# 1. Create a function that uses the residual fits as input to compute potential density
# 2. Perform CV for potential density for 4 model types (nig_cor, nig_indep, gauss_cor, gauss_indep)
# 3. Compute bivariate log-likelihood, only available for gauss_cor models ####

## HELPER FUNCTIONS ####
library(ngme2)
library(fmesher)
library(inlabru)
library(Matrix)
library(MASS)
library(R.matlab)
library(gsw)
library(mvtnorm)


## HELPER FUNCTIONS####
# Extract parameter functions
extract_ngme_param_gauss <- function(result) {
  params <- list(
    Rho = result$replicates[[1]]$models[["spde"]]$theta_K[1],
    Sigma1 = exp(result$replicates[[1]]$models[["spde"]]$theta_K[2]),
    Sigma2 = exp(result$replicates[[1]]$models[["spde"]]$theta_K[3]),
    Theta = 0,
    Kappa1 = exp(result$replicates[[1]]$models[["spde"]]$theta_K[4]),
    Kappa2 = exp(result$replicates[[1]]$models[["spde"]]$theta_K[5]),
    sigma_e1 = exp(result$replicates[[1]]$noise$theta_sigma[1]),
    sigma_e2 = exp(result$replicates[[1]]$noise$theta_sigma[2]),
    rho_e = ifelse(
      !is.numeric(result$replicates[[1]]$noise$rho) ||
        length(result$replicates[[1]]$noise$rho) == 0,
      0,
      as.double(result$replicates[[1]]$noise$rho)
    )
  )
  return(params)
}

extract_ngme_param_nig <- function(result) {
  params <- list(
    Theta = 0,
    Rho = result$replicates[[1]]$models[["spde"]][["theta_K"]][1],
    Sigma1 = exp(result$replicates[[1]]$models[["spde"]][["theta_K"]][2]),
    Sigma2 = exp(result$replicates[[1]]$models[["spde"]][["theta_K"]][3]),
    Kappa1 = exp(result$replicates[[1]]$models[["spde"]][["theta_K"]][4]),
    Kappa2 = exp(result$replicates[[1]]$models[["spde"]][["theta_K"]][5]),
    Mu1 = result$replicates[[1]]$models[["spde"]]$noise$theta_mu[1],
    Mu2 = result$replicates[[1]]$models[["spde"]]$noise$theta_mu[2],
    Nu1 = exp(result$replicates[[1]]$models[["spde"]]$noise$theta_nu)[1],
    Nu2 = exp(result$replicates[[1]]$models[["spde"]]$noise$theta_nu)[2],
    sigma_e1 = exp(result$replicates[[1]][["noise"]][["theta_sigma"]][1]),
    sigma_e2 = exp(result$replicates[[1]][["noise"]][["theta_sigma"]][2]),
    rho_e = ifelse(
      !is.numeric(result$replicates[[1]]$noise$rho) ||
        length(result$replicates[[1]]$noise$rho) == 0,
      0,
      as.double(result$replicates[[1]]$noise$rho)
    )
  )
  return(params)
}

# Track problematic grids
problematic_grids <- list()

# 4. Fit cv Metric for potential density ####
fit_cv_metric <- function(GridID, modelType, presLevel, gibbs_samples = 500, burnin = 0) {
  library(ngme2)
  library(gsw)
  cat("POTENTIAL DENSITY CV - Grid:", GridID, "| Model:", modelType, "\n")
  
  # STEP 1: Load fitted model
  base_dir <- "~/Results"
  save_directory <- file.path(base_dir, as.character(presLevel), modelType)
  save_filename <- file.path(save_directory, paste0(modelType, "_", GridID, "_new.RData"))
  
  if (!file.exists(save_filename)) {
    cat("[SKIP] Model not found for Grid", GridID, "-", modelType, "\n")
    return(NULL)
  }
  # Ыkip if already done
  cv_save_filename <- file.path(save_directory, paste0("cv_", modelType, "_", GridID, "_pdensity.RData"))
  if (file.exists(cv_save_filename)) {
    cat("[SKIP] CV already done for Grid", GridID, "-", modelType, "\n")
    return(NULL)
  }
  
  cat("[1/5] Loading fitted model...\n")
  load(save_filename)
  cat("        Loaded\n\n")
  
  # STEP 2: Get test indices from model data (10×10 box)
  cat("[2/5] Identifying test locations in 10×10 box...\n")
  data$Lat <- data$data$Lat
  data$Long <- data$data$Long
  predLatMin <- data$predLat[1]
  predLatMax <- data$predLat[2]
  predLongMin <- data$predLong[1]
  predLongMax <- data$predLong[2]
  
  # Get indices of (Psal,Temp) observations in prediction box
  idx_year_cv <- which(data$Lat > predLatMin & data$Lat < predLatMax &
                         data$Long > predLongMin & data$Long < predLongMax)
  idx_year_cv <- idx_year_cv[idx_year_cv <= data$n] # Only Psal indices
  
  cat("      Test locations(per field):", length(idx_year_cv), "\n")
  cat("      Model total obs(per field):", data$n, "\n\n")
  
  # Create paired test/train lists
  test_list <- lapply(idx_year_cv, function(idx) c(idx, idx + data$n))
  train_list <- lapply(idx_year_cv, function(idx) {
    l1 <- setdiff(idx_year_cv, idx)
    l2 <- l1 + data$n
    c(l1, l2)
  })
  
  # STEP 3: Get RG means directly from data (no matching needed!)
  cat("[3/5] Extracting RG means from data...\n")
  
  rg_mean_psal <- data$data$estimatedMean[idx_year_cv]
  rg_mean_temp <- data$data$estimatedMean[idx_year_cv + data$n]
  long_vals <- data$data$Long[idx_year_cv]
  lat_vals <- data$data$Lat[idx_year_cv]
  
  cat("        Extracted", length(idx_year_cv), "RG means\n\n")
  
  # STEP 4: Create metric function
  cat("[4/5] Creating potential density metric...\n")
  metric <- local({
    rg_psal <- rg_mean_psal
    rg_temp <- rg_mean_temp
    longs <- long_vals
    lats <- lat_vals
    pres <- presLevel
    n_locs <- length(rg_psal)
    current_idx <- 0
    
    function(data_metric) {
      current_idx <<- current_idx + 1
      if (current_idx > n_locs) current_idx <<- 1
      i <- current_idx
      
      # Only process data_metric$y (named vector with Psal and Temp)
      res_psal <- data_metric$y["Psal"]
      res_temp <- data_metric$y["Temp"]
      
      # Add RG means
      actual_psal <- res_psal + rg_psal[i]
      actual_temp <- res_temp + rg_temp[i]
      
      # Compute potential density
      SA <- gsw_SA_from_SP(actual_psal, pres, longs[i], lats[i])
      pd <- gsw_pot_rho_t_exact(SA, actual_temp, pres, 0)
      
      # Return just the scalar
      list(
        y = as.numeric(pd),
        label = "pot_density"
      )
    }
  })
  cat("        Metric created\n")
  cat("      Transform: Residual + RGmean -> AbsSal  -> PotDensity\n\n")
  cat("      Gibbs samples:", gibbs_samples, "| N_sim: 2\n")
  start.time <- Sys.time()
  
  cv <- cross_validation(
    result_new,
    print = T,
    type = "custom",
    test_idx = test_list,
    train_idx = train_list,
    n_gibbs_samples = gibbs_samples,
    n_burnin = burnin,
    N_sim = 2,
    parallel = T,
    cores_layer1 = 5,
    cores_layer2 = 1,
    thining_gap = 0,
    metric = metric
  )
  
  time.taken <- round(Sys.time() - start.time, 2)
  
  cat("\n CV COMPLETED in", time.taken, "\n\n")
  cat("RESULTS:\n")
  print(cv$mean.scores)
  
  # Save
  save_filename <- file.path(save_directory, paste0("cv_", modelType, "_", GridID, "_pdensity.RData"))
  save(cv, time.taken, file = save_filename)
  
  return(cv)
}

# Example usage
for (i in 1:386) {
  fit_cv_metric(i, "nig_cor", 10, gibbs_samples = 1000, burnin = 200)
}

for (i in 1:386) {
  fit_cv_metric(i, "nig_indep", 10, gibbs_samples = 1000, burnin = 200)
}

for (i in 1:386) {
  fit_cv_metric(i, "gauss_indep", 10, gibbs_samples = 500, burnin = 0)
}

for (i in 1:386) {
  fit_cv_metric(i, "gauss_cor", 10, gibbs_samples = 500, burnin = 0)
}
