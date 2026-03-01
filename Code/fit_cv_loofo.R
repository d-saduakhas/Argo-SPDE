# Fit-CV using the LOOFO (leave one floatID out)

# Load libraries
library(R.matlab)
# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(fmesher)
library(ngme2)
library(MASS)


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

fit_cv_loofo <- function(GridID, modelType, presLevel, N_gibbs, n_burnin, N_sim) {
  library(ngme2)
  base_dir <- "~/Results"
  save_directory <- file.path(base_dir, as.character(presLevel), modelType)
  save_filename <- file.path(save_directory, paste0(modelType, "_", GridID, "_new.RData"))
  load_cvfilename_new <- file.path(save_directory, paste0("cv_", modelType, "_", GridID, "_loofo.RData"))
  
  if (file.exists(load_cvfilename_new)) {
    return(NULL)
  }
  
  print(paste("Starting CV for iGrid:", GridID, "modelType:", modelType, "pressure:", presLevel))
  load(save_filename) # loads data, result_new, etc.
  
  # --- sanity check that floatID data aligns with fitted data ---
  floatID <- data$data
  floatID$group <- rep(c("Psal", "Temp"), each = length(floatID$Y) / 2)
  if (!is.data.frame(floatID) || !all(c("Y", "Lat", "Long", "group") %in% names(floatID))) {
    stop("get_floatID() did not return expected columns.")
  }
  
  if (!all(floatID$Y == data$data$Y) ||
      !all(floatID$Lat == data$data$Lat) ||
      !all(floatID$Long == data$data$Long) ||
      !all(floatID$group == data$group)) {
    # Instead of stopping let's extract the matching rows only and proceed. Use the fitted data as reference and match the floatID data to it.
    cat("Data/order mismatch detected between get_floatID() and fitted data. Aligning floatID data to fitted data...\n")
    matched_indices <- match(data$data$Y, floatID$Y)
    floatID <- floatID[matched_indices, ]
    cat("Alignment completed.\n")
  }
  # ---------------- end sanity check ----------------
  
  # original indexing helpers
  data$Lat <- data$data$Lat
  data$Long <- data$data$Long
  predLatMin <- data$predLat[1]
  predLatMax <- data$predLat[2]
  predLongMin <- data$predLong[1]
  predLongMax <- data$predLong[2]
  
  # Psal rows inside prediction box, in ORIGINAL index space (1..data$n)
  idx_psal_box <- which(
    data$Lat > predLatMin & data$Lat < predLatMax &
      data$Long > predLongMin & data$Long < predLongMax &
      data$group == "Psal"
  )
  
  # their paired Temp rows (n-offset)
  idx_temp_box <- idx_psal_box + data$n
  
  # unique FloatIDs present in the box (based on Psal rows)
  fids_in_box <- unique(floatID$FloatID[idx_psal_box])
  print(paste("Unique FloatIDs in prediction box:", length(fids_in_box)))
  
  
  # TEST: all rows (Psal+Temp) for one FloatID, but only those inside the box
  test_list <- lapply(fids_in_box, function(fid) {
    idx_psal_fid <- idx_psal_box[floatID$FloatID[idx_psal_box] == fid]
    c(idx_psal_fid, idx_psal_fid + data$n)
  })
  
  # TRAIN: complement within the box (original indexing)
  all_in_box <- c(idx_psal_box, idx_temp_box)
  train_list <- lapply(test_list, function(te) setdiff(all_in_box, te))
  
  start.time <- Sys.time()
  
  # USE UPDATED MODEL IN CROSS VALIDATION
  cv <- cross_validation(
    result_new,
    print = TRUE,
    parallel = T,
    type = "custom",
    test_idx = test_list,
    train_idx = train_list,
    thining_gap = 0,
    n_gibbs_samples = N_gibbs,
    cores_layer1 = 6,
    cores_layer2 = 1,
    n_burnin = n_burnin,
    N_sim = N_sim
  )
  
  end.time <- Sys.time()
  time.taken.cv <- round(end.time - start.time, 2)
  
  # fixed: save only the objects, not the path variable
  save(cv, time.taken.cv, file = load_cvfilename_new)
  return(cv)
}

grids <- c(1:386)
presLevel <- c(10, 1000)

for(p in presLevel){
  for(i in grids){
    tryCatch({
      fit_cv_loofo(i, "gauss_indep", p, 500, 0, 2)
      fit_cv_loofo(i, "gauss_cor", p, 500, 0, 2)
      fit_cv_loofo(i, "nig_indep", p, 1000, 200, 2)
      fit_cv_loofo(i, "nig_cor", p, 1000, 200, 2)
    }, error = function(e) {
      message(sprintf("Error in iGrid %d: %s", i, e$message))
    })
  }
}
