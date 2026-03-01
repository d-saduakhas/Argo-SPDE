library(ngme2)
library(fmesher)
library(inlabru)
library(ggplot2)
library(gridExtra)
library(viridis)

# Function to create mesh
create_mesh <- function(n_obs) {
    set.seed(42)
    long <- runif(n_obs, 0, 20)
    lat <- runif(n_obs, 0, 20)
    loc_2d_mesh <- cbind(long, lat)
    bnd <- spoly(data.frame(x = c(0, 20, 20, 0), y = c(0, 0, 20, 20)))
    mesh <- fm_mesh_2d_inla(
        boundary = bnd,
        loc.domain = loc_2d_mesh,
        max.edge = c(1, 5),
        cutoff = 0.1,
        offset = c(0.1, -0.1),
        min.angle = 20,
        max.n.strict = c(2500, 250)
    )
    return(list(mesh = mesh, loc_2d_mesh = loc_2d_mesh))
}

## OLD version 
# # Function to simulate data using version 0.6.0
# simulate_data <- function(param, loc_2d_mesh, mesh, n_obs, n_rep) {
#     group <- c(rep("Temp", n_rep * n_obs), rep("Psal", n_rep * n_obs))
#     repl <- rep(rep(1:n_rep, each = n_obs), 2)
#     Long <- rep(loc_2d_mesh[, 1], 2 * n_rep)
#     Lat <- rep(loc_2d_mesh[, 2], 2 * n_rep)
#     true_model <- f(
#         ~ Long + Lat,
#         mesh = mesh,
#         group = group,
#         replicate = repl,
#         model = "bv_matern_nig",
#         sd1 = param$sigma[1],
#         sd2 = param$sigma[2],
#         rho = param$rho,
#         theta = param$theta,
#         sub_models = list(
#             Temp = list(model = "matern", theta_K = log(param$kappa[1])),
#             Psal = list(model = "matern", theta_K = log(param$kappa[2]))
#         ),
#         noise = list(
#             Temp = noise_nig(sigma = 1, mu = param$mu[1], nu = param$nu[1]),
#             Psal = noise_nig(sigma = 1, mu = param$mu[2], nu = param$nu[2])
#         ),
#         control = control_f(numer_grad = TRUE),
#         name = "spde",
#         eval = TRUE
#     )
#     W <- simulate(true_model)[[1]]
#     # Note we have to generate the random noise independently for each replicate
#     sd_1 <- param$sigma_eps[1]
#     sd_2 <- param$sigma_eps[2]
#     rho_e <- param$rho_eps
#     Cov_same_idx <- matrix(c(sd_1^2, rho_e * sd_1 * sd_2, rho_e * sd_1 * sd_2, sd_2^2), nrow = 2)
#     Cov_measurement <- Cov_same_idx %x% diag(n_obs * n_rep)
#     L <- t(chol(Cov_measurement))
#     e <- L %*% rnorm(2 * n_obs * n_rep)
#     Y <- W + as.numeric(e)
#     data <- data.frame(
#         Long = rep(loc_2d_mesh[, 1], 2 * n_rep),
#         Lat = rep(loc_2d_mesh[, 2], 2 * n_rep),
#         Y = Y
#     )
#     return(list(data = data, group = group, repl = repl, Y = Y))
# }

# NEW Version for ngme2 0.7.1
simulate_data <- function(param, loc_2d_mesh, mesh, n_obs, n_rep) {
  # Simulate each replicate independently, then stack
  W_all <- c()
  
  # Per-replicate group (no replicates dimension)
  group_single <- c(rep("Temp", n_obs), rep("Psal", n_obs))
  Long_single <- rep(loc_2d_mesh[, 1], 2)
  Lat_single <- rep(loc_2d_mesh[, 2], 2)
  loc_single <- cbind(Long_single, Lat_single)
  
  for (r in 1:n_rep) {
    sim_model <- f(
      loc_single,
      mesh = mesh,
      group = group_single,
      model = "bv_matern_nig",  # or "bv_matern_normal" for Gaussian
      sd1 = param$sigma[1],
      sd2 = param$sigma[2],
      rho = param$rho,
      theta = param$theta,
      sub_models = list(
        Temp = list(model = "matern", theta_K = log(param$kappa[1])),
        Psal = list(model = "matern", theta_K = log(param$kappa[2]))
      ),
      noise = list(
        Temp = noise_nig(sigma = 1, mu = param$mu[1], nu = param$nu[1]),
        Psal = noise_nig(sigma = 1, mu = param$mu[2], nu = param$nu[2])
      ),
      control = control_f(numer_grad = TRUE),
      name = "spde",
      eval = TRUE
    )
    W_r <- simulate(sim_model)[[1]]
    W_all <- c(W_all, W_r)
  }
  
  # Now W_all has length 2 * n_obs * n_rep
  # Structure: [Temp_rep1, Psal_rep1, Temp_rep2, Psal_rep2, ...]
  # We need to reorganize to: [Temp_rep1, Temp_rep2, ..., Psal_rep1, Psal_rep2, ...]
  # to match the group/repl structure
  
  # Extract per-replicate blocks and reorganize
  W_temp <- c()
  W_psal <- c()
  for (r in 1:n_rep) {
    start_idx <- (r - 1) * 2 * n_obs + 1
    W_temp <- c(W_temp, W_all[start_idx:(start_idx + n_obs - 1)])
    W_psal <- c(W_psal, W_all[(start_idx + n_obs):(start_idx + 2 * n_obs - 1)])
  }
  W_stacked <- c(W_temp, W_psal)
  
  # Build group and replicate vectors
  group <- c(rep("Temp", n_rep * n_obs), rep("Psal", n_rep * n_obs))
  repl <- rep(rep(1:n_rep, each = n_obs), 2)
  
  # Correlated measurement noise
  sd_1 <- param$sigma_eps[1]
  sd_2 <- param$sigma_eps[2]
  rho_e <- param$rho_eps
  Cov_same_idx <- matrix(c(sd_1^2, rho_e * sd_1 * sd_2, 
                           rho_e * sd_1 * sd_2, sd_2^2), nrow = 2)
  
  # Generate noise per replicate (paired Temp/Psal at same location)
  e_all <- c()
  for (r in 1:n_rep) {
    Cov_r <- Cov_same_idx %x% diag(n_obs)
    L <- t(chol(Cov_r))
    e_r <- L %*% rnorm(2 * n_obs)  # [temp_noise, psal_noise]
    e_all <- c(e_all, as.numeric(e_r))
  }
  # e_all structure: [temp_r1, psal_r1, temp_r2, psal_r2, ...]
  # Reorganize to match group structure
  e_temp <- c()
  e_psal <- c()
  for (r in 1:n_rep) {
    start_idx <- (r - 1) * 2 * n_obs + 1
    e_temp <- c(e_temp, e_all[start_idx:(start_idx + n_obs - 1)])
    e_psal <- c(e_psal, e_all[(start_idx + n_obs):(start_idx + 2 * n_obs - 1)])
  }
  e_stacked <- c(e_temp, e_psal)
  
  Y <- W_stacked + e_stacked
  
  data <- data.frame(
    Long = rep(loc_2d_mesh[, 1], 2 * n_rep),
    Lat = rep(loc_2d_mesh[, 2], 2 * n_rep),
    Y = Y
  )
  return(list(data = data, group = group, repl = repl, Y = Y))
}

# Function to run ngme model
run_ngme_model <- function(data, group, repl, mesh, alpha, n_obs, n_rep) {
    B_sigma <- matrix(0, nrow = 2 * n_obs * n_rep, ncol = 2)
    B_sigma[group == "Psal", 1] <- 1
    B_sigma[group == "Temp", 2] <- 1
    out_cor <- ngme(
        Y ~ 0 + f(
            ~ Long + Lat,
            mesh = mesh,
            model = "bv_matern_nig",
            name = "bv",
            sub_models = list(
                Temp = list(model = "matern"),
                Psal = list(model = "matern")
            ),
            noise = list(Temp = noise_nig(), Psal = noise_nig())
        ),
        replicate = repl,
        data = data,
        group = group,
        family = noise_normal(
            corr_measurement = TRUE,
            index_corr = c(1:(n_obs * n_rep), 1:(n_obs * n_rep)),
            B_sigma = B_sigma,
            theta_sigma = c(0, 0)
        ),
        control_opt = control_opt(
            converge_eps = 1e-4,
            optimizer = precond_sgd(stepsize = alpha),
            stop_points = 50,
            iterations = 5000,
            n_parallel_chain = 2,
            max_num_threads = 20,
            rao_blackwellization = TRUE,
            print_check_info = FALSE
        ),
        debug = FALSE
    )
    return(out_cor)
}
# Function to run ngme model
run_ngme_model_ind <- function(data, group, repl, mesh, alpha, n_obs, n_rep) {
    B_sigma <- matrix(0, nrow = 2 * n_obs * n_rep, ncol = 2)
    B_sigma[group == "Temp", 1] <- 1
    B_sigma[group == "Psal", 2] <- 1
    out_cor <- ngme(
        Y ~ 0 + f(
            ~ Long + Lat,
            mesh = mesh,
            model = "bv_matern_nig",
            name = "bv",
            sub_models = list(
                Temp = list(model = "matern"),
                Psal = list(model = "matern")
            ),
            noise = list(Temp = noise_nig(), Psal = noise_nig())
        ),
        replicate = repl,
        data = data,
        group = group,
        family = noise_normal(
            corr_measurement = FALSE,
            B_sigma = B_sigma,
            theta_sigma = c(0, 0)
        ),
        control_opt = control_opt(
            converge_eps = 1e-4,
            optimizer = precond_sgd(stepsize = alpha),
            stop_points = 50,
            iterations = 5000,
            n_parallel_chain = 2,
            max_num_threads = 20,
            rao_blackwellization = TRUE,
            print_check_info = FALSE
        ),
        debug = FALSE
    )
    return(out_cor)
}



# Main function to run simulations and save results
run_simulations <- function(rhos, rhos_eps, n_obs, n_rep) {
    mesh_data <- create_mesh(n_obs)
    mesh <- mesh_data$mesh
    loc_2d_mesh <- mesh_data$loc_2d_mesh
    folder_path <- "/ibex/scratch/saduakd/rho_simulation_nig/"
    # inspired by parameters of grid 328
    for (rho in rhos) {
        for (rho_eps in rhos_eps) {
            param <- list(
                kappa = c(2, 1.5),
                sigma = c(0.5, 1),
                sigma_eps = c(0.5, 0.5),
                nu = c(0.3, 0.05),
                mu = c(1, 0.5),
                rho = rho,
                rho_eps = rho_eps,
                theta = 0
            )

            for (i in 1:2) {
                simulation_data <- simulate_data(param, loc_2d_mesh, mesh, n_obs, n_rep)
                save(simulation_data, mesh, loc_2d_mesh, param, n_obs, n_rep, file = paste0(folder_path, "simulation_data_", rho, "_rho_eps_", rho_eps, "_", i, ".RData"))
                # Estimate with the correlated model
                result_cor <- run_ngme_model(simulation_data$data, simulation_data$group, simulation_data$repl, mesh, 0.5, n_obs, n_rep)
                saveRDS(result_cor, file = paste0(folder_path, "result_rho_", rho, "_rho_eps_", rho_eps, "_cor_", i, ".rds"))
                # Estimate with the independent model
                result_ind <- run_ngme_model_ind(simulation_data$data, simulation_data$group, simulation_data$repl, mesh, 0.5, n_obs, n_rep)
                saveRDS(result_ind, file = paste0(folder_path, "result_rho_", rho, "_rho_eps_", rho_eps, "_ind_", i, ".rds"))
            }
        }
    }
}

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
rho <- as.numeric(args[1])
rho_eps <- as.numeric(args[2])
n_obs <- as.numeric(args[3])
n_rep <- as.numeric(args[4])

# Run the simulation with command-line arguments
run_simulations(rho, rho_eps, n_obs, n_rep)

# # Run simulations
# rhos <- c( 0.7, 0.05, 0.2,-0.7, -0.2, -0.05, 0)
# rhos_eps <- c(-0.8, -0.4, -0.1, 0, 0.1, 0.4, 0.8)
# run_simulations(rhos, rhos_eps, n_obs = 1000, n_rep = 1)
