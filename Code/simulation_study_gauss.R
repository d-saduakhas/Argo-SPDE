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

# Function to simulate data
simulate_data <- function(param, loc_2d_mesh, mesh, n_obs, n_rep) {
    group <- c(rep("Temp", n_rep * n_obs), rep("Psal", n_rep * n_obs))
    repl <- rep(rep(1:n_rep, each = n_obs), 2)
    Long <- rep(loc_2d_mesh[, 1], 2 * n_rep)
    Lat <- rep(loc_2d_mesh[, 2], 2 * n_rep)
    true_model <- f(
        ~ Long + Lat,
        mesh = mesh,
        group = group,
        replicate = repl,
        model = "bv_matern_normal",
        sd1 = param$sigma[1],
        sd2 = param$sigma[2],
        rho = param$rho,
        theta = param$theta,
        sub_models = list(
            Temp = list(model = "matern", theta_K = log(param$kappa[1])),
            Psal = list(model = "matern", theta_K = log(param$kappa[2]))
        ),
        noise = list(Temp = noise_normal(sigma = 1), Psal = noise_normal(sigma = 1)),
        control = control_f(numer_grad = TRUE),
        name = "spde",
        eval = TRUE
    )
    W <- simulate(true_model)[[1]]
    # Note we have to generate the random noise independently for each replicate
    sd_1 <- param$sigma_eps[1]
    sd_2 <- param$sigma_eps[2]
    rho_e <- param$rho_eps
    Cov_same_idx <- matrix(c(sd_1^2, rho_e * sd_1 * sd_2, rho_e * sd_1 * sd_2, sd_2^2), nrow = 2)
    Cov_measurement <- Cov_same_idx %x% diag(n_obs * n_rep)
    L <- t(chol(Cov_measurement))
    e <- L %*% rnorm(2 * n_obs * n_rep)
    Y <- W + as.numeric(e)
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
    B_sigma[group == "Temp", 1] <- 1
    B_sigma[group == "Psal", 2] <- 1
    out_cor <- ngme(
        Y ~ 0 + f(
            ~ Long + Lat,
            mesh = mesh,
            model = "bv_matern_normal",
            name = "bv",
            sub_models = list(
                Temp = list(model = "matern"),
                Psal = list(model = "matern")
            ),
            noise = list(Temp = noise_normal(), Psal = noise_normal())
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
            stop_points = 50,
            optimizer = adam(stepsize = alpha),
            iterations = 5000,
            n_parallel_chain = 4,
            max_num_threads = 28,
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
            model = "bv_matern_normal",
            name = "bv",
            sub_models = list(
                Temp = list(model = "matern"),
                Psal = list(model = "matern")
            ),
            noise = list(Temp = noise_normal(), Psal = noise_normal())
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
            stop_points = 50,
            optimizer = adam(stepsize = alpha),
            iterations = 5000,
            n_parallel_chain = 4,
            max_num_threads = 28,
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
    folder_path <- "/ibex/scratch/saduakd/rho_simulation/"
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

            for (i in 1:10) {
                simulation_data <- simulate_data(param, loc_2d_mesh, mesh, n_obs, n_rep)
                save(simulation_data, mesh, loc_2d_mesh, param, n_obs, n_rep, file = paste0(folder_path, "simulation_data_", rho, "_rho_eps_", rho_eps, "_", i, ".RData"))
                # Estimate with the correlated model
                result_cor <- run_ngme_model(simulation_data$data, simulation_data$group, simulation_data$repl, mesh, 0.01, n_obs, n_rep)
                saveRDS(result_cor, file = paste0(folder_path, "result_rho_", rho, "_rho_eps_", rho_eps, "_cor_", i, ".rds"))
                # Estimate with the independent model
                result_ind <- run_ngme_model_ind(simulation_data$data, simulation_data$group, simulation_data$repl, mesh, 0.01, n_obs, n_rep)
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
# #rhos <- c( 0, 0.05, 0.2, 0.7,-0.7, -0.2, -0.05)
# rhos <- c( 2)
# rhos_eps <- c(-0.8, -0.4, -0.1, 0, 0.1, 0.4, 0.8)
# run_simulations(rhos, rhos_eps, n_obs = 1000, n_rep = 10)
