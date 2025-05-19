# Compute the QQ plots for the standardized marginal residuals of the bivariate models
# *Only points that are within the grid are considered, no extensions, only the CV points are included!!
# Functions:
# * extract_ngme_param_gauss(result): Extract the estimated parameters from the NGME model, bv_matern_model
# * simulate_gauss_cor(n_obs_per_rep, n_replicate, ngme_out, mesh, loc, meas_noise): Simulate the bivariate Gaussian model, AW + meas_noise
# * simulate_nig_cor(n_obs_per_rep, n_replicate, ngme_out, mesh, loc, meas_noise): Simulate the bivariate NIG model
# * process_grid_model(grid_id, model_type): Compute the QQ plots for the given grid and model type
# * save.qq.plot(data, sim, name, main, lim): Save the QQ plot
# * save.cdf.plot(data, sim, name, main, lim): Save the CDF plot
# * save.pp.plot(data, sim, name, main, lim): Save the PP plot

rm(list = ls())
library(ngme2)
library(fmesher)
library(inlabru)
library(Matrix)
library(MASS)

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
        rho_e = ifelse(is.numeric(result[["replicates"]][[1]][["noise"]][["rho"]]), 0, result[["replicates"]][[1]][["noise"]][["rho"]])
    )

    required_parameters <- c("Rho", "Sigma1", "Sigma2", "Kappa1", "Kappa2", "sigma_e1", "sigma_e2", "rho_e")
    for (param in required_parameters) {
        if (!param %in% names(params) || is.null(params[[param]])) {
            params[[param]] <- NA # Assign NA if the parameter is missing or NULL
        }
    }

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
        rho_e = ifelse(is.numeric(result[["replicates"]][[1]][["noise"]][["rho"]]), 0, result[["replicates"]][[1]][["noise"]][["rho"]])
    ) # by default 0 if nig_indep

    required_parameters <- c("Theta", "Rho", "Sigma1", "Sigma2", "Kappa1", "Kappa2", "Mu1", "Mu2", "Nu1", "Nu2", "sigma_e1", "sigma_e2", "rho_e")
    for (param in required_parameters) {
        if (!param %in% names(params) || is.null(params[[param]])) {
            params[[param]] <- NA # Assign NA if the parameter is missing or NULL
        }
    }

    return(params)
}

simulate_gauss_cor <- function(n_obs_per_rep, n_replicate = 1, ngme_out = NULL, mesh = NULL, loc = NULL, meas_noise = FALSE, replicate = NULL) {
    if (!is.null(ngme_out)) {
        param <- extract_ngme_param_gauss(ngme_out)
        sigma <- c(param$Sigma1, param$Sigma2)
        rho <- param$Rho
        rho_e <- param$rho_e
        sigma_e <- c(param$sigma_e1, param$sigma_e2)
        kappa <- c(param$Kappa1, param$Kappa2)
    }
    group_per_rep <- c(rep("Psal", n_obs_per_rep / 2), rep("Temp", n_obs_per_rep / 2))
    idx_per_rep <- c(1:(n_obs_per_rep / 2), 1:(n_obs_per_rep / 2))

    repl <- replicate
    true_model <- f(loc,
        model = "bv_matern_normal", rho = rho, sd1 = sigma[1], sd2 = sigma[2],
        sub_models = list(
            Psal = list(model = "matern", theta_K = log(kappa[1])),
            Temp = list(model = "matern", theta_K = log(kappa[2]))
        ),
        mesh = mesh, group = group_per_rep, replicate = repl,
        noise = list(Psal = noise_normal(sigma = 1), Temp = noise_normal(sigma = 1)),
        control = control_f(numer_grad = TRUE)
    )
    W <- simulate(true_model)[[1]]
    Cov_same_idx <- matrix(c(
        sigma_e[1]^2, sigma_e[1] * sigma_e[2] * rho_e,
        sigma_e[1] * sigma_e[2] * rho_e, sigma_e[2]^2
    ), 2)
    Cov_measurement <- Cov_same_idx %x% diag(n_obs_per_rep * n_replicate / 2)
    L <- t(chol(Cov_measurement))
    e <- L %*% rnorm(n_obs_per_rep * n_replicate)
    Y <- W + as.numeric(e)

    formula <- Y ~ 0 + f(~ Long + Lat,
        model = "bv_matern_normal", rho = rho, sd1 = sigma[1], sd2 = sigma[2],
        sub_models = list(
            Psal = list(model = "matern", theta_K = log(kappa[1])),
            Temp = list(model = "matern", theta_K = log(kappa[2]))
        ),
        mesh = mesh, group = group_per_rep, replicate = repl,
        noise = list(Psal = noise_normal(sigma = 1), Temp = noise_normal(sigma = 1)),
        control = control_f(numer_grad = TRUE), name = "spde"
    )

    group <- rep(group_per_rep, n_replicate)
    # repl <- rep(1:n_replicate, each = n_obs_per_rep)

    return(list(
        Y = Y, loc = loc, formula = formula, true_model = true_model,
        group_per_rep = group_per_rep, rho = rho, group = group, repl = repl, mesh = mesh))
}

simulate_nig_cor <- function(n_obs_per_rep, n_replicate = 1, ngme_out = NULL, mesh = NULL, loc = NULL, meas_noise = FALSE, replicate = NULL) {
    if (!is.null(ngme_out)) {
        param <- extract_ngme_param_nig(ngme_out)
        sigma <- c(param$Sigma1, param$Sigma2)
        rho <- param$Rho
        rho_e <- param$rho_e
        sigma_e <- c(param$sigma_e1, param$sigma_e2)
        kappa <- c(param$Kappa1, param$Kappa2)
        nu <- c(param$Nu1, param$Nu2)
        mu <- c(param$Mu1, param$Mu2)
        theta <- 0
    }
    group_per_rep <- c(rep("Psal", n_obs_per_rep / 2), rep("Temp", n_obs_per_rep / 2))
    idx_per_rep <- c(1:(n_obs_per_rep / 2), 1:(n_obs_per_rep / 2))
    repl <- replicate
    true_model <- f(loc,
        model = "bv_matern_nig", rho = rho, theta = 0, sd1 = sigma[1], sd2 = sigma[2],
        sub_models = list(
            Psal = list(model = "matern", theta_K = log(kappa[1])),
            Temp = list(model = "matern", theta_K = log(kappa[2]))
        ),
        mesh = mesh, group = group_per_rep, replicate = repl,
        noise = list(
            Psal = noise_nig(sigma = 1, mu = mu[1], nu = nu[1]),
            Temp = noise_nig(sigma = 1, mu = mu[2], nu = nu[2])
        ),
        control = control_f(numer_grad = TRUE)
    )

    W <- simulate(true_model)[[1]]
    Cov_same_idx <- matrix(c(
        sigma_e[1]^2, sigma_e[1] * sigma_e[2] * rho_e,
        sigma_e[1] * sigma_e[2] * rho_e, sigma_e[2]^2
    ), 2)
    Cov_measurement <- Cov_same_idx %x% diag(n_obs_per_rep * n_replicate / 2)
    L <- t(chol(Cov_measurement))
    e <- L %*% rnorm(n_obs_per_rep * n_replicate)
    Y <- W + as.numeric(e)

    formula <- Y ~ 0 + f(~ Long + Lat,
        model = "bv_matern_nig", rho = rho, theta = theta, sd1 = sigma[1], sd2 = sigma[2],
        sub_models = list(
            Psal = list(model = "matern", theta_K = log(kappa[1])),
            Temp = list(model = "matern", theta_K = log(kappa[2]))
        ),
        mesh = mesh, group = group_per_rep, replicate = repl,
        noise = list(
            Psal = noise_nig(sigma = 1, mu = mu[1], nu = nu[1]),
            Temp = noise_nig(sigma = 1, mu = mu[2], nu = nu[2])
        ),
        control = control_f(numer_grad = TRUE), name = "spde"
    )

    group <- rep(group_per_rep, n_replicate)

    return(list(
        Y = Y, loc = loc, formula = formula, true_model = true_model,
        group_per_rep = group_per_rep, rho = rho, group = group, repl = repl, mesh = mesh))
}

confidence.bands <- function(ngme_out, loc, mesh, n.iter, nSampleSize, gauss = FALSE, meas_noise = TRUE, repl = NULL) {
    Y.res <- list()
    U.res <- list()
    for (k in 1:n.iter) {
        # 1. simulate one field realization
        cat(k)
        if (gauss) {
            if (meas_noise) { # correlated meas_noise
                sim <- simulate_gauss_cor(n_obs_per_rep = dim(loc)[1], ngme_out = ngme_out, loc = loc, mesh = mesh, meas_noise = TRUE, repl = repl)
            } else { # uncorrelated meas_noise
                sim <- simulate_gauss_cor(n_obs_per_rep = dim(loc)[1], ngme_out = ngme_out, loc = loc, mesh = mesh, meas_noise = FALSE, repl = repl)
            }
        } else { # for NIG
            # the replicate structure is not capture correctly
            if (meas_noise) {
                sim <- simulate_nig_cor(n_obs_per_rep = dim(loc)[1], ngme_out = ngme_out, loc = loc, mesh = mesh, meas_noise = TRUE, repl = repl)
            } else {
                sim <- simulate_nig_cor(n_obs_per_rep = dim(loc)[1], ngme_out = ngme_out, loc = loc, mesh = mesh, meas_noise = FALSE, repl = repl)
            }
        }
        # To replace the past data with new simulated data, do estimation FALSE
        if (gauss) {
            param <- extract_ngme_param_gauss(ngme_out)
        } else {
            param <- extract_ngme_param_nig(ngme_out)
            nu <- c(param$Nu1, param$Nu2)
            mu <- c(param$Mu1, param$Mu2)
            theta <- param$Theta
        }
        sigma <- c(param$Sigma1, param$Sigma2)
        rho <- param$Rho
        rho_e <- param$rho_e
        sigma_e <- c(param$sigma_e1, param$sigma_e2) # Order of the fields - Temp Psal
        kappa <- c(param$Kappa1, param$Kappa2)
        repl_full <- sim$repl
        group_per_rep <- sim$group_per_rep
        formula <- sim$formula
        mesh <- sim$mesh
        Long <- sim$loc[, 1]
        Lat <- sim$loc[, 2]
        n <- dim(sim$loc)[1] / 2 # Each field has n total observations

        result <- ngme(formula,
            replicate = repl_full, group = group_per_rep,
            data = data.frame(Y = sim$Y, Long = loc[, 1], Lat = loc[, 2]),
            control_opt = control_opt(estimation = FALSE),
            family = noise_normal(
                theta_sigma = log(sigma_e),
                rho = rho_e,
                corr_measurement = meas_noise,
                index_corr = rep(1:n, 2), # both index_corr & B_sigma are structured for case when two fields are stacked: Y = c(field1, field2)
                B_sigma = do.call(rbind, replicate(1, matrix(c(rep(1, n), rep(0, n), rep(0, n), rep(1, n)), 2 * n, 2), simplify = FALSE))
            )
        )

        # 2.Predict the simulated model
        pred.data <- predict(result,
            map = list(spde = loc), group = group_per_rep,
            sampling_size = nSampleSize
        )

        # Goal: simulated_data - mean predicted data for the process only without the measurement noise involoved
        Y.res[[k]] <- unlist(sim$Y)
        U.res[[k]] <- pred.data$mean
    }
    return(list(
        Y = Y.res,
        U = U.res
    ))
}

### PLOT FUNCTIONS ###
create_qq_plot <- function(res_data, sim_data, plot_name, plot_title) {
    library(ggplot2)
    # Create a data frame with data and simulation quantiles
    qq_data <- data.frame(
        DataQuantiles = sort(res_data),
        SimulationQuantiles = sort(sim_data[[1]])
    )

    for (i in 2:length(sim_data)) {
        new_data <- data.frame(
            DataQuantiles = sort(res_data),
            SimulationQuantiles = sort(sim_data[[i]])
        )
        qq_data <- rbind(qq_data, new_data)
    }

    p <- ggplot(qq_data, aes(x = DataQuantiles, y = SimulationQuantiles)) +
        geom_point(alpha = 0.5, color = "grey") +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        labs(title = plot_title, x = "Data Quantiles", y = "Simulation Quantiles") +
        theme_minimal()
    ggsave(plot_name, plot = p, width = 8, height = 8, units = "in")
}

save.qq.plot <- function(data, sim, name, main, lim = c(0.001, 0.999)) {
    # pdf(name, width = 5, height = 5 , pointsize = 13, family = "Helvetica")
    png(name, width = 1000, height = 1000, pointsize = 13, family = "Helvetica")
    xlim <- quantile(data, probs = lim, na.rm = TRUE)
    ylim <- quantile(unlist(sim[[1]]), probs = lim)
    ylim[1] <- min(ylim[1], xlim[1])
    ylim[2] <- max(ylim[2], xlim[2])
    plot(sort(data), sort(sim[[1]]),
        type = "l", col = "grey", xlim = xlim, ylim = ylim,
        xlab = "Data quantiles", ylab = "Simulation quantiles", main = main
    )
    for (i in 2:length(sim)) {
        lines(sort(data), sort(sim[[i]]), col = "grey")
    }
    lines(sort(data), sort(data), col = 1, type = "l")
    dev.off()
}

save.cdf.plot <- function(data, sim, name, main, lim = c(0.001, 0.999)) {
    # pdf(name, width = 5, height = 5 , pointsize = 13, family = "Helvetica")
    png(name, width = 1000, height = 1000, pointsize = 13, family = "Helvetica")
    xlim <- quantile(data, probs = lim)
    plot(ecdf(sim[[1]]), col = "grey", xlim = xlim, main = main)
    for (i in 2:length(sim)) {
        plot(ecdf(sim[[i]]), col = "grey", add = TRUE)
    }
    plot(ecdf(data), col = 1, add = TRUE, do.points = FALSE)
    dev.off()
}
save.pp.plot <- function(data, sim, name, main, lim = c(0.001, 0.999)) {
    # pdf(name, width = 5, height = 5 , pointsize = 13, family = "Helvetica")
    png(name, width = 1000, height = 1000, pointsize = 13, family = "Helvetica")
    q <- quantile(data, probs = c(0.01, 0.99))
    x <- seq(from = q[1], to = q[2], length.out = 1000)
    c1 <- ecdf(data)
    c2 <- ecdf(sim[[1]])

    plot(c1(x), c2(x),
        col = "grey", xlim = c(0, 1), main = main, type = "l",
        xlab = "Data percentiles", ylab = "Simulation percentiles"
    )
    for (i in 2:length(sim)) {
        c2 <- ecdf(sim[[i]])
        lines(c1(x), c2(x), col = "grey")
    }
    lines(c1(x), c1(x), col = 1)
    dev.off()
}

presLevel <- 300

model_types <- c("gauss_indep", "gauss_cor", "nig_indep", "nig_cor")

setwd("Results/300/QQ_plots/")

process_grid_model <- function(grid_id, model_type) {
    file_path <- paste0("Results/300/", model_type, "/", model_type, "_", grid_id, "_new.RData")
    # Check if already done
    if (file.exists(paste0("Grid", grid_id, "_", model_type, "_results.RData"))) {
        cat("QQ for Grid", grid_id, "and model", model_type, "already done Skipping...\n")
        return()
    } else {
        cat("Processing grid", grid_id, "and model", model_type)
    }
    # cat("Processing grid", grid_id,"and model", model_type)
    # Check if the file exists
    if (!file.exists(file_path)) {
        cat("File for Grid", grid_id, "and model", model_type, "does not exist. Skipping...\n")
        return()
    }
    load(file_path)
    for (i in 1:result_new$n_repls) {
        result_new$replicates[[i]]$models[[1]]$operator$bv_theta <- result_new$replicates[[i]]$models[[1]]$operator$theta_K[1]
    }

    data$Lat <- data$data$Lat
    predLatMin <- data$predLat[1]
    predLatMax <- data$predLat[2]
    data$Long <- data$data$Long
    predLongMin <- data$predLong[1]
    predLongMax <- data$predLong[2]

    # Subset the data points within the prediction Grid only, Only CV-ed data points will be used!!!
    idx_year_cv <- which(data$Lat > predLatMin & data$Lat < predLatMax & data$Long > predLongMin & data$Long < predLongMax & data$group == "Psal")
    n_pred <- length(idx_year_cv) # number of points for 1 field to be predicted for all replicates
    coord.prd <- cbind(data$Long[idx_year_cv], data$Lat[idx_year_cv]) # locations of points for 1 field
    coord.prd = rbind(coord.prd,coord.prd) # for two fields
    n <- n_pred
    obs_data <- data.frame(
        Salinity = data$data$Y[idx_year_cv],
        Temperature = data$data$Y[idx_year_cv + data$n]
    )

    pred.data <- predict(result_new, map = list(spde = rbind(coord.prd, coord.prd)), group = c(rep("Psal", n), rep("Temp", n)), sampling_size = 1000) # making the predictions

    # simulated_data - mean predicted data for the process only without the measurement noise involved
    bb <- pred.data$mean # opposite way given Salinity Temp
    res <- c(obs_data$Salinity, obs_data$Temperature) - bb # salinity,temperature
    mesh <- data$mesh
    loc <- coord.prd # for two fields
    repl = c(data$repl[idx_year_cv],data$repl[idx_year_cv+data$n])d

    if (model_type == "nig_cor") {
        cb <- confidence.bands(ngme_out = result_new, loc = loc, mesh = data$mesh, n.iter = 20, nSampleSize = 1000, repl = repl)
    } else if (model_type == "nig_indep") {
        cb <- confidence.bands(ngme_out = result_new, loc = loc, mesh = data$mesh, n.iter = 20, nSampleSize = 1000, meas_noise = FALSE, repl = repl)
    } else if (model_type == "gauss_cor") {
        cb <- confidence.bands(ngme_out = result_new, loc = loc, mesh = data$mesh, n.iter = 20, gauss = TRUE, nSampleSize = 1000, repl = repl)
    } else if (model_type == "gauss_indep") {
        cb <- confidence.bands(ngme_out = result_new, loc = loc, mesh = data$mesh, n.iter = 20, gauss = TRUE, nSampleSize = 1000, meas_noise = FALSE, repl = repl)
    }
    # order of Y salinity, temperature
    Y <- (lapply(1:length(cb$Y), function(i) cb$Y[[i]] - cb$U[[i]]))
    n <- dim(loc)[1] / 2
    print(n)
    res_psal <- res[1:n]
    Y_psal <- (lapply(1:length(Y), function(i) Y[[i]][1:n]))
    res_temp <- res[(n + 1):(2 * n)]
    Y_temp <- (lapply(1:length(Y), function(i) Y[[i]][(n + 1):(2 * n)]))

    results <- list(data_res = res, sim_res = Y, bb = bb, cb = cb, obs_data = obs_data)
    save(results, file = paste0("Grid", grid_id, "_", model_type, "_results.RData"))

    plot_name_prefix <- paste0("Grid", grid_id, "_", model_type, "_temp")
    save.qq.plot(res_temp, Y_temp, paste0(plot_name_prefix, "QQ.png"), paste0(model_type, " Temp"))
    # save.pp.plot(res_temp, Y_temp, paste0(plot_name_prefix, "PP.png"), paste0(model_type, " Temp"))
    # save.cdf.plot(res_temp, Y_temp, paste0(plot_name_prefix, "CDF.png"), paste0(model_type, " Temp"))
    # create_qq_plot(res_temp, Y_temp, paste0(plot_name_prefix, "QQggplot.png"), "Gauss Cor Temp QQ Plot")

    plot_name_prefix <- paste0("Grid", grid_id, "_", model_type, "_psal")
    save.qq.plot(res_psal, Y_psal, paste0(plot_name_prefix, "QQ.png"), paste0(model_type, " Psal"))
    # save.pp.plot(res_psal, Y_psal, paste0(plot_name_prefix, "PP.png"), paste0(model_type, " Psal"))
    # save.cdf.plot(res_psal, Y_psal, paste0(plot_name_prefix, "CDF.png"), paste0(model_type, " Psal"))
    # create_qq_plot(res_psal, Y_psal,  paste0(plot_name_prefix, "QQggplot.png"), "Gauss Cor Psal QQ Plot")
}

problematic_grids <- list() 

for (i in 6:390) {
    tryCatch(
        {
            process_grid_model(i, "gauss_cor")
            process_grid_model(i, "gauss_indep")
            process_grid_model(i, "nig_cor")
            process_grid_model(i, "nig_indep")
        },
        error = function(e) {
            # Record the problematic grid and error message
            problematic_grids[[length(problematic_grids) + 1]] <- list(grid = i, error = e$message)
            message(sprintf("Error in grid %d: %s", i, e$message))
        }
    )
}

if (length(problematic_grids) > 0) {
    print("Problematic grids encountered:")
    print(problematic_grids)
} else {
    print("No errors encountered.")
}
