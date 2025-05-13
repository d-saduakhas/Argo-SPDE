load_and_prepare_data <- function(iGrid, presLevel) {
    # Load libraries
    library(R.matlab)
    library(fmesher)
    library(ngme2)
    library(MASS)
    library(inlabru)

    # Load data
    S <- readMat(paste0("Data/residual_", as.character(presLevel), "_01.mat"))
    interpLatYear <- S$interpLat
    interpLongYear <- S$interpLong
    interpJulDayYear <- S$interpJulDay
    interpTempYear <- S$differenceTemp
    interpPsalYear <- S$differencePsal
    # general conditions for the data
    if (presLevel == 10) {
        drop_id <- which(interpTempYear < -15 | interpTempYear > 10 | interpPsalYear < -4 | interpPsalYear > 4)
        interpTempYear <- interpTempYear[-drop_id]
        interpLatYear <- interpLatYear[-drop_id]
        interpLongYear <- interpLongYear[-drop_id]
        interpPsalYear <- interpPsalYear[-drop_id]
        interpJulDayYear <- interpJulDayYear[-drop_id]
    }
    if (presLevel == 300) {
        drop_id <- which(interpPsalYear < -5)
        interpTempYear <- interpTempYear[-drop_id]
        interpLatYear <- interpLatYear[-drop_id]
        interpLongYear <- interpLongYear[-drop_id]
        interpPsalYear <- interpPsalYear[-drop_id]
        interpJulDayYear <- interpJulDayYear[-drop_id]
    }
    if (presLevel == 1000) {
        drop_id <- which(interpPsalYear < -10 | interpTempYear < -10)
        interpTempYear <- interpTempYear[-drop_id]
        interpLatYear <- interpLatYear[-drop_id]
        interpLongYear <- interpLongYear[-drop_id]
        interpPsalYear <- interpPsalYear[-drop_id]
        interpJulDayYear <- interpJulDayYear[-drop_id]
    }
    # Data preparation
    month <- 1
    startYear <- 2007
    endYear <- 2020
    nYear <- endYear - startYear + 1
    n_rep <- nYear
    grid_table_ext <- readMat("Data/grid_equal.mat")
    Grid <- grid_table_ext$Grid
    # Prediction Point of approximately 10x10 degrees, shifted by 20 degrees to match RG masked data
    predLatMin <- Grid[iGrid, 1]
    predLatMax <- Grid[iGrid, 2]
    predLongMin <- Grid[iGrid, 3] + 20
    predLongMax <- Grid[iGrid, 4] + 20

    # The box of coordinates used for parameter estimation
    windowSize <- 5
    latMin <- predLatMin - windowSize
    latMax <- predLatMax + windowSize
    longMin <- predLongMin - windowSize
    longMax <- predLongMax + windowSize

    interpLatAggr <- vector("list", nYear)
    interpLongAggr <- vector("list", nYear)
    interpTempAggr <- vector("list", nYear)
    interpPsalAggr <- vector("list", nYear)
    idx_cv <- vector("list", nYear)
    # Calculate the number of leap years between 0000 and 1970
    leap_years <- sum((0:1969) %% 4 == 0) - sum((0:1969) %% 100 == 0) + sum((0:1969) %% 400 == 0)

    # Calculate the total number of days between MATLAB's and R's reference dates
    offset <- 1970 * 365 + leap_years
    for (iYear in startYear:endYear) {
        startJulDay <- as.numeric(as.Date(paste0(iYear, "-01-01"), format = "%Y-%m-%d")) + offset
        endJulDay <- as.numeric(as.Date(paste0(iYear, "-12-31"), format = "%Y-%m-%d")) + offset

        # Filter based on interpJulDay
        idx <- which(interpLatYear > latMin & interpLatYear < latMax & interpLongYear > longMin & interpLongYear < longMax & interpJulDayYear >= startJulDay & interpJulDayYear <= endJulDay)

        interpLatAggr[[iYear - startYear + 1]] <- interpLatYear[idx]
        interpLongAggr[[iYear - startYear + 1]] <- interpLongYear[idx]
        interpTempAggr[[iYear - startYear + 1]] <- interpTempYear[idx]
        interpPsalAggr[[iYear - startYear + 1]] <- interpPsalYear[idx]

        # idx for data inside the prediction box of 10 by 10 degrees - for cv later
        idx_cv[[iYear - startYear + 1]] <- which(interpLatYear > predLatMin & interpLatYear < predLatMax & interpLongYear > predLongMin & interpLongYear < predLongMax & interpJulDayYear >= startJulDay & interpJulDayYear <= endJulDay)
    }
    n_pred <- sum(sapply(idx_cv, length))
    n <- sum(sapply(interpTempAggr, length))
    if (n_pred < 100) {
        return(list(error = "Not Enough Data")) # not enough data
    }

    years_vector <- unlist(sapply(startYear:endYear, function(y) rep(y, length(interpLatAggr[[y - 2006]]))))

    data <- data.frame(
        Y = c(unlist(interpPsalAggr), unlist(interpTempAggr)), # ordered alphabetically
        Long = rep(unlist(interpLongAggr), 2),
        Lat = rep(unlist(interpLatAggr), 2),
        Year = rep(years_vector, 2) - 2006
    )
    # Generate mesh with fmesher
    max_edge <- 1
    bnd <- spoly(data.frame(x = c(longMin, longMax, longMax, longMin), y = c(latMin, latMin, latMax, latMax)))
    mesh <- fm_mesh_2d_inla(
        boundary = bnd,
        loc.domain = cbind(data$Long, data$Lat),
        max.edge = c(1, 5) * max_edge,
        cutoff = 0.1,
        offset = c(0.1, -0.1),
        min.angle = 20,
        max.n.strict = c(3000, 300) # maximum precision
    )
    return(list(
        data = data,
        mesh = mesh,
        predLong = c(predLongMin, predLongMax),
        predLat = c(predLatMin, predLatMax),
        group = c(rep("Psal", n), rep("Temp", n)),
        repl = c(years_vector, years_vector),
        n = n
    ))
}

# Fit NIG independent model, i.e. independent measurement error
fit_nig_indep <- function(data, mesh, group, repl, n, n_opt, result_gauss = NULL) {
    result <- ngme(
        Y ~ 0 + f(~ Long + Lat,
            mesh = mesh,
            model = "bv_matern_nig", theta = 0, fix_bv_theta = T,
            sub_models = list(Psal = "matern", Temp = "matern"),
            noise = list(Psal = noise_nig(), Temp = noise_nig()),
            control = control_f(numer_grad = T),
            name = "spde",
            eval = T
        ),
        replicate = repl,
        data = data,
        group = group,
        control_ngme = control_ngme(
            n_gibbs_samples = n_opt[4]
        ),
        family = noise_normal(
            theta_sigma = c(0, 0),
            corr_measurement = F,
            index_corr = rep(1:n, 2), # both index_corr & B_sigma are structured for case when two fields are stacked: Y = c(field1, field2)
            B_sigma = do.call(rbind, replicate(1, matrix(c(rep(1, n), rep(0, n), rep(0, n), rep(1, n)), 2 * n, 2), simplify = FALSE))
        ),
        control_opt = control_opt(
            converge_eps = 1e-05,
            optimizer = precond_sgd(stepsize = 0.3),
            burnin = n_opt[1],
            estimation = T,
            iterations = n_opt[2],
            n_parallel_chain = n_opt[3],
            verbose = F,
            print_check_info = T,
            max_num_threads = 28,
            stop_points = 200,
            std_lim = 0.1,
            trend_lim = 0.01,
            max_relative_step = 0.25,
            max_absolute_step = 0.25,
            n_slope_check = 4,
            rao_blackwellization = TRUE
        ),
        start = result_gauss
    )

    return(result)
}

# Fit NIG correlated model, i.e. correlated measurement error
fit_nig_cor <- function(data, mesh, group, repl, n, n_opt, result_gauss = NULL) {
    result <- ngme(
        Y ~ 0 + f(~ Long + Lat,
            mesh = mesh,
            model = "bv_matern_nig", theta = 0, fix_bv_theta = T,
            sub_models = list(Psal = "matern", Temp = "matern"),
            noise = list(Psal = noise_nig(), Temp = noise_nig()),
            control = control_f(numer_grad = T),
            name = "spde",
            eval = T
        ),
        replicate = repl,
        data = data,
        group = group,
        control_ngme = control_ngme(
            n_gibbs_samples = n_opt[4]
        ),
        family = noise_normal(
            theta_sigma = c(0, 0),
            corr_measurement = T,
            index_corr = rep(1:n, 2), # both index_corr & B_sigma are structured for case when two fields are stacked: Y = c(field1, field2)
            B_sigma = do.call(rbind, replicate(1, matrix(c(rep(1, n), rep(0, n), rep(0, n), rep(1, n)), 2 * n, 2), simplify = FALSE))
        ),
        control_opt = control_opt(
            converge_eps = 1e-05,
            optimizer = precond_sgd(stepsize = 0.3),
            burnin = n_opt[1],
            estimation = T,
            iterations = n_opt[2],
            n_parallel_chain = n_opt[3],
            verbose = F,
            print_check_info = T,
            max_num_threads = 28,
            stop_points = 200,
            std_lim = 0.1,
            trend_lim = 0.01,
            max_relative_step = 0.25,
            max_absolute_step = 0.25,
            n_slope_check = 4,
            rao_blackwellization = TRUE
        ),
        start = result_gauss
    )
    return(result)
}

fit_model_new <- function(iGrid, modelType, presLevel) {
    print(paste("Running with iGrid:", iGrid, "modelType:", modelType, "pressure:", presLevel))

    # Check if the model fit exists, if not, run the model
    base_dir <- "Results"
    save_directory <- file.path(base_dir, as.character(presLevel), modelType)
    save_filename <- file.path(save_directory, paste0(modelType, "_", iGrid, "_new.RData"))
    # UNCOMMENT to enable check of exisiting model fits
    if (file.exists(save_filename)) {
        print(paste("Model already fit for iGrid:", iGrid, "modelType:", modelType, "pressure:", presLevel))
        return(NULL)
    }

    if (!file.exists(save_filename)) {
        # n_opt structure: 1.burn-in samples 2.total iterations 3.number of chains 4. Gibbs samples per iteration
        n_opt <- c(100, 10000, 2, 3)
        start.time <- Sys.time()
        print(paste("New fit for iGrid:", iGrid, "modelType:", modelType, "pressure:", presLevel))
        data <- load_and_prepare_data(iGrid, presLevel)
        # Check if there's not enough data
        if (is.list(data) && !is.null(data$error) && data$error == "Not Enough Data") {
            print(paste("Not enough data for GridID:", iGrid))
            return(NULL)
        }
        # load the corresponding Gaussian model fit as the starting point
        if (modelType == "nig_indep") {
            file_name <- file.path(base_dir, as.character(presLevel), "gauss_indep", paste0("gauss_indep", "_", iGrid, "_new.RData"))
            if (file.exists(file_name)) {
                env <- new.env()
                load(file_name, env)
                result_gauss <- env$result_new
                print("Gauss model loaded as pre-estimate.")
            } else {
                cat("The file", file_name, "does not exist.\n")
                result_gauss <- NULL
            }
            result_new <- fit_nig_indep(data$data, data$mesh, data$group, data$repl, data$n, n_opt, result_gauss)
        } else if (modelType == "nig_cor") {
            file_name <- file.path(base_dir, as.character(presLevel), "gauss_cor", paste0("gauss_cor", "_", iGrid, "_new.RData"))
            if (file.exists(file_name)) {
                env <- new.env()
                load(file_name, env)
                result_gauss <- env$result_new
                print("Gauss model loaded as pre-estimate.")
            } else {
                cat("The file", file_name, "does not exist.\n")
                result_gauss <- NULL
            }
            result_new <- fit_nig_cor(data$data, data$mesh, data$group, data$repl, data$n, n_opt, result_gauss)
        }
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        save(result_new, data, time.taken, file = save_filename)
    }
}

fit_cv <- function(GridID, modelType, presLevel) {
    library(ngme2)
    # Perform cross-validation
    base_dir <- "Results"
    save_directory <- file.path(base_dir, as.character(presLevel), modelType)
    save_filename <- file.path(save_directory, paste0(modelType, "_", GridID, "_new.RData"))

    load_cvfilename_new <- file.path(save_directory, paste0("cv_", modelType, "_", GridID, ".RData"))
    if (file.exists(load_cvfilename_new)) {
        print(paste("Model CV already done for iGrid:", GridID, "modelType:", modelType, "pressure:", presLevel))
        return(NULL)
    }
    print(paste("Starting CV for iGrid:", GridID, "modelType:", modelType, "pressure:", presLevel))
    load(save_filename)
    data$Lat <- data$data$Lat
    predLatMin <- data$predLat[1]
    predLatMax <- data$predLat[2]
    data$Long <- data$data$Long
    predLongMin <- data$predLong[1]
    predLongMax <- data$predLong[2]

    idx_year_cv <- which(data$Lat > predLatMin & data$Lat < predLatMax & data$Long > predLongMin & data$Long < predLongMax & data$group == "Psal")

    test_list <- lapply(idx_year_cv, function(idx) {
        return(c(idx, idx + data$n))
    })
    train_list <- lapply(idx_year_cv, function(idx) {
        l1 <- setdiff(idx_year_cv, idx)
        l2 <- l1 + data$n
        return(c(l1, l2))
    })
    start.time <- Sys.time()

    cv <- cross_validation(result_new,
        print = T,
        type = "custom", test_idx = test_list, train_idx = train_list,
        n_gibbs_samples = 200,
        n_burnin = 200,
        N = 2, parallel = T,
        cores_layer1 = 7,
        cores_layer2 = 1
    )

    end.time <- Sys.time()
    time.taken.cv <- round(end.time - start.time, 2)
    time.taken.cv

    save_filename <- file.path(save_directory, paste0("cv", modelType, "_", GridID, "_new.RData"))
    save(result_new, cv, data, time.taken.cv, file = save_filename)
    return(cv)
}

# When this script is run, it will take in command line arguments
if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    GridID <- as.numeric(args[1])
    modelType <- args[2]
    presLevel <- as.numeric(args[3])
    # Execute the model
    fit_model_new(GridID, modelType, presLevel)
    fit_cv(GridID, modelType, presLevel)
}
