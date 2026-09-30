# Set of functions used to create / handle / analyse models
##### Libraries ##### ---------------------------------------------------------
library(sf)
library(cli)
library(coda)
library(Hmsc)
library(tidyr)
library(dplyr)
library(abind)
library(ggplot2)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters


##### Data ##### --------------------------------------------------------------
# A function to turn a list of lists into a dataframe with combinations as rows.
# ARGS:
#   - combinations: a list of lists.
# RETURNS:
#   - a table where each row is a different combination of the given lists.
build_param_grid <- function(combinations) {
    # create a grid, each line is a combination of variables
    # that way, we loop on this table instead of having nested loops
    tables_list <- vector("list", length(combinations))

    for (i in seq_along(combinations)) {
        tables_list[[i]] <- expand.grid(
            combinations[[i]], stringsAsFactors = FALSE)
    }

    table <- do.call(rbind, tables_list)
    table <- unique(table)   # in case blocks overlap, e.g. the (none, 0, none) row
    rownames(table) <- NULL

    return(table |> distinct())
}

# A function to turn the row of a dataframe into a folder path for a model run.
# ARGS:
#   - folder: a string, the path to a parent folder.
#   - params: a row of a dataframe.
#   - k_fold: the k_fold ID.
# RETURNS:
#   - a folder path.
make_run_path <- function(folder, params, k_fold){
    row_string <- tolower(
        paste(names(params), unlist(params), sep = "-", collapse = "_")
    )
    
    return(file.path(folder, paste0("model_", row_string, "_", k)))
}


# A function to automatise verbose interpretation of diagnostic vectors.
# ARGS:
#   - vector: a vector of values to analyse.
#   - bad: a numeric. The threshold under which values reveal a bad fit.
#   - good: a numeric. The threshold over which values reveal a good fit.
#   - order: a string. Indicates if lower is better ("low_better", default) 
#       or higer is bette (high_better).
#   - mode: a sting. Indicates if showing full analysis ("full", default) 
#       or a quick 1-line summary ("quick").
# RETURNS:
#   - NULL, prints cli messages in console
interpret_diagnostics <- function(
    vector, 
    good, 
    bad = NULL, 
    mode = "full",
    order = "low_better") {
    
    if (order == "low_better") {
        n_good <- sum(vector < good)    

        if (is.null(bad)) {
            n_bad <- sum(vector > good)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector > good) & (vector < bad))
            n_bad <- sum(vector > bad)    
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: <", good," (good), ",
                    good, "-", bad, ", (acceptable), >",
                    bad, " (bad)"))
            }
        }
        
    } else if (order == "high_better") {
        n_good <- sum(vector > good)    

        if (is.null(bad)) {
            n_bad <- sum(vector < good)  
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good)"))
            }

            n_acceptable <- 0
        } else {
            n_acceptable <- sum((vector < good) & (vector > bad))
            n_bad <- sum(vector < bad)
            if (mode == "full") {
                cli_alert_info(paste0(
                    "Rule of thumb: >", good," (good), ",
                    good, "-", bad, ", (acceptable), <",
                    bad, " (bad)"))   
            }
        }
     
    }
    else {
        stop(paste0("Mode should be one of 'low_better' or 'high_better'. ",
        "Got '", mode, "' ."))
    }
    
    if (mode == "full") {
        if (n_good > 0) {
            cli_alert_success(paste0(
            "- Number of 'good' estimates: ", n_good, 
            " (", round(100*n_good/length(vector), 2), "% of given values)"
            ))
        }
        if ((n_acceptable > 0) & !(is.null(bad))) {
            cli_alert_info(paste0(
            "- Number of 'acceptable' estimates: ", n_acceptable, 
            " (", round(100*n_acceptable/length(vector), 2), "% of given values)"
            ))  
        }
        if (n_bad > 0) {
            cli_alert_info(paste0(
            "- Number of 'bad' estimates: ", n_bad, 
            " (", round(100*n_bad/length(vector), 2), "% of given values)"
            ))
        }
    } else if (mode == "quick") {
        if (is.null(bad)) {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        } else {
            cli_alert_info(paste0(
                n_good, " (", round(100*n_good/length(vector), 2), "%) are good, ",
                n_acceptable, " (", round(100*n_acceptable/length(vector), 2), "%) are acceptable ",
                "and ", n_bad ," (", round(100*n_bad/length(vector), 2),"%) are bad."))
        }
        
    }

}

##### HMSC ##### --------------------------------------------------------------
# A function to prepare dataset for Hmsc training, outputs a model template
# ARGS:
#   - subdataset: a dataframe. Must contain columns listed in x_cols and y_cols.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - y_cols: a list of strings. The columns containing species occurrences.
#   - formula: a formula for the Hmsc model. Based on names in x_cols/y_cols.
#   - random_effect: whether to add 
#       a random effect to Hmsc ("points", "carres" or "spatial") 
#       or not ("none", default).
# RETURNS:
#   - the initialised HMSC object.
prepare_hmsc_training <- function(
        subdataset, 
        x_cols, y_cols, 
        formula,
        id_column = "row_id",
        random_effect = "none"){
    # Format y and x data to a format that Hmsc accepts.
    ydata <- as.matrix(subdataset[y_cols]) 
    xdata <- as.data.frame(setNames(
        lapply(x_cols, function(col) subdataset[[col]]),
        x_cols)
    )
    studyDesign <- data.frame(
        spatial = as.factor(subdataset[[id_column]])
    )

    # Depending on the random effect, the model has to be set up differently
    # FYI: here, we use occurence data. For that, a "probit" distribution is
    # the most logical option (but others are available such as "normal",
    # "poisson", "lognormal poisson")
    if (random_effect == "none") {
        # No effect: easy, do not add anything to the model.
        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            XFormula = formula, 
            distr = "probit",           
            studyDesign = studyDesign)
        
    }  else if (random_effect == "spatial") {
        # Random effect as "spatial": we consider the coordinates of each point
        # sampled. The random effect is a function of the distance between the
        # points. Takes more time to compute but (usually), yields to better 
        # results in prediction.
        # convert coordinates to metric
        coords_sf <- st_as_sf(
            subdataset, 
            coords = c("lon", "lat"), 
            crs = 4326)
        coords_proj <- st_transform(coords_sf, crs = 2154)

        # associate coordinate to the name of each sampled point
        xy <- st_coordinates(coords_proj)
        rownames(xy) <- as.character(subdataset[[id_column]])
        colnames(xy) <- c("longitude_grid_2154", "latitude_grid_2154")

        # There should be one observation per point (because i_point_annee)
        # but we must make sure or else Hmsc will fail.
        if (any(duplicated(xy))) {
            stop("Duplicate coordinates found across distinct points; check LON/LAT data.")
            # xy <- xy[!duplicated(rownames(xy)), ] # quick fix for this error
        }
        
        # There are two possibilities : map the entire grid of points or use
        # nearest neighbour approximation. After ~1000 points its better to use
        # the approximation for faster computation.
        if (nrow(xy) < 1000) {
            rL.spatial <- HmscRandomLevel(sData = xy)
        } else {
            rL.spatial <- HmscRandomLevel(
                sData = xy, sMethod = "NNGP", nNeighbours = 10)
        }
        rL.spatial = setPriors(
            rL.spatial, nfMin =1,  nfMax = 1) # limit number of latent variables

        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            ranLevels = list("spatial" = rL.spatial),
            XFormula = formula, 
            distr = "probit",
            studyDesign = studyDesign)

    } else {
        stop(paste0(
            "Random effect must be one of 'none', 'spatial'. Got ", 
            random_effect))
    }

    cli_alert_success("Created Hmsc object!\n\n")
    return(hmsc_object)
}

# A function to fit a Hmsc model and save results.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - save_to: a string. The path where the file will be saved (end with .rds).
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - nsamples: a numeric. The number of samples to collect
#   - ntransient: a numeric. The number of steps to wait for before 
#       collecting samples.
#   - freq_verbose: a numeric (default is 100). The frequency of verbose 
#       messages, you get one every "freq_verbose" step.
#   - allow_parallel: a boolean (default is TRUE). 
#       If nchains >1, allows to compute chains in parallel for faster fitting. 
#       Removes verbose messages during fitting.
# RETURNS:
#   - the fitted HMSC object.
fitting_hmsc <- function(
    hM, 
    nchains,
    thin,
    nsamples,
    ntransient,
    freq_verbose = 100,
    allow_parallel = TRUE,
    save_to = NULL) {
    # Start clock to get runtime
    start <- Sys.time()
    cli_alert_info(paste0("Started fitting at: ", start))

    # When HMSC runs chains in parallel, it does not show progress bar
    # Since I was expecting some update messages during the runs at first, 
    # its probably better to print some warning messages.
    if (!is.null(freq_verbose) & (nchains > 1) & allow_parallel) {
        cli_alert_warning(paste0("Cannot display fitting progress when running",
        " chains in parallel."))
        cli_alert_warning("Set `allow_parallel` to `FALSE` to see progress.")
    }

    # Depending on allow_parallel = TRUE/FALSE, overwrite nParallel
    fitted.hmsc <-  sampleMcmc(
        hM, 
        thin = thin, 
        samples = nsamples, 
        transient = ntransient, 
        nChains = nchains, 
        nParallel = if (allow_parallel) nchains else 1, 
        updater=list(GammaEta=FALSE),
        verbose = freq_verbose)
    
    # Print runtime
    stop <- Sys.time()
    cli_alert_info(paste0("Completed fitting at: ", stop))
    cli_alert_info(paste0("Time elapsed: ", round(stop-start, 2)))

    # Save run 
    if (!is.null(save_to)) {
        cli_alert_info("Saving model...")
        saveRDS(fitted.hmsc, file = save_to)
        cli_alert_success("Model saved!\n\n")    
    }
    return(fitted.hmsc)
}

# A function that automates the process of making predictions with Hmsc models.
# ARGS:
#   - hM: a Hmsc fitted model object.
#   - df: a dataframe with columns "point", id_column
#       and names in x_variables.
#   - x_variables: a list of strings. The names of columns to keep in data.
# RETURNS:
#   - see Hmsc::predict() returned object
predict_hmsc <- function(hM, df, x_variables, id_column = "row_id") {
    # Format explanatory variables to Hmsc expected format
    XData <- as.data.frame(setNames(
        lapply(x_variables, function(col) df[[col]]),
        x_variables))

    if ("spatial" %in% names(hM$ranLevels)) {
        # Format coordinates associated to each point
        coords_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
        coords_proj <- st_transform(coords_sf, crs = 2154)
        xy_new <- st_coordinates(coords_proj)
        rownames(xy_new) <- as.character(df[[id_column]])
        colnames(xy_new) <- c("longitude_grid_2154", "latitude_grid_2154")

        # Use Gradient (instead of studyDesign), for spatial gradients
        Gradient <- prepareGradient(
            hM,
            XDataNew = XData,
            sDataNew = list(spatial = xy_new)
        )

        # Make prediction on new dataset
        return(predict(hM, Gradient = Gradient, expected = TRUE))

    } else {
        # Format study design (grouping of samples together)
        studyDesign <- data.frame(
            spatial = as.factor(df[[id_column]])
        )

        # Make prediction on new dataset
        return(predict(
            hM, XData = XData, studyDesign = studyDesign, expected = TRUE))
    }
}

# A function that computes the uncertainty of a Hmsc model on its predictions.
# Uncertainty can be defined in many different ways. Here, uncertainty refers
# to the average of the standard deviation accross observed species
# ARGS:
#   - hM: a Hmsc fitted model object.
#   - df: a dataframe with columns "point", id_column
#       and names in x_variables.
#   - x_variables: a list of strings. The names of columns to keep in data.
# RETURNS:
#   - The measured uncertainty
get_uncertainty_hmsc <- function(hM, df, x_cols, id_column = "row_id") {
    # Get predictions
    predicted_occurrences <- predict_hmsc(
        hM = hM, df = df, x_variables = x_cols, id_column = id_column)
    
    # predicted_occurrences contains a list of length = number of samples.
    # for each sample, we get a matrix of n_obs x n_species
    sd_point_sp <- apply(simplify2array(predicted_occurrences), c(1,2), sd)
    uncertainty_per_point <- as_tibble(sd_point_sp) |> 
        mutate(average_sd = rowMeans(across(everything())))
    
    return(uncertainty_per_point$average_sd)
}

# A function to display convergence diagnostics for a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - nchains: a numeric. The number of chains to run.
#   - thin: a numeric. The number of steps between each recording of a sample.
#   - save_folder: a string. The path where plotted PDFs will be saved.
# RETURNS:
#   - NULL: saves objects along the way
convergence_hmsc <- function(hM, nchains, thin, save_folder) {
    # Hmsc object can be (and should be) converted to the commonly used Coda
    # format in Bayesian statistics.
    coda_outputs <- convertToCodaObject(hM)
    
    # Summary plots (not necessary here but might be useful one day):
    # MCMCsummary(object = coda_outputs$Beta, round = 2) 
    # MCMCplot(object = coda_outputs$Beta)

    # HMSC has many parameters, in our case we are only interested in:
    #   - Beta: fixed effects
    #   - Omega: random variation in co-occurence
    for (param in c("Beta", "Omega")) {
        # skip if name not in parameters (e.g. no omega if no random effect)
        if (!(param %in% names(coda_outputs))) {
            next
        } 
        cli_alert_info(paste0("*** [Parameter: ", param, "] ***"))

        ### Convergence diagnostics
        # Convert Omega output to match Beta format
        if (param == "Omega") {
            chains <- coda_outputs$Omega[[1]]  # 3D array: [iter, sp, sp]
        } else {
            chains <- coda_outputs[[param]]
        }
        
        ## Traceplot (Rhat and effective size)      
        cli_alert_info("Computation of traceplots, can take some time...")
        # Fast version:
        tryCatch({ # avoids debugger mode
            MCMCtrace(
                object = chains,
                pdf = TRUE,
                filename = file.path(
                    save_folder, 
                    paste0(param, "_all_traceplots.pdf")),
                ind = TRUE,
                open_pdf = FALSE,
                plot = TRUE,
                Rhat = TRUE, 
                n.eff = TRUE, 
                type = "both" # explicitly request trace + density
            )}, 
        error = function(e) {
            message("Error in MCMCtrace: ", e$message)
        })
        # # Slower but "more beautiful" version
        # traceplots <- ggplot_custom_MCMCtrace(
        #     coda_object = coda_fitted_model$Beta,
        #     show_Rhat = TRUE,
        #     show_Neff = TRUE)       
        # for (i_plot in seq_along(traceplots)) {
        #     standardised_ggplot_save(
        #         figure = (traceplots[[i_plot]]$trace + traceplots[[i_plot]]$density), 
        #         save_path = file.path(
        #             save_folder, 
        #             paste0("Beta_", i_plot, "_traceplot.pdf")))
        # }
        cli_alert_info("Traceplots saved!")


        ## Effective size
        cli_alert_info("-> Effective size:")
        eff_size <- effectiveSize(chains)
        eff_size_plot <- ggplot_bars(
                as.data.frame(eff_size), "eff_size", 
                breaks = ceiling(c(0, seq(100, max(eff_size), length.out = 19))),
                underlayers = list(
                    bad = annotate("rect", 
                        xmin = -Inf, xmax = 100, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 100, xmax = 400, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 400, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
        print(eff_size_plot)
        interpret_diagnostics(
            eff_size, bad = 100, good = 400, order = "high_better", mode = "quick")
        standardised_ggplot_save(
            figure = eff_size_plot, 
            save_path = file.path(
                save_folder, 
                paste0("hist_neff_", param, ".pdf")))

        ## Gelman-Rubin convergence diagnostic
        if (nchains > 1) {
            cli_alert_info("-> Gelman-Rubin convergence diagnostic:")
            psrf <- gelman.diag(
                chains,  multivariate = FALSE)$psrf[, "Point est."]
            psrf_plot <- ggplot_bars(
                as.data.frame(psrf), "psrf", bins = 10,                
                underlayers = list(
                    good = annotate("rect", 
                        xmin = -Inf, xmax = 1.05, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[2], alpha = 0.5),
                    acceptable = annotate("rect", 
                        xmin = 1.05, xmax = 1.1, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[5], alpha = 0.5),
                    good = annotate("rect", 
                        xmin = 1.1, xmax = Inf, ymin = -Inf, ymax = Inf,
                        fill = PALETTE[1], alpha = 0.5),
                    captions = labs(caption = "Green: good, Orange: acceptable, Red: bad."))
                )
            print(psrf_plot)
            interpret_diagnostics(psrf, bad=1.1, good=1.05, mode = "quick")
            standardised_ggplot_save(
                figure = psrf_plot, 
                save_path = file.path(
                    save_folder, 
                    paste0("hist_psrf_", param, ".pdf")))

            # ## Geweke diagnostic
            # cli_alert_info("-> Geweke convergence diagnostic:")
            # # Rule of thumb: <2 (no proof of non-convergence), >2 (monitor convergence closely)
            # geweke_estimates_all <- geweke.diag(chains)
            # for (i in seq(nchain(chains))) {
            #     cli_alert_info(paste0("[Chain ", i, "]"))
            #     interpret_diagnostics(geweke_estimates_all[[i]][[1]], good=2, mode = "quick")
            # }

            ## Autocorrelation
            cli_alert_info("-> Autocorrelation:")
            for (lag in c(50)) {
                autocorr_estimates <- autocorr.diag(
                    chains, lags=c(lag/thin))
                cli_alert_info(paste0("[Lag ", lag, "]:"))
                interpret_diagnostics(autocorr_estimates, good=0.1, mode = "quick")
            }
            
        }
    }
    cat("\n")
    return(NULL)
}

# A function to display XX and XY associations, as well as Variance partitions
# after fitting a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - save_folder: a string. The path where plotted PDFs will be saved.
#   - x_groups_cats: a list of integers. For each x variable, 
#       a number assigning it to a group (for variance partitioning)
#   - x_groups_names: a list of strings. 
#       A label for each number in x_groups_cat.
#   - supportLevel: a numeric between 0 and 1. 
#       The minimum confidence to display results (default is 0.95)
# RETURNS:
#   - NULL, saves plots to save_folder without printing them.
analyses_hmsc <- function(
        hM, save_folder, x_groups_cats, x_groups_names, supportLevel = 0.95) {
    # X-Y associations
    for (param in c("Beta", "Omega")) {
        if (is.null(hM$ranLevels) & (param=="Omega")){
            next
        }
        post_association = getPostEstimate(hM, parName = param)
        XY_grid <- ggplot_custom_plotBeta(
            hM, post = post_association, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = XY_grid, 
            save_path = file.path(save_folder, paste0(param, "_XY_associations.pdf")))
    }

    if (!is.null(hM$ranLevels)) {
        rand_XX_grid <- ggplot_custom_random_corr_associations(
            hM, supportLevel = supportLevel)
        standardised_ggplot_save(
            figure = rand_XX_grid, 
            save_path = file.path(save_folder, "random_XX_associations.pdf"))
    }
    
    # Variance partitionning 
    vp = computeVariancePartitioning(
        hM, 
        group = x_groups_cats, # c(1,2,2)
        groupnames = x_groups_names) # c("habitat","climate"))
    variance_bars <- ggplot_custom_plotVariancePartitioning(hM, VP = vp)
    standardised_ggplot_save(
        figure = variance_bars, 
        save_path = file.path(save_folder, "variance_partitioning.pdf"))
}

# A function that mimciks Hmsc::evaluateModelFit but can also work on 
# non-training data.
# ARGS :
#   - hM: a Hmsc fitted model object.
#   - y: a matrix of species observation ("ground truth").
#   - predY: the predictions made by the model.
# RETURNS:
#   - a list of scores
evaluateModelFitCustom <- function(hM, Y, predY) {

    ns <- ncol(Y) # number of samples per observation/species
    mPredY <- apply(predY, c(1, 2), mean)  # mean prediction per obs/species
    sdPredY <- apply(predY, c(1, 2), sd) # sd prediction per obs/species


    # Initialise metrics to compute
    RMSE <- rep(NA, ns)     # RMSE (the lower the better)
    AUC <- rep(NA, ns)      # AUC (the closer to 1, the better)
    TjurR2 <- rep(NA, ns)   # Tjur R² (% of variance explained)
    SD <- rep(NA, ns)       # Standard Deviation (only uses predicted values)

    # For each sample
    for (j in seq_len(ns)) {
        sel <- !is.na(Y[, j])      # extract observations/species
        obs <- Y[sel, j]           # get observed value
        pred <- mPredY[sel, j]     # get predicted value
        predSD <- sdPredY[sel, j]  # get predicted sd

        # compute RMSE / MSE
        RMSE[j] <- sqrt(mean((obs - pred)^2))
        

        # compute variance (only obs needed)
        SD[j] <- mean(predSD)

        # compute AUC (only meaningful if both 0s and 1s present)
        if (length(unique(obs)) == 2) {
            AUC[j] <- as.numeric(pROC::auc(obs, pred, quiet = TRUE))
        }

        # compute Tjur R2: difference in mean predicted probability between
        # presences and absences
        if (length(unique(obs)) == 2) {
            TjurR2[j] <- mean(pred[obs == 1]) - mean(pred[obs == 0])
        }
    }

    names(RMSE) <- names(SD) <- names(AUC) <- names(TjurR2) <- colnames(Y)
    return(list(RMSE = RMSE, AUC = AUC, TjurR2 = TjurR2, SD = SD))
}

# A function to display convergence diagnostics for a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - data_set: a dataframe or SpatVector. 
#       Must contain columns listed in x_cols and y_cols.
#       Automatically removes rows containing NAs.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - sp_cols: a list of strings. The columns containing species occurrences.
# RETURNS:
#   - a list of scores
evaluate_hmsc_performances <- function(hM, data_set, x_cols, sp_cols) {
    if (inherits(data_set, "SpatVector")) {
        data_set <- as.data.frame(data_set)
    } else if (!inherits(data_set, "data.frame")) {
        stop("Class of data_set is not recognised.")
    }

    if (any(is.na(data_set))) {
        # cli_alert_warning("NAs detected: rows with NAs will be excluded.")
        data_set <- data_set |> drop_na()
    }
    
    local_preds_list <- predict_hmsc(
        hM = hM, 
        df = data_set, 
        x_variables = x_cols)
    local_preds <- abind(local_preds_list, along = 3) # model predictions
    local_Y <- as.matrix(data_set[sp_cols])           # actual observations
    
    # Extract metric by comparing predictions and observed values
    evaluateModelFitCustom(hM = hM, Y = local_Y, predY = local_preds)
}