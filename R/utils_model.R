# Set of functions used to create / handle / analyse models
##### Libraries ##### ---------------------------------------------------------
library(sf)
library(cli)
library(coda)
library(Hmsc)
library(abind)
library(tidyr)
library(dplyr)
library(abind)
library(ggplot2)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # Import global parameters


##### Data ##### --------------------------------------------------------------
# A function that converts the "lon" / "lat" columns of a dataframe (WGS84 
# degrees, EPSG:4326) to metric coordinates (PROJECT_CRS).
# ARGS:
#   - df: a dataframe with columns lon_c, lat_c and id_column.
#   - id_column: a string. The column with unique id of the points.
#   - lon_c, lat_c: strings. The names of the longitude/latitude columns.
# RETURNS:
#   - a matrix (n points x 2) of coordinates in meters.
lonlat_to_grid_xy <- function(
    df, id_column = "row_id", lon_c = "lon", lat_c = "lat") {
    # Check that the dataframe is in EPSG:4326 coordinates (approx.)
    if (!(all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[lon_c])) &
        all(as.data.frame(df[lon_c]) <= FRANCE_BBOX_4326["xmax"])   &
        all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[lat_c]))   &
        all(as.data.frame(df[lat_c]) <= FRANCE_BBOX_4326["ymax"])    
    )) {
        stop(paste0(
            "lonlat_to_grid_xy() assumes that your data ",
            "is in EPSG:4326 over France's extent. It looks like 'df' is in ",
            "another CRS, or on a different area. Please check your data."
        ))
    }
    
    xy <- st_coordinates(st_transform(
        st_as_sf(as.data.frame(df), coords = c(lon_c, lat_c), crs = 4326), 
        PROJECT_CRS))
    rownames(xy) <- as.character(df[[id_column]])
    colnames(xy) <- paste0(c("longitude_grid_", "latitude_grid_"), PROJECT_CRS_TAG)
    return(xy)
}

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
make_run_path <- function(folder, params, k_fold) {
    row_string <- tolower(
        paste(names(params), unlist(params), sep = "-", collapse = "_")
    )
    
    return(file.path(folder, paste0("model_", row_string, "_k", k_fold)))
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
#   - random_effect: adds a random effect to Hmsc (default is "none").
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
            studyDesign = NULL)
    
    } else if (random_effect == "site") {
        # Random effect as a "site" effect:
        #   An unstructured sample-unit effect with no direct spatial structure, 
        #   it only captures residual species associations.
        if (anyDuplicated(subdataset[[id_column]])) {
            stop("id_column must be unique per point for a site-level random effect.")
        }
        studyDesign <- data.frame(site = as.factor(subdataset[[id_column]]))

        rL.site <- HmscRandomLevel(units = levels(studyDesign$site))
        rL.site <- setPriors(rL.site, nfMin = 1,  nfMax = 2)

        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata,
            XFormula = formula,
            ranLevels = list(site = rL.site),
            distr = "probit",
            studyDesign = studyDesign)        
        
    } else if (random_effect == "spatial") {
        # Random effect as a "spatial" effect:
        #   The random effect is a function of the distance between the
        #   points. Takes more time to compute but (usually), yields to better 
        #   results in prediction.

        studyDesign <- data.frame(spatial = as.factor(subdataset[[id_column]]))

        # convert coordinates (lon/lat in degrees from obs are in EPSG 4326) 
        # to metric coordinates (PROJECT_CRS)
        xy <- lonlat_to_grid_xy(subdataset, id_column = id_column)

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
            rL.spatial, nfMin =1,  nfMax = 2) # limit number of latent variables

        hmsc_object <- Hmsc(
            Y = ydata, XData = xdata, 
            ranLevels = list("spatial" = rL.spatial),
            XFormula = formula, 
            distr = "probit",
            studyDesign = studyDesign)

    } else {
        stop(paste0(
            "Random effect must be one of 'none', 'site' or 'spatial'. Got ", 
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

    if ("site" %in% names(hM$ranLevels)) {
        studyDesign <- data.frame(id = as.factor(df[[id_column]]))
        names(studyDesign) <- names(hM$ranLevels)

        return(predict(
            hM, XData = XData, studyDesign = studyDesign, expected = TRUE))
        
    } else if ("spatial" %in% names(hM$ranLevels)) {
        # convert coordinates (lon/lat in degrees from obs are in EPSG 4326) 
        # to metric coordinates (PROJECT_CRS)
        xy_new <- lonlat_to_grid_xy(df, id_column = id_column)

        # Use Gradient (instead of studyDesign), for spatial gradients
        Gradient <- prepareGradient(
            hM,
            XDataNew = XData,
            sDataNew = list(spatial = xy_new)
        )

        # Make prediction on new dataset
        return(predict(hM, Gradient = Gradient, expected = TRUE))

    } else {
        # No random effect or unknow, make predictions without design
        return(predict(
            hM, XData = XData, studyDesign = NULL, expected = TRUE))
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

# A function that computes the autocorrelation of every variable of an MCMC
# object at ONE lag, averaged over chains. It gives the same values as 
# coda::autocorr.diag(chains, lags = lag), but much faster.
# ARGS:
#   - chains: a coda mcmc.list (or mcmc) object.
#   - lag: an integer. The lag, in number of recorded samples (not iterations).
# RETURNS:
#   - a named numeric vector, one autocorrelation per variable.
autocorr_at_lag <- function(chains, lag) {
    if (abs(lag - round(lag)) > 1e-8) {
        warning("Lag (", lag, ") is not a whole number of samples: rounded.")
    }
    lag <- round(lag)
    chains <- coda::as.mcmc.list(chains)

    per_chain <- lapply(chains, function(chain) {
        x <- as.matrix(chain)
        n <- nrow(x)
        if (lag >= n) stop("Lag (", lag, ") must be smaller than the number of samples (", n, ").")
        x <- sweep(x, 2, colMeans(x), "-")
        colSums(x[seq_len(n - lag), , drop = FALSE] * x[(lag + 1):n, , drop = FALSE]) / 
            colSums(x^2)
    })
    Reduce(`+`, per_chain) / length(per_chain)
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
        # print(eff_size_plot)
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
            # print(psrf_plot)
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
                autocorr_estimates <- autocorr_at_lag(chains, lag = lag / thin)
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


    # metrics on predictions
    RMSE <- rep(NA, ns)     # RMSE (the lower the better)
    AUC <- rep(NA, ns)      # AUC (the closer to 1, the better)
    TjurR2 <- rep(NA, ns)   # Tjur R² (% of variance explained)
    SD <- rep(NA, ns)       # Standard Deviation (only uses predicted values)

    # For each sample
    for (j in seq_len(ns)) {
        sel <- !is.na(Y[, j])      # extract observations/species
        obs <- Y[sel, j]           # get observed value
        pred <- mPredY[sel, j]     # get predicted value

        # compute RMSE / MSE
        RMSE[j] <- sqrt(mean((obs - pred)^2))       

        # compute variance (only obs needed)
        SD[j] <- mean(sdPredY[sel, j])

        # If presence and absence data:
        #   - compute AUC
        #   - compute TjurR2: difference in mean predicted probability between
        #                     presences and absences
        if (length(unique(obs)) == 2) {
            AUC[j] <- as.numeric(pROC::auc(obs, pred, quiet = TRUE))
            TjurR2[j] <- mean(pred[obs == 1]) - mean(pred[obs == 0])
        }
    }

    names(RMSE) <- names(SD) <- names(AUC) <- names(TjurR2) <- colnames(Y)
    list(RMSE = RMSE, AUC = AUC, TjurR2 = TjurR2, SD = SD)
}

# A function that extracts the posterior draws of Beta from a Hmsc model and
# puts them, together with the design matrix, on a standardised scale of X
# (mean 0 and sd 1 for each covariate, as seen by the model at fitting time).
# By defaylt, Hmsc returns Beta on the ORIGINAL scale of X.
# ARGS:
#   - hM: a fitted Hmsc object.
#   - x_cols: a vector of strings, the names of the explanatory variables.
#   - scaling_data: a dataframe with columns x_cols, or NULL (default). The 
#       means and SDs used to standardise X are computed on this dataset. When 
#       NULL, the training data of the model are used. Use a COMMON dataset to
#       compare models trained on different sites: the value of Ds, for 
#       example, shifts by a constant that depends on the scaling used.
# RETURNS: a list
#   - beta_std: array (p x n_species x n_draws) of posterior draws.
#   - standardise: function turning a dataframe (with columns x_cols) into a
#       standardised design matrix (n sites x p).
#   - post: the pooled posterior samples (see Hmsc::poolMcmcChains).
get_standardised_posterior <- function(hM, x_cols, scaling_data = NULL) {
    post <- Hmsc::poolMcmcChains(hM$postList)
    beta_orig <- simplify2array(lapply(post, function(s) s$Beta))
    p <- dim(beta_orig)[1]
    n_draws <- dim(beta_orig)[3]
    if (n_draws <= p) {
        stop("Not enough posterior samples (", n_draws, ") to estimate a ", 
             p, " x ", p, " posterior covariance.")
    }

    # Centering / scaling of X (from the training data of the model)
    i_int <- which(hM$covNames == "(Intercept)")
    X_scaling <- if (is.null(scaling_data)) {
        hM$X
    } else {
        model.matrix(hM$XFormula, data = as.data.frame(scaling_data[x_cols]))
    }
    mu <- colMeans(X_scaling)
    sdv <- apply(X_scaling, 2, sd)
    mu[i_int] <- 0
    sdv[i_int] <- 1
    sdv[sdv == 0] <- 1
    standardise <- function(df) {
        X <- model.matrix(hM$XFormula, data = as.data.frame(df[x_cols]))
        stopifnot(identical(colnames(X), hM$covNames))
        sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
    }

    beta_std <- beta_orig * sdv        # recycled along the first dimension
    beta_std[i_int, , ] <- beta_orig[i_int, , ] + 
        apply(beta_orig[-i_int, , , drop = FALSE] * mu[-i_int], c(2, 3), sum)

    list(beta_std = beta_std, standardise = standardise, post = post)
}

# A function that measures, for each species, how precisely a fitted Hmsc 
# model knows its parameters and its predictions. All of them are computed
# from the posterior of Beta on the standardised scale of X (see 
# get_standardised_posterior), and are therefore comparable between models.
# ARGS:
#   - hM: a fitted Hmsc object.
#   - data_set: a dataframe with columns x_cols (no NAs).
#   - x_cols: a vector of strings, the names of the explanatory variables.
#   - actionable: a vector of strings, the parameters of interest 
#       (subset of x_cols).
#   - scaling_data: see get_standardised_posterior(). Pass the same dataset for
#       all the models you want to compare.
# RETURNS:
#   - a list of named vectors (one value per species).
#       - BetaSD: mean posterior SD of all coefficients (including the intercept).
#       - BetaSD_actionable: same, for the actionable coefficients only.
#       - Ds: Ds-criterion on the actionable coefficients, log det of the Schur 
#             complement of the posterior precision matrix of Beta (= minus the log 
#             det of the posterior covariance of the actionable coefficients).
#             The HIGHER the better. Does not depend on data_set.
#       - I: I-criterion, average over the sites of data_set of the (delta-method)
#            variance of the predicted probability, phi(eta)^2 x' Cov(Beta) x. 
#            The LOWER the better.
get_posterior_design_metrics <- function(
        hM, data_set, x_cols, actionable, scaling_data = NULL) {
    idx_act <- match(actionable, hM$covNames)
    if (anyNA(idx_act)) {
        stop("`actionable` contains variables absent from the model.")
    }

    sp <- get_standardised_posterior(hM, x_cols, scaling_data)
    beta_std <- sp$beta_std
    ns <- dim(beta_std)[2]
    X <- sp$standardise(data_set)
    eta <- X %*% apply(beta_std, c(1, 2), mean)

    BetaSD <- BetaSD_actionable <- Ds <- I <- rep(NA_real_, ns)
    for (j in seq_len(ns)) {
        post_cov <- cov(t(beta_std[, j, ]))
        BetaSD[j] <- mean(sqrt(diag(post_cov)))
        BetaSD_actionable[j] <- mean(sqrt(diag(post_cov)[idx_act]))
        Ds[j] <- -as.numeric(determinant(
            post_cov[idx_act, idx_act, drop = FALSE], logarithm = TRUE)$modulus)
        Xw <- X * dnorm(eta[, j])
        I[j] <- -mean(rowSums((Xw %*% post_cov) * Xw))
    }
    out <- list(
        BetaSD = BetaSD, BetaSD_actionable = BetaSD_actionable, Ds = Ds, I = I)
    lapply(out, function(v) setNames(v, hM$spNames))
}

# Ratio of posterior SDs of Beta between a retrained model and a reference 
# model, on the standardised X scale.
# ARGS:
#   - hM_new, hM_ref: fitted Hmsc objects (same species, same covariates).
#   - x_cols: the explanatory variables (see get_standardised_posterior).
#   - actionable: the actionable parameters (subset of x_cols).
#   - scaling_data: a COMMON dataset for standardising X in both models.
# RETURNS: a list of named vectors (one value per species):
#   - BetaSD_ratio: geometric mean over ALL coefficients (incl. intercept) of 
#       SD_new / SD_ref. 
#   - BetaSD_actionable_ratio: same, over the actionable coefficients only.
#   1 = no change; < 1 = the new model is more precise than the reference.
get_beta_sd_ratio <- function(
        hM_new, hM_ref, x_cols, actionable, scaling_data = NULL) {
    idx_act <- match(actionable, hM_new$covNames)
    if (anyNA(idx_act)) {
        stop("`actionable` contains variables absent from the model.")
    }
    stopifnot(
        identical(hM_new$covNames, hM_ref$covNames),
        identical(hM_new$spNames, hM_ref$spNames))

    posterior_sd <- function(hM) { # matrix (p x n_species)
        beta_std <- get_standardised_posterior(hM, x_cols, scaling_data)$beta_std
        apply(beta_std, c(1, 2), sd)
    }
    log_ratio <- log(posterior_sd(hM_new) / posterior_sd(hM_ref))

    out <- list(
        BetaSD_ratio = exp(colMeans(log_ratio)),
        BetaSD_actionable_ratio = exp(
            colMeans(log_ratio[idx_act, , drop = FALSE])))
    lapply(out, function(v) setNames(v, hM_new$spNames))
}

# A function to display convergence diagnostics for a Hmsc model.
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - data_set: a dataframe or SpatVector. 
#       Must contain columns listed in x_cols and y_cols.
#       Automatically removes rows containing NAs.
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - sp_cols: a list of strings. The columns containing species occurrences.
#   - actionable: a vector of strings, the parameters of interest (subset of 
#       x_cols), used for the "actionable" scores. 
#   - scaling_data: a dataframe (or NULL). Dataset used to standardise X for the
#       parameter-based scores (BetaSD, Ds, I). Use the same one for all the 
#       models that will be compared (see get_standardised_posterior).
# RETURNS:
#   - a list of scores (one named vector, with one value per species, per score)
evaluate_hmsc_performances <- function(
        hM, data_set, x_cols, sp_cols, actionable = ACTIONABLE_VARIABLES,
        scaling_data = NULL) {
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
    c(
        evaluateModelFitCustom(
            hM = hM, Y = local_Y, predY = local_preds),
        get_posterior_design_metrics(
            hM = hM, data_set = data_set, x_cols = x_cols, 
            actionable = actionable, scaling_data = scaling_data)
    )
}

# A function that extracts everything needed to score candidate sites with a 
# Fisher-information based criterion from a HMSC model. 
# Everything is computed on a standardised X (mean 0, sd 1 for each covariate).
# ARGS:
#   - hM: a fitted Hmsc object (probit model, no spatial random level).
#   - newdata: a dataframe of candidate sites with columns listed in x_cols.
#       Must not contain NAs in x_cols.
#   - x_cols: a vector of strings, the names of the explanatory variables.
#   - scaling_data: a dataframe with columns x_cols, or NULL (default). The 
#       means and SDs used to standardise X are computed on this dataset. When 
#       NULL, the training data of the model are used. Use a COMMON dataset to
#       compare models trained on different sites: the value of Ds, for 
#       example, shifts by a constant that depends on the scaling used.
# RETURNS: a list
#   - X: matrix (n_candidates x p), standardised design matrix of candidates.
#   - M0: array (p x p x n_species), current information matrices.
#   - w: matrix (n_candidates x n_species), information weights.
#   - beta_hat: matrix (p x n_species), posterior mean of Beta (standardised X).
get_hmsc_information <- function(hM, newdata, x_cols, scaling_data = NULL) {
    if (!inherits(hM, "Hmsc")) {
        stop("'hM' must be a fitted Hmsc object.")
    }
    # Spatial random levels need the full covariance between sites: not implemented
    for (r in seq_len(hM$nr)) {
        if (isTRUE(hM$rL[[r]]$sDim > 0)) {
            stop("Information-based criteria are implemented for models with no ",
                 "random effect or non-spatial ('site') random effects only.")
        }
    }

    # Posterior of Beta on the standardised scale of X
    sp <- get_standardised_posterior(hM, x_cols, scaling_data)
    beta_std <- sp$beta_std
    standardise <- sp$standardise
    p <- dim(beta_std)[1]
    ns <- dim(beta_std)[2]
    beta_hat <- apply(beta_std, c(1, 2), mean)
    post <- sp$post

    # Current information = posterior precision of each species' coefficients
    M0 <- array(NA_real_, dim = c(p, p, ns))
    for (j in seq_len(ns)) {
        post_cov <- cov(t(beta_std[, j, ]))
        M0[, , j] <- solve(post_cov + diag(1e-10 * mean(diag(post_cov)), p))
    }

    # Weights of candidates, attenuated by non-spatial random effects
    X <- standardise(newdata)
    w <- probit_weights(X %*% beta_hat)
    g <- rep(0, ns)
    for (r in seq_len(hM$nr)) {
        g <- g + rowMeans(sapply(post, function(s) {
            colSums(matrix(s$Lambda[[r]], ncol = ns)^2)
        }), na.rm = TRUE)
    }
    w <- 1 / (1 / w + matrix(g, nrow(w), ns, byrow = TRUE))

    list(
        X = X, 
        M0 = M0, 
        w = w, 
        beta_hat = beta_hat, 
        standardise = standardise,
        spNames = hM$spNames)
}

# A to compute weights for a probit regression model
# ARGS:
#   - eta: a vector of numerics. Typically \(\eta = X\beta\), 
#          with \(X\) the design matrix and \(\beta\) the coefficient vector.
#   - eps: a very small positive float. used to avoid solving problems.
# RETURNS:
#   - a vector of weights
probit_weights <- function(eta, eps = 1e-8) {
    p <- pmin(pmax(pnorm(eta), eps), 1 - eps)
    dnorm(eta)^2 / (p * (1 - p))
}

# A function that computes the matrices B_j for the I-criterion: the average, 
# over a target set of sites, of phi(eta)^2 x x' (delta-method gradient of the 
# predicted probability), one matrix per species.
# ARGS:
#   - info: the output of get_hmsc_information().
#   - target: a dataframe of target sites (the "design space": where we want 
#       good predictions), with the same columns as the candidates.
# RETURNS:
#   - an array (p x p x n_species)
get_prediction_moments <- function(info, target) {
    Xt <- info$standardise(target)
    eta <- Xt %*% info$beta_hat
    p <- ncol(Xt)
    ns <- ncol(eta)
    B <- array(NA_real_, dim = c(p, p, ns))
    for (j in seq_len(ns)) {
        Xw <- Xt * dnorm(eta[, j])
        B[, , j] <- crossprod(Xw) / nrow(Xt)
    }
    B
}

# A function that *greedily* selects sites (aka selects one at a time) 
# maximising an information criterion.
# At each step, the exact gain of every remaining candidate is computed in 
# closed form (Sherman-Morrison), the best one is added, and the information 
# matrices are updated. Gains are summed over species.
#   - Ds: gain = -log(1 - w a_S' S a_S / (1 + w x'M^-1 x)), 
#       S = Schur complement of the parameters of interest.
#   - I:  gain = w a' B a / (1 + w x'M^-1 x) (decrease of tr(M^-1 B)).
#   with a = M^-1 x.
# ARGS:
#   - M0: array (p x p x n_species) of current information matrices.
#   - X: matrix (n_candidates x p) of candidates' (standardised) covariates.
#   - w: matrix (n_candidates x n_species) of information weights.
#   - n_select: an integer. The number of sites to select.
#   - criterion: "ds" or "i".
#   - idx_params: (Ds only) integer positions, among the columns of X, of the 
#       parameters of interest (all other columns are nuisance parameters).
#   - B: (I only) array (p x p x n_species) from get_prediction_moments().
#   - species_weights: a vector (length n_species) of weights, default is 
#       equal weights for all species.
# RETURNS:
#   - the row indices (in X) of the selected candidates, in selection order.
greedy_design <- function(
        M0, X, w, n_select, criterion = c("ds", "i"), 
        idx_params = NULL, B = NULL, species_weights = NULL) {
    
    criterion <- match.arg(criterion)
    ns <- dim(M0)[3]
    if (is.null(species_weights)) species_weights <- rep(1 / ns, ns)
    if (n_select > nrow(X)) stop("Cannot select more sites than candidates.")
    if (criterion == "ds") {
        if (length(idx_params) == 0) stop("Ds-optimality needs `idx_params`.")
    } else if (is.null(B)) {
        stop("I-optimality needs `B`.")
    }

    M <- M0
    cand_left <- seq_len(nrow(X))
    selected <- integer(0)

    for (k in seq_len(n_select)) {
        Xc <- X[cand_left, , drop = FALSE]
        gain <- numeric(length(cand_left))

        for (j in seq_len(ns)) {
            Minv <- solve(M[, , j])
            br <- Minv %*% t(Xc)                                                 # M^-1 * r
            q <- 1 / w[cand_left, j] + colSums(t(Xc) * br)                       # 1/w + x'M^-1x
            if (criterion == "ds") {
                schur_comp <- solve(Minv[idx_params, idx_params, drop = FALSE]) # M_{1|2}
                a <- br[idx_params, , drop = FALSE]                              # M^{-1} x
                ratio <- colSums(a * (schur_comp %*% a)) / q
                gain_j <- -log(pmax(1 - ratio, 1e-12))
            } else {
                gain_j <- colSums(br * (B[, , j] %*% br)) / q
            }
            gain <- gain + species_weights[j] * gain_j
        }

        best <- which.max(gain)
        x_best <- X[cand_left[best], ]
        for (j in seq_len(ns)) {
            M[, , j] <- M[, , j] + w[cand_left[best], j] * tcrossprod(x_best)
        }
        selected <- c(selected, cand_left[best])
        cand_left <- cand_left[-best]
    }
    selected
}

# A function that computes the *expected* (approximated) Ds and I criteria of a 
# design: the values predicted by the reference model once a set of sites is 
# added to its training data. These are the quantities that greedy_design() 
# optimises, as opposed to the "realised" ones measured on the re-trained model 
# (see get_posterior_design_metrics).
# ARGS:
#   - info: output of get_hmsc_information(), computed on the added sites 
#       (so that info$X and info$w describe them), with the reference model.
#   - actionable: a vector of strings, the parameters of interest for the 
#       Ds-criterion.
#   - target: a dataframe of target sites for the I-criterion (columns x_cols, 
#       no NAs), or NULL to skip I. Use the dataset on which the realised I 
#       is computed to compare them.
# RETURNS:
#   - a list of named vectors (one value per species): Ds_expected, I_expected.
get_expected_design_metrics <- function(info, actionable, target = NULL) {
    idx_act <- match(actionable, colnames(info$X))
    if (anyNA(idx_act)) {
        stop("`actionable` contains variables absent from the model.")
    }
    ns <- dim(info$M0)[3]
    B <- if (!is.null(target)) get_prediction_moments(info, target) else NULL

    Ds <- I <- rep(NA_real_, ns)
    for (j in seq_len(ns)) {
        # M_new = M0 + sum_i w_ij x_i x_i'
        M_new <- info$M0[, , j] + crossprod(info$X, info$X * info$w[, j])
        Minv <- solve(M_new)
        Ds[j] <- as.numeric(determinant(
            solve(Minv[idx_act, idx_act, drop = FALSE]), 
            logarithm = TRUE)$modulus)
        if (!is.null(B)) {
            I[j] <- -sum(Minv * B[, , j])   # tr(Minv B), both symmetric
        }
    }
    out <- list(Ds_expected = Ds)
    if (!is.null(B)) out$I_expected <- I
    lapply(out, function(v) setNames(v, info$spNames))
}


##### KL Divergence ##### -----------------------------------------------------
# A function that converts probabilities to the logit scale, after clamping
# them away from 0 and 1 (which would give infinite values).
# ARGS:
#   - p: a numeric (vector, matrix or array). Probabilities in [0, 1].
#   - eps: a numeric. Probabilities are clamped to [eps, 1 - eps].
# RETURNS:
#   - A numeric of the same shape as p. log(p / (1 - p)) of the clamped values.
logit_safe <- function(p, eps = 1e-6) {
    p <- pmin(pmax(p, eps), 1 - eps)
    log(p / (1 - p))
}

# A function that computes the closed-form KL divergence between two univariate
# Gaussians, KL(P || Q) with P = N(mu_p, sd_p^2) and Q = N(mu_q, sd_q^2).
# ARGS:
#   - mu_p, sd_p: numerics. Mean and standard deviation of P.
#   - mu_q, sd_q: numerics. Mean and standard deviation of Q.
#       (sd_p and sd_q must be strictly positive.)
# RETURNS:
#   - A numeric of the same shape as the inputs. The KL(P || Q), in nats.
kl_gaussian <- function(mu_p, sd_p, mu_q, sd_q) {
    log(sd_q / sd_p) + (sd_p^2 + (mu_p - mu_q)^2) / (2 * sd_q^2) - 0.5
}

# A function that computes the KL divergence between two Bernoulli
# distributions, KL(Bern(p) || Bern(q)).
# Compared to kl_gaussian, it ignores uncertainty, but needs no distributional 
# assumption, so it is a good robustness check.
# ARGS:
#   - p: a numeric. Success probability of the first distribution (P).
#   - q: a numeric. Success probability of the second distribution (Q).
#   - eps: a numeric. Probabilities are clamped to [eps, 1 - eps] to avoid
#       infinite values when p or q is exactly 0 or 1.
# RETURNS:
#   - A numeric of the same shape as the inputs. The KL(P || Q), in nats.
kl_bernoulli <- function(p, q, eps = 1e-6) {
    p <- pmin(pmax(p, eps), 1 - eps)
    q <- pmin(pmax(q, eps), 1 - eps)
    p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
}

# A function that summarises the posterior distribution of the predictions of
# an Hmsc model, for each (site, species) pair. Used to avoid storing every
# MCMC draw (and to compute the reference model's summary only once).
# ARGS:
#   - preds_list: a list of length n_draws, each element a matrix of
#       n_sites x n_species (output of predict_hmsc()).
#   - eps: a numeric. Probabilities are clamped before the logit
#       transformation (see logit_safe()).
# RETURNS:
#   - A list of three n_sites x n_species matrices:
#       - p_mean: posterior mean of the predicted probability.
#       - lg_mean: posterior mean of the predicted probability, logit scale.
#       - lg_sd: posterior standard deviation, logit scale.
summarise_prediction_posterior <- function(preds_list, eps = 1e-6) {
    arr <- abind(preds_list, along = 3)          # site x species x draws
    lg <- logit_safe(arr, eps)
    list(
        p_mean = apply(arr, c(1, 2), mean),
        lg_mean = apply(lg, c(1, 2), mean),
        lg_sd = apply(lg, c(1, 2), sd)
    )
}

# A function that computes the KL divergence between the posterior predictions
# of a reference Hmsc model and a retrained Hmsc model, on a common set of
# evaluation sites. D_KL are averaged over sites, giving one value per species.
# ARGS:
#   - hM_ref: a fitted Hmsc model object. The reference model.
#   - hM_new: a fitted Hmsc model object. The retrained model (same
#       random-effect structure as hM_ref).
#   - df: a dataframe of evaluation sites, with columns lon, lat, row_id
#       and names in x_cols.
#   - x_cols: a list of strings. The names of the explanatory variables.
#   - direction: a string. D_KL is directional, use "new_vs_ref" for 
#       KL(new || ref), "ref_vs_new" for KL(ref || new), and "symmetric" for 
#       the mean of both. Outputs all three by defaults.
#   - sd_floor: a numeric. Lower bound on logit-scale standard deviations
#       (avoids infinite KL when a posterior is almost degenerate).
#   - eps: a numeric. Clamping of probabilities (see logit_safe()).
#   - ref_summary: NULL or the output of summarise_prediction_posterior() for
#       hM_ref on df. The reference is the same for all runs of a given fold:
#       compute it once and pass it here to save time.
# RETURNS:
#   - A named list of two named numeric vectors (one value per species), in
#       the same format as evaluateModelFitCustom():
#       - KL_logit_gauss
#       - KL_bernoulli
kl_between_models <- function(
        hM_ref, hM_new, df, x_cols,
        direction = c("new_vs_ref", "ref_vs_new", "symmetric"),
        sd_floor = 1e-3, eps = 1e-6, ref_summary = NULL) {
    direction <- match.arg(direction)

    # Check that models used a "probit" distribution, this function was not
    # made to handle other distributions.
    if (hM_ref$call$distr != "probit") {
        stop("kl_between_models can only handle HMSC models that used parameter `distr = 'probit'`")
    }

    if (is.null(ref_summary)) {
        ref_summary <- summarise_prediction_posterior(
            predict_hmsc(hM = hM_ref, df = df, x_variables = x_cols), eps)
    }
    new_summary <- summarise_prediction_posterior(
        predict_hmsc(hM = hM_new, df = df, x_variables = x_cols), eps)

    sd_ref <- pmax(ref_summary$lg_sd, sd_floor)
    sd_new <- pmax(new_summary$lg_sd, sd_floor)

    kl_g_new_ref <- kl_gaussian(
        new_summary$lg_mean, sd_new, ref_summary$lg_mean, sd_ref)
    kl_g_ref_new <- kl_gaussian(
        ref_summary$lg_mean, sd_ref, new_summary$lg_mean, sd_new)
    kl_b_new_ref <- kl_bernoulli(new_summary$p_mean, ref_summary$p_mean, eps)
    kl_b_ref_new <- kl_bernoulli(ref_summary$p_mean, new_summary$p_mean, eps)

    pick <- function(new_ref, ref_new) {
        switch(direction,
            new_vs_ref = new_ref,
            ref_vs_new = ref_new,
            symmetric = (new_ref + ref_new) / 2)
    }
    kl_g <- pick(kl_g_new_ref, kl_g_ref_new)
    kl_b <- pick(kl_b_new_ref, kl_b_ref_new)

    sp_names <- hM_new$spNames
    list(
        KL_logit_gauss = setNames(colMeans(kl_g), sp_names),
        KL_bernoulli = setNames(colMeans(kl_b), sp_names)
    )
}