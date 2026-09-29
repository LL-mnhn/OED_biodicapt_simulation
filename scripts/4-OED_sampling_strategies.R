# This script:
#   - Generates observation datasets from simulations
#   - Fits HMSC models to observations
#   - Tries several sampling strategies 
#   - Searches for the strategy that reduces variance the most

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(dplyr)
library(geosphere) # distm function

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_model.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_simul.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

if (authorise_overwrite(RESULTS_PATH)) {
    unlink(RESULTS_PATH, recursive = TRUE)
    dir.create(RESULTS_PATH)
} else {
    stop("User refused to overwrite previous model's runs.")
}


### Variables
FORMULA <- reformulate(X_VARIABLES)


### MCMC
NSAMPLES <- 1500 # mcmc will stop after saving that much samples
THIN <- 1 # number of steps between each recording
NTRANSIENT <- 0.5*NSAMPLES*THIN # burn-in iterations
NCHAINS <- 3

### Effects to test
# Each list is a combination of parameters to use to train the model
EFFECTS <- list(
    # Reference effects (ALWAYS the first element)
    list(
        R_EFFECTS = c("none"), 
        NEW_SAMPLE_SIZES = c(0),
        STRATEGIES = c("none")
    ),
    # Add new samples according to different strategies
    list(
        R_EFFECTS = c("none"), 
        NEW_SAMPLE_SIZES = c(10, 25, 50, 100, 150),
        STRATEGIES = c(
            "business-as-usual", "gap-filling", "simplified-I-optimality", 
            "2-part-simplified-I-optimality", "5-part-simplified-I-optimality", 
            "10-part-simplified-I-optimality")
    )
)
REFERENCE_EFFECTS <- EFFECTS[1][[1]] # Always the first element (isolate it)
saveRDS(EFFECTS, file.path(RESULTS_PATH, "effects.rds"))


##### Local functions ##### ---------------------------------------------------
prepare_datasets <- function(detection_probability = 0.95) {
    # load simulated presence_absence maps (dont forget to unwrap SpatRasters)
    pa_obs <- unwrap_simulations(readRDS(
        file.path(SIMULATE_PATH, "simulations.rds")))$sp_observations

    status_msg <- cli_status(paste0(
        "[k-fold {1}/{K_FOLDS}]: Observation of sp {1}/{length(NAME_SP_SIMUL)}..."))
    
    # load observation locations
    all_datasets <- vector("list", K_FOLDS)
    for (k in 1:K_FOLDS) {
        # initialise locations (useless columns will be removed)
        train_df <- fix_names(vect(BIODICAPT_OBS_FULL))
        new_pool_df <- fix_names(vect(ENI500_OBS_FULL))
        test_df <- fix_names(vect(STOC_OBS_FULL))
        
        for (sp_i in 1:length(NAME_SP_SIMUL)) {
            cli_status_update(
                status_msg, 
                paste0(
                    "[k-fold {k}/{K_FOLDS}]: Observation of sp {sp_i}/{length(NAME_SP_SIMUL)}...")
            )
            # Use CustomSampleOccurences for fixed locations
            train_df[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                train_df,
                detection.probability = detection_probability
            )
            new_pool_df[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                new_pool_df,
                detection.probability = detection_probability
            )
            test_df[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                test_df,
                detection.probability = detection_probability
            )
        }

        # clean columns
        all_datasets[[k]]$train <- train_df |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "lon", "lat"))) |>
            mutate(row_id = paste0("id_train_", 1:nrow(train_df)))
        all_datasets[[k]]$new_pool <- new_pool_df |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "lon", "lat"))) |>
            mutate(row_id = paste0("id_new_pool_", 1:nrow(new_pool_df)))
        all_datasets[[k]]$test <- test_df |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "lon", "lat"))) |>
            mutate(row_id = paste0("id_test_", 1:nrow(test_df)))
    }

    cli_status_clear(status_msg)
    cli_alert_success("Simulated datasets are ready!")

    return(all_datasets)
}

extended_training_design <- function(
    strat, base_subset, extension_subset, new_sample_size, variables, hM) {
    # dataset changes depending on chosen strat
    if ((strat == "none") || (new_sample_size == 0)) {
        # No addition of data to base_subet
        extended_set <- base_subset
    } else if (strat == "business-as-usual") {
        # Business-as-usual: random sampling
        cli_alert_info(
            paste0("Pulling ", new_sample_size, " new random samples..."))
        extended_set <- bind_rows(
            base_subset,
            extension_subset |> slice_sample(n = new_sample_size))
        
    } else if (strat == "gap-filling") {
        # Gap-filling: select points that are the farthest from current design
        cli_alert_info(paste0("Searching for the most distant points..."))
        extended_set <- base_subset  # no need for cbind() here at all

        # Select one point at a time 
        # (re-compute distances each time a point is added)
        for (i in 1:new_sample_size) {
            dist_matrix <- distm(
                extended_set[, c("lon", "lat")],
                extension_subset[, c("lon", "lat")],
                fun = distGeo)
            idx_highest_min_distance <- order(
                apply(dist_matrix, 2, min), decreasing = TRUE)[1]
            extended_set <- bind_rows(
                extended_set,
                extension_subset[idx_highest_min_distance, ])
        }
    } else if (strat == "simplified-I-optimality") {
        # This is the method used in doi.org/10.1111/2041-210X.14355
        # In short: sample where uncertainty is the highest.
        # Here, we create a "layer" and pick the most uncertain samples at once
        cli_alert_warning("Running simplified-I-optimality...")

        # TODO: We exclude base_subset but ACTUALLY it might be interesting
        # to resample on already sampled positions to improve performance...
        uncertainty <- get_uncertainty_hmsc(
            hM = hM, 
            df = extension_subset, 
            x_cols = variables) 
        idx_most_uncertain <- order(
            uncertainty, decreasing = TRUE)[1:new_sample_size]
        extended_set <- bind_rows(
                base_subset,
                extension_subset[idx_most_uncertain, ])
    
    } else {
        stop(paste0("Unknown strategy: '", strat, "'."))
    }
    return(extended_set)
}

optimise_model_training <- function(
        training_set, remaining_pool,
        params, 
        run_k_fold, run_folder_path,
        x_variables = X_VARIABLES, y_species = NAME_SP_SIMUL,
        ref_model_effects = REFERENCE_EFFECTS) {
    # Verify data formats
    if (inherits(training_set, "SpatVector")) {
        training_set <- as.data.frame(training_set)
    } else if (!inherits(training_set, "data.frame")) {
        stop("Class of training_set is not recognised.")
    }

    if (inherits(remaining_pool, "SpatVector")) {
        remaining_pool <- as.data.frame(remaining_pool)
    } else if (!inherits(remaining_pool, "data.frame")) {
        stop("Class of remaining_pool is not recognised.")
    }
    
    
    strategy <- params$STRATEGIES
    # if strategy is X_parts, divide the number of samples and fit as many
    # times as necessary
    pattern_n_parts <- "([^_]+)-part-"
    regex_parts <- regmatches(strategy, regexec(pattern_n_parts, strategy))
    n_parts <- ifelse(
        is.na(regex_parts[[1]][2]), 1, as.numeric(regex_parts[[1]][2]))
    strategy_alone <- substring(
        strategy, 
        ifelse((n_parts == 1), 1, nchar(regex_parts[[1]][1])+1), 
        nchar(strategy))
    n_new_samples_per_part <- split_evenly(params$NEW_SAMPLE_SIZES, n_parts)
    
    if (grepl("simplified-I-optimality", strategy)) {
        # when using uncertainty model, we need to have an initial model to
        # compute the uncertainty on. For that, we can use the BASE_MODEL
        previous_model <- readRDS(
            file.path(
                make_run_path(
                    folder = dirname(run_folder_path), 
                    params = ref_model_effects,
                    k_fold = run_k_fold
                ),  
                "chains.rds"))
    } else {
        previous_model <- NULL
    }

    # If loop runs n times (n>1), the model used to compute the second part
    # of the dataset becomes the model computed on part n-1.
    # By default:
    #   - training data subset is overwritten during each loop 
    #       (as we append new points to it)
    #   - new_pool data subset is renamed as remaining_pool and overwritten 
    #       each loop (as we remove the points added to the training set)
    # This *could* be a problem if we want to sample the same points on
    # repeated experiments (e.g. during different years but at the same place). 
    for (part in 1:n_parts) {
        # 0. Setting up dataset
        training_set <- extended_training_design(
            strat = strategy_alone,
            base_subset = training_set, 
            extension_subset = remaining_pool,
            new_sample_size = n_new_samples_per_part[part],
            variables = x_variables,
            hM = previous_model)
        
        # 1. Fit model
        cat("\n")
        cli_alert_info(paste0("1.", part, ". Fitting model..."))
        base_model <- prepare_hmsc_training(
            subdataset = training_set,
            x_cols = x_variables, 
            y_cols = y_species,
            formula = FORMULA,
            random_effect = params$R_EFFECTS)
        fitted_model <- fitting_hmsc(
            hM = base_model, 
            save_to = file.path(run_folder_path, "chains.rds"),
            nchains = NCHAINS,
            thin = THIN,
            nsamples = NSAMPLES,
            ntransient = NTRANSIENT,
            freq_verbose = (NSAMPLES*2 + NTRANSIENT)/10,
            allow_parallel = TRUE)
        
        # previous model becomes the one we just computed
        previous_model <- readRDS(
            file.path(run_folder_path, "chains.rds"))
    }

    return(list(fitted_model = fitted_model, base_dataset = training_set))
}


##### Prepare datasets ##### --------------------------------------------------
param_grid <- build_param_grid(EFFECTS)
total_loops <- nrow(param_grid) * K_FOLDS
datasets <- prepare_datasets()
# save for later usage
saveRDS(
    wrap_simulations(datasets), 
    file = file.path(RESULTS_PATH, "datasets.rds"))



##### Running model ##### -----------------------------------------------------
cli_alert_info("------------ Fitting models ------------\n\n")
for (k in seq(K_FOLDS)) {
    for (p in seq_len(nrow(param_grid))) {
        cli_alert_info("---- [k-fold {k}/{K_FOLDS}]: Run {p} of {nrow(param_grid)} ----")
        
        # In each run, we use simulated data:
        #   - we train on "BIODICAPT" simulations
        #   - we use "500ENI" as a new pool of observations 
        #   - we use "STOC" simulations as test data

        # 0. Load parameters
        cli_alert_info("Loading parameters...")
        local_parameters <- param_grid[p, ]
        local_path_results <- make_run_path(
            folder = RESULTS_PATH, 
            params = local_parameters,
            k_fold = k
        )   
        if (file.exists(local_path_results)) {
            stop("Error, model folder already exists!")
        } else {
            # create
            dir.create(local_path_results, recursive = TRUE)
            # save parameters for this run (same as in name but easier to access)
            local_parameters_extended <- c(as.list(local_parameters), list(k = k))
            saveRDS(
                local_parameters_extended, 
                file.path(local_path_results, "local_parameters.rds"))

        }


        # 1. Fit model
        cli_alert_info("Fitting model...")
        opt_model <- optimise_model_training(
            training_set = datasets[[k]]$train, 
            remaining_pool = datasets[[k]]$new_pool,
            params = local_parameters, 
            run_k_fold = k, run_folder_path = local_path_results) 
        # Save extended training set 
        saveRDS(
            grep("train", opt_model$base_dataset$row_id, value = TRUE, invert = TRUE), 
            file.path(local_path_results, "ID_new_samples.rds"))
    
        # 2. Analysis of convergence
        cli_alert_info("Compiling convergence diagnostics...")
        . <- convergence_hmsc(
            hM = opt_model$fitted_model, 
            nchains = NCHAINS, 
            thin = THIN, 
            save_folder = local_path_results)
        
        cli_alert_success("Run complete!\n\n")
    }
}

cli_alert_success("Simulated species are ready!\n\n")