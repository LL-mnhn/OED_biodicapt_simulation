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
set.seed(42) # reproducibility 

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
        R_EFFECTS = c("site"),  # Best according to outputs_2026-11-16
        NEW_SAMPLE_SIZES = c(0),
        STRATEGIES = c("none")
    ),
    # Add new samples according to different strategies
    list(
        R_EFFECTS = c("site"), 
        NEW_SAMPLE_SIZES = c(10, 25, 50, 100, 150),
        STRATEGIES = c(
            "random-sampling", "by-region", "gap-filling", 
            "max-uncertainty", "ds-optimality", "i-optimality")
    )
)
REFERENCE_EFFECTS <- EFFECTS[1][[1]] # Always the first element (isolate it)
saveRDS(EFFECTS, file.path(RESULTS_PATH, "effects.rds"))


##### Local functions ##### ---------------------------------------------------
# A function that prepares datasets (samples random observations on each site)
prepare_datasets <- function() {
    # load simulated presence_absence maps (dont forget to unwrap SpatRasters)
    simulations <- unwrap_simulations(readRDS(
        file.path(SIMULATE_PATH, "simulations.rds")))
    env_pca <- readRDS(file.path(EXPLORE_PATH, "pca_data.rds"))

    status_msg <- cli_status(paste0(
        "[k-fold {1}/{K_FOLDS}]: Observation of sp {1}/{length(NAME_SP_SIMUL)}..."))
    
    # load observation locations
    all_datasets <- vector("list", K_FOLDS)
    for (k in 1:K_FOLDS) {
        # initialise locations (useless columns will be removed)
        biodicapt_vec <- fix_names(vect(BIODICAPT_OBS_TRF))
        eni500_vec <- fix_names(vect(ENI500_OBS_TRF))
        stoc_vec <- fix_names(vect(STOC_OBS_TRF))
        
        for (sp_i in 1:length(NAME_SP_SIMUL)) {
            sp_name <- paste0("sp_", sp_i)
            cli_status_update(
                status_msg, 
                paste0(
                    "[k-fold {k}/{K_FOLDS}]: Observation of sp {sp_i}/{length(NAME_SP_SIMUL)}...")
            )

            # Use CustomSampleOccurences for fixed locations
            biodicapt_vec[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                biodicapt_vec,
                detection.probability = PROBABILITY_TP,
                error.probability = PROBABILITY_FP
            )
            eni500_vec[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                eni500_vec,
                detection.probability = PROBABILITY_TP,
                error.probability = PROBABILITY_FP
            )
            stoc_vec[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
                pa_obs[[k]][[paste0("sp_", sp_i)]],
                stoc_vec,
                detection.probability = PROBABILITY_TP,
                error.probability = PROBABILITY_FP
            )

            # # Use customSampleOccurrencesLocationFree for all datasets
            # # (designed for anonymised data but works for all)
            # biodicapt_vec[[sp_name]] <- customSampleOccurrencesLocationFree(
            #     vector_stack = biodicapt_vec[X_VARIABLES], 
            #     precomputed_pca = env_pca$pca, 
            #     sp_pca_coords = simulations$sp_params[[k]][[sp_name]], 
            #     sp_pa_params = simulations$sp_pa_params[[k]][[sp_name]], 
            #     detection.probability = PROBABILITY_TP,
            #     error.probability = PROBABILITY_FP
            # )
            # eni500_vec[[sp_name]] <- customSampleOccurrencesLocationFree(
            #     vector_stack = eni500_vec[X_VARIABLES], 
            #     precomputed_pca = env_pca$pca, 
            #     sp_pca_coords = simulations$sp_params[[k]][[sp_name]], 
            #     sp_pa_params = simulations$sp_pa_params[[k]][[sp_name]], 
            #     detection.probability = PROBABILITY_TP,
            #     error.probability = PROBABILITY_FP
            # )
            # stoc_vec[[sp_name]] <- customSampleOccurrencesLocationFree(
            #     vector_stack = stoc_vec[X_VARIABLES], 
            #     precomputed_pca = env_pca$pca, 
            #     sp_pca_coords = simulations$sp_params[[k]][[sp_name]], 
            #     sp_pa_params = simulations$sp_pa_params[[k]][[sp_name]], 
            #     detection.probability = PROBABILITY_TP,
            #     error.probability = PROBABILITY_FP
            # )
        }

        # clean columns
        all_datasets[[k]]$train <- biodicapt_vec |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "lon", "lat"))) |>
            mutate(row_id = paste0("id_biodicapt_", 1:nrow(biodicapt_vec)))
        all_datasets[[k]]$new_pool <- eni500_vec |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "region", "lon", "lat"))) |>
            mutate(row_id = paste0("id_eni500_", 1:nrow(eni500_vec)))
        all_datasets[[k]]$test <- stoc_vec |>
            select(all_of(c(X_VARIABLES, NAME_SP_SIMUL, "lon", "lat"))) |>
            mutate(row_id = paste0("id_stoc_", 1:nrow(stoc_vec)))

        if( (sum(is.na(as.data.frame(all_datasets[[k]]$train))) + 
            sum(is.na(as.data.frame(all_datasets[[k]]$new_pool))) + 
            sum(is.na(as.data.frame(all_datasets[[k]]$test)))) > 0) {
            stop(paste0("[Stopping process]: NAs detected in observations."))
        }
    }

    cli_status_clear(status_msg)
    cli_alert_success("Simulated datasets are ready!")

    return(all_datasets)
}

# A function that adds new sites (taken from a pool of candidates) to a design.
# ARGS:
#   - strat: a string, the strategy used to choose new sites. One of:
#       "none", "random-sampling", "gap-filling", "max-uncertainty",
#       "ds-optimality", "i-optimality".
#   - base_subset: a dataframe, the current design (training set).
#   - extension_subset: a dataframe, the candidates (must have a `row_id`).
#   - new_sample_size: an integer, the number of sites to add.
#   - variables: a vector of strings, the names of the explanatory variables.
#   - hM: a Hmsc model fitted on base_subset (needed by max-uncertainty, 
#       ds-optimality and i-optimality, ignored otherwise).
#   - actionable: a vector of strings, the variables of interest (only used 
#       when strat="ds-optimality"). Other coefficients are nuisance parameters.
#   - target_set: a dataframe of sites where predictions should be as precise
#       as possible (only used when strat="i-optimality"). Default is NULL: 
#       the candidates are used as a proxy for the region of interest.
# RETURNS:
#   - base_subset with the selected rows of extension_subset appended.
extended_training_design <- function(
    strat, base_subset, extension_subset, new_sample_size, variables, hM,
    actionable = ACTIONABLE_VARIABLES, target_set = NULL, random_seed = 42) {
    strat <- tolower(strat)

    # dataset changes depending on chosen strat
    if ((strat == "none") || (new_sample_size == 0)) {
        # No addition of data to base_subet
        extended_set <- base_subset
    } else if (strat == "random-sampling") {
        # Random sampling: choose points at random
        cli_alert_info(
            paste0("Pulling ", new_sample_size, " new random samples..."))
        candidates <- extension_subset[
            complete.cases(extension_subset[variables]), ]
        
        # WARNING : fixed random draw everytime so that the selected sites are
        #           always the same ones. (Its ok because its not the strategy
        #           that should perform well anyway).
        idx_random <- withr::with_seed(
            random_seed, sample.int(nrow(candidates)))[seq_len(new_sample_size)]
        extended_set <- bind_rows(base_subset, candidates[idx_random, ])
        
    } else if (strat == "by-region") {
        # By region: 
        #   - pick a region
        #   - pick points at random inside it
        #   - when all points are selected, pick a new region
        #   - repeat
        idx_candidates <- c()
        for (region_name in ORDER_REGIONS) {
            if (length(idx_candidates) < new_sample_size) {
                regional_subset <- extension_subset |> 
                    filter(region == region_name)

                # select a random number of ids
                select_n_rows <- min(
                    nrow(regional_subset), 
                    new_sample_size - length(idx_candidates))
                idx_candidates <- c(
                    idx_candidates, 
                    withr::with_seed(
                        random_seed, 
                        sample(regional_subset$row_id))[seq_len(select_n_rows)]
                )
            }
        }
        extended_set <- bind_rows(
            base_subset, 
            extension_subset |> filter(row_id %in% idx_candidates))

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
    } else if (strat == "max-uncertainty") {
        # This is the method used in doi.org/10.1111/2041-210X.14355
        # In short: sample where uncertainty is the highest.
        # Here, we create a "layer" and pick the most uncertain samples at once
        cli_alert_warning("Running max-uncertainty...")

        # Warning: We exclude base_subset but, later, it might actually be 
        # interesting to resample on positions that were already sampled
        uncertainty <- get_uncertainty_hmsc(
            hM = hM, 
            df = extension_subset, 
            x_cols = variables) 
        idx_most_uncertain <- order(
            uncertainty, decreasing = TRUE)[1:new_sample_size]
        extended_set <- bind_rows(
                base_subset,
                extension_subset[idx_most_uncertain, ])
    
    } else if (strat %in% c("ds-optimality", "i-optimality")) {
        # "Greedy" selection: building the set of new sites one at a time.
        # The selection uses a criteria based on the Fisher information
        # matrix of each species, approximated from the posterior of `hM`:
        #   1. Start from the current information matrix (reference model)
        #   2. Score every remaining candidate.
        #   3. Add the best candidate and updated the information matrix
        #   4. Repeat (without redrawing of the same candidates)
        #
        # Notes:
        #   - Hmsc standardises X by default but returns Beta back on the 
        #     original scale of X, so XScale = FALSE is NOT needed. 
        #     get_hmsc_information() takes care of using X and Beta on the same
        #     scale (and of standardising them for numerical stability).
        #   - Works with no random effect and with "site" random effects (which
        #     only attenuate the information of each observation). A spatial 
        #     random effect is currently not supported (stops with an error).
        #   - Criteria are summed over species (with equal weights).
        #   - I-optimality ignores correlation between species residuals and 
        #     approximates probit with a linear model.
        if (is.null(hM)) {
            stop("Strategy '", strat, "' needs a fitted model (`hM`).")
        }
        cli_alert_info(paste0("Searching for the best ", new_sample_size, 
            " points according to ", strat, "..."))
        
        # candidates with missing covariates cannot be scored (never selected)
        candidates <- extension_subset[complete.cases(extension_subset[variables]), ]
        info <- get_hmsc_information(hM, candidates, variables)

        if (strat == "ds-optimality") {
            idx_params <- match(actionable, hM$covNames)
            if (anyNA(idx_params)) {
                stop("`actionable` contains variables absent from the model.")
            }
            idx_selected <- greedy_design(
                M0 = info$M0, X = info$X, w = info$w, 
                n_select = new_sample_size, 
                criterion = "ds", idx_params = idx_params)
        } else {
            if (is.null(target_set)) {
                target_set <- candidates
            }
            target_set <- target_set[complete.cases(target_set[variables]), ]
            idx_selected <- greedy_design(
                M0 = info$M0, X = info$X, w = info$w, 
                n_select = new_sample_size, 
                criterion = "i", B = get_prediction_moments(info, target_set))
        }

        extended_set <- bind_rows(
            base_subset,
            candidates[idx_selected, ])

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
        ref_model_effects = REFERENCE_EFFECTS,
        actionable = ACTIONABLE_VARIABLES, target_set = NULL) {
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
    
    if (    tolower(strategy_alone) %in% 
            c("max-uncertainty", "ds-optimality", "i-optimality")
    ) {
        # these strategies need an initial model to compute uncertainty / 
        # information on. For that, we can use the BASE_MODEL, fitted with the 
        # same random effects as the current run.
        previous_model <- readRDS(
            file.path(
                make_run_path(
                    folder = dirname(run_folder_path), 
                    params = modifyList(
                        ref_model_effects, 
                        list(R_EFFECTS = params$R_EFFECTS)),
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
            hM = previous_model,
            actionable = actionable,
            target_set = target_set)
        # points added to training_set must be removed from remaining_pool
        remaining_pool <- remaining_pool |>
            filter(!(row_id %in% training_set$row_id))
        
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
        # Save new_pool IDs
        saveRDS(
            grep("id_eni500", opt_model$base_dataset$row_id, value = TRUE, invert = FALSE), 
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

cli_alert_success("HMSC models are ready!\n\n")