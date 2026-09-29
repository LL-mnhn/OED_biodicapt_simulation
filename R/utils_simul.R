# Set of utility functions used to help species simulations (virtualspecies)

##### Liraries ##### ----------------------------------------------------------
library(cli)
library(terra)


##### Global functions ##### --------------------------------------------------
# A function to "wrap" SpatRasters, which enables saving them.
# By default, they are stored using C++ external pointers 
# which cannot be accessed after saving.
# ARGS:
#   - x: a list or an object containing SpatRasters
# RETURNS:
#   - the same object with PackedSpatRaster replacing SpatRasters
wrap_simulations <- function(simulations) {
  rapply(
    simulations,
    function(x) if (inherits(x, "SpatRaster")) terra::wrap(x) else x,
    how = "replace"
  )
}

# A function to "unwrap" SpatRasters which were saved to PackedSpatRaster.
# By default, they are stored using C++ external pointers 
# which cannot be accessed after saving.
# ARGS:
#   - x: a list or an object containing PackedSpatRaster
# RETURNS:
#   - the same object with SpatRasters replacing PackedSpatRaster
unwrap_simulations <- function(simulations_wrapped) {
  rapply(
    simulations_wrapped,
    function(x) if (inherits(x, "PackedSpatRaster")) terra::unwrap(x) else x,
    how = "replace"
  )
}


# A function that mimics virtualspecies::sampleOccurrences by sampling 
# occurences from fixed positions given by a SpatVector, instead of sampling
# random locations.
# ARGS:
#   - vsp_map: a SpatRaster object or the output list from generateSpFromFun,
#       generateSpFromPCA, generateRandomSp, convertToPA or limitDistribution.
#       The raster must contain values of 0 or 1 (or NA).
#   - df: a SpatVector containing the point location to be sampled.
#   - detection.probability: a numeric value between 0 and 1, 
#       corresponding to the probability of detection of the species.
# RETURNS:
#   - a vector of observations/detections (contains values of 0/NA or 1)
customSampleOccurrences <- function(vsp_map, df, detection.probability) {
    if (!inherits(df, "SpatVector")) {
        stop("'df' is not a SpatVector object.")
    }

    # get presence absence
    if (inherits(vsp_map, "list") && inherits(vsp_map, "virtualspecies")) {
        pa_rast <- vsp_map$pa.raster 
    } else if (inherits(vsp_map, "SpatRaster")) {
        pa_rast <- vsp_map
    } else {
        stop("'vsp_map' is not recognised.")
    }
    

    # extract PA at given points
    pa_vals <- terra::extract(pa_rast, df)[, 2]

    # add detection errors afterwards
    detected <- ifelse(
        pa_vals == 1,
        rbinom(length(pa_vals), 1, prob = detection.probability),
        0
    )

    return(detected)
}

# A function that generates random observations of virtualspecies generated
# using a set of environmental variables (summarised in a PCA) and prevalences.
# ARGS:
#   - raster_stack: a multi-layer SpatRaster of environmental variables.
#       Each layer must be a numeric.
#   - precomputed_pca: a dudi.pca object (computed from raster_stack values).
#   - n_sp: an integer >0. The number of species to simulate.
#   - prevalences: a set of values between 0 and 1. The target frequency of 
#       occurence for each generated species. 
#       Must have length(prevalences) = n_sp.
# RETURNS: a list
#   - sp_params: means and sds for the ellipsoids of each species
#   - k_sp_obs: 
simulate_from_PCA <- function(
        raster_stack, precomputed_pca, n_sp, prevalences, 
        k_folds = K_FOLDS, seed = 42) {
    # Check parameters
    
    set.seed(seed)
    k_simulations <- list(sp_params = list(), sp_observations = list())
    status_msg <- cli_status(
        "[k-fold {1}/{k_folds}]: Generating species {1} of {n_sp}...")

    for (k in 1:k_folds) {
        k_sp_params <- list()
        k_sp_obs <- list()

        for (sp_i in 1:n_sp) {
            cli_status_update(
                status_msg, 
                "[k-fold {k}/{k_folds}]: Generating species {sp_i} of {n_sp}...")

            # Pick random coordinates in PCA to generate sp suitability
            suitability <- suppressMessages(generateSpFromPCA(
                raster.stack = raster_stack,
                pca = precomputed_pca,
                niche.breadth = "narrow",
                sample.points = TRUE,
                nb.points = 25000, # default is 10000
                plot = FALSE
            ))

            # Store the coordinates that were picked for the ellipsoid
            k_sp_params[[paste0("sp_", sp_i)]] <- list(
                means = suitability$details$means,
                sds = suitability$details$sds
            )

            # Determine presence-absence based on target species prevalence
            # (virtualspecies auto solves alpha and beta to match prevalence)
            k_sp_obs[[paste0("sp_", sp_i)]] <- suppressMessages(convertToPA(
                suitability,
                PA.method = "probability",
                prob.method = "logistic",
                species.prevalence = prevalences[sp_i], 
                plot = FALSE))$pa.raster
            
        }
        k_simulations$sp_params[[k]] <- k_sp_params
        k_simulations$sp_observations[[k]] <- k_sp_obs
    }

    cli_status_clear(status_msg)
    cli_alert_success("Simulated species are ready!")

    return(k_simulations)
}
