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
wrap_simulations <- function(data) {
  rapply(
    data,
    function(x) if (inherits(x, c("SpatRaster", "SpatVector"))) terra::wrap(x) else x,
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
unwrap_simulations <- function(data_wrapped) {
  rapply(
    data_wrapped,
    function(x) if (inherits(x, c("PackedSpatRaster", "PackedSpatVector"))) terra::unwrap(x) else x,
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
#   - error.probability: a numeric value between 0 and 1, 
#       corresponding to the probability of detecting a species that is not 
#       there in reality.
# RETURNS:
#   - a vector of observations/detections (contains values of 0/NA or 1)
customSampleOccurrences <- function(
    vsp_map, df, detection.probability, error.probability) {
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
        rbinom(length(pa_vals), 1, prob = error.probability)
    )

    return(detected)
}

# A function that mimics virtualspecies::convertToPA but uses a specific set
# of sites to compute prevalences, isntead of the whole raster
# ARGS:
#   - vspsuit_rast_map: a SpatRaster object.
#       A suitability layer from virtualspecies (0-1 values, and NAs)
#   - sites: a SpatVector containing the point location to be sampled.
#   - target: a numeric between 0 and 1, the target proportion.
#   - alpha: a numeric, defines the slope of the logistic curve.
#   - detection: a numeric value between 0 and 1 corresponding to the 
#       probability of detection of the species.
# RETURNS:
#   - a vector of observations/detections (contains values of 0/NA or 1)
calibrate_beta_to_sites <- function(
        suit_rast, 
        sites, 
        target,
        alpha = -0.05, 
        detection = 1) {
    s <- terra::extract(suit_rast, sites)[, 2]
    s <- s[!is.na(s)]
    f <- function(beta) mean(1 / (1 + exp((s - beta) / alpha))) * detection - target
    uniroot(f, interval = c(-1, 2))$root
}

# A function that generates random observations of virtualspecies generated
# using a set of environmental variables (summarised in a PCA) and prevalences.
# ARGS:
#   - raster_stack: a multi-layer SpatRaster of environmental variables.
#       Each layer must be a numeric.
#   - precomputed_pca: a dudi.pca object (computed from raster_stack values).
#   - n_sp: an integer >0. The number of species to simulate.
#   - sites: a SpatVector, the locations to sample prevalences.
#       Default is NULL (uses the whole raster).
#   - prevalences: a set of values between 0 and 1. The target frequency of 
#       occurence for each generated species. 
#       Must have length(prevalences) = n_sp.
# RETURNS: a list
#   - sp_params: means and sds for the ellipsoids of each species
#   - k_sp_obs: 
simulate_from_PCA <- function(
        raster_stack, precomputed_pca, n_sp, prevalences, sites = NULL,
        k_folds = K_FOLDS, max_draws = 99, min_suit = 0.25, seed = 42) {
    # Global parameter
    alpha = -0.05
    
    # Check parameters
    set.seed(seed)
    k_simulations <- list(
        sp_params = list(), 
        sp_suitability = list(),
        sp_probability = list(),
        sp_pa_params = list(),
        sp_observations = list())
    status_msg <- cli_status(
        "[k-fold {1}/{k_folds}]: Generating species {1} of {n_sp}...")

    for (k in 1:k_folds) {
        k_sp_params <- list()
        k_sp_suit <- list()
        k_sp_prob <- list()
        k_sp_conv <- list()
        k_sp_obs <- list()
        
        for (sp_i in 1:n_sp) {
            cli_status_update(
                status_msg, 
                "[k-fold {k}/{k_folds}]: Generating species {sp_i} of {n_sp}...")

            # With a set of given prevalences, suitability has to reach
            # a minimum level, or it will just be a random draw of presence-
            # absence in each cell, and the map will just be noise.
            target <- prevalences[sp_i]
            accepted <- FALSE
            draw <- 0
            beta <- NA
            while (!accepted && (draw < max_draws)) {
                accepted <- TRUE
                draw <- draw + 1
                cli_status_update(
                    status_msg, 
                    "[k-fold {k}/{k_folds}]: Generating species {sp_i} of {n_sp} (draw: {draw}/{max_draws})...")

                # Pick random coordinates in PCA to generate sp suitability
                suitability <- suppressMessages(generateSpFromPCA(
                    raster.stack = raster_stack,
                    pca = precomputed_pca,
                    niche.breadth = "any",
                    sample.points = TRUE,
                    nb.points = 25000, # default is 10000
                    plot = FALSE
                ))

                # Suitability values where prevalence is evaluated
                if (is.null(sites)) {
                    s <- as.vector(terra::values(suitability$suitab.raster))
                } else {
                    s <- terra::extract(suitability$suitab.raster, sites)[, 2]
                }
                s <- s[!is.na(s)]

                # If suitability of the marginal cell (the least suitable cell
                # that can still qualify as "presence" data under the chosen 
                # target prevalence) is inferior to threshold, re-draw
                if (quantile(s, 1 - target) < min_suit) {
                    accepted <- FALSE
                    next
                }

                # If the proportion of cells with high suitability (> 50%) is
                # clearly larger than the target prevalence needs, then the 
                # distribution risks being very sparse, re-draw.
                if (mean(s > 0.5) > 3 * target) {    
                    accepted <- FALSE
                    next
                }

                # Compute beta manually for sites
                if (!is.null(sites)) {
                    # Use a custom function that computes prevalence on sites
                    # An infeasible niche makes uniroot fail: just redraw
                    beta <- tryCatch(
                        calibrate_beta_to_sites(
                            suit_rast = suitability$suitab.raster, 
                            sites = sites, 
                            target = target,
                            alpha = alpha),
                        error = function(e) NA)
                    accepted <- !is.na(beta)
                }
            }
            if (!accepted) {
                cli_alert_warning("[k-fold {k}/{k_folds}]: Exceeded max number of draws for {sp_i} of {n_sp}.")
                stop("Prevalence could not be fitted to given target.")
            }

            # Store the coordinates that were picked for the ellipsoid
            k_sp_params[[paste0("sp_", sp_i)]] <- list(
                means = suitability$details$means,
                sds = suitability$details$sds
            )

            # Determine presence-absence based on target species prevalence
            if (is.null(sites)){
                # use virtualspecies default function
                results_PA <- suppressMessages(convertToPA(
                    suitability,
                    PA.method = "probability",
                    prob.method = "logistic",
                    alpha = alpha,
                    species.prevalence = target, 
                    plot = FALSE))
                k_sp_suit[[paste0("sp_", sp_i)]] <- results_PA$suitab.raster
                k_sp_prob[[paste0("sp_", sp_i)]] <- results_PA$probability.of.occurrence
                k_sp_conv[[paste0("sp_", sp_i)]] <- results_PA$PA.conversion
                k_sp_obs[[paste0("sp_", sp_i)]] <- results_PA$pa.raster

            } else {
                prob_rast <- 1 / (1 + exp((suitability$suitab.raster - beta) / alpha))
                pa_rast <- terra::app(prob_rast, function(p) rbinom(length(p), 1, p))
                
                k_sp_suit[[paste0("sp_", sp_i)]] <- suitability$suitab.raster
                k_sp_prob[[paste0("sp_", sp_i)]] <- prob_rast
                k_sp_conv[[paste0("sp_", sp_i)]] <- list(alpha = alpha, beta = beta)
                k_sp_obs[[paste0("sp_", sp_i)]]  <- pa_rast
            }
        }
        k_simulations$sp_params[[k]] <- k_sp_params
        k_simulations$sp_suitability[[k]] <- k_sp_suit
        k_simulations$sp_probability[[k]] <- k_sp_prob
        k_simulations$sp_pa_params[[k]] <- k_sp_conv
        k_simulations$sp_observations[[k]] <- k_sp_obs
    }

    cli_status_clear(status_msg)
    cli_alert_success("Simulated species are ready!")

    return(k_simulations)
}
