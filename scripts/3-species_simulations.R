# This script:
#   - uses raster with environmental data in multiple layers
#   - uses shapefiles with environmental and location data
#   - generates virtual species distriutions from credible parameters

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(terra)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggExtra)
library(tidyterra)
library(virtualspecies)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_analysis.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_simul.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

REMOVE_LAYERS_FROM_MASTER <- c("CLC", "distance to NODATA")
CORRELATED_LAYERS_IN_MASTER <- list(  # keep first name for each vector in list
    c("CHELSA_tas", "Elevation")
) 
SHOW_SP <- 3 # picks this number of species in simulations to show their maps


if (authorise_overwrite(SIMULATE_PATH)) {
    unlink(SIMULATE_PATH, recursive = TRUE)
    dir.create(SIMULATE_PATH)
} else {
    stop("User refused to overwrite previous simulations.")
}


##### Local functions ##### ---------------------------------------------------
plot_virtualspecies_rasters <- function() {
    idx_sp <- sample(seq(length(NAMES_SPECIES)), SHOW_SP)
    layer_names <- c("sp_suitability", "sp_probability", "sp_observations")

    for (sp in idx_sp) {
        # the 3 rasters for this species
        rasters <- lapply(
            layer_names, 
            function(l) simulations[[l]][[1]][[paste0("sp_", sp)]])

        plots <- lapply(seq_along(rasters), function(i) {
            p <- ggplot() +
            geom_spatraster(data = rasters[[i]]) +
            ggtitle(layer_names[i]) +
            scale_fill_stepsn(
                colours = GG_TERRAIN_PALETTE,
                n.breaks = 9,
                limits = c(0, 1)
            )
            my_custom_ggplot_theme(p, with_palette = FALSE)
        })

        combined <- wrap_plots(plots, nrow = 1, guides = "collect")
        print(combined)
        
        standardised_ggplot_save(
            combined, 
            file.path(SIMULATE_PATH, paste0("simulation_sp", sp,"_k1.pdf")))
    }
}




##### Load datasets ##### -----------------------------------------------------
cli_alert_info("Loading datasets...")
# environmental data for virtual species (downsampled version)
envs_raster_subres <- get_master_raster(
    mode = "res",
    rm.lyr = REMOVE_LAYERS_FROM_MASTER, 
    rm.cor.lyr = CORRELATED_LAYERS_IN_MASTER)

# species + environmental dataset (a sort of biased "ground truth")
stoc_df <- fix_names(vect(STOC_OBS_FULL))

# Other data precomputed in data-analysis
env_pca <- readRDS(file.path(EXPLORE_PATH, "pca_data.rds"))
stoc_occurences <- read_csv(
    file.path(EXPLORE_PATH, "stoc_occurences_in_full_dataset.csv"),
    show_col_types = FALSE)

cli_alert_success("Datasets loaded!\n\n")


##### Simulations ##### -------------------------------------------------------
cli_alert_info("Simulation of species distributions...")

### 1. Use base virtualspecies::generateSpFromPCA function 
#       - This function does all the work for us, we just need to specify
#           species.prevalence to mimic our observed species prevalence
simulations <- simulate_from_PCA(
    raster_stack = envs_raster_subres, 
    precomputed_pca = env_pca, 
    n_sp = length(NAMES_SPECIES), 
    prevalences = stoc_occurences$frequency)

# Saving dont forget to wrap any SpatRaster objects inside the list)
saveRDS(
    wrap_simulations(simulations), 
    file = file.path(SIMULATE_PATH, "simulations.rds"))
plot_virtualspecies_rasters()

# compare occurences with stoc_df
stoc_simulations <- list()
for (k in 1:K_FOLDS) {
    local_simul_obs_df <- stoc_df
    
    for (sp_i in 1:length(NAMES_SPECIES)) {
        # Use CustomSampleOccurences for fixed locations
        local_simul_obs_df[[paste0("sp_", sp_i)]] <- customSampleOccurrences(
            simulations$sp_observations[[k]][[paste0("sp_", sp_i)]],
            stoc_df,
            detection.probability = 0.95
        )

    }

    stoc_simulations[[k]] <- local_simul_obs_df[paste0("sp_", seq(1:length(NAMES_SPECIES)))]
}

sim_rank_occurence_plot <- rank_occurence_curves(
    stoc_df,
    fix_names(NAMES_SPECIES),
    stoc_simulations, 
    paste0("sp_", seq(1:length(NAMES_SPECIES))))
print(sim_rank_occurence_plot$roc_plot)
standardised_ggplot_save(
    sim_rank_occurence_plot$roc_plot, 
    file.path(SIMULATE_PATH, "sim_rank-occurence_stoc.pdf"))
