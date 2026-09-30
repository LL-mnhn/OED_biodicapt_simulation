# This script:
#   - uses pre-processed datasets
#   - analyses the datasets distributions

##### Libraries ##### ---------------------------------------------------------
library(sf)
library(cli)
library(terra)
library(tidyr)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggExtra)
library(stringr)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_analysis.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

EXPLORE_COLUMNS <- c(
    "clc", "chelsa_hurs", "chelsa_pr", "chelsa_tas", "ndvi", "elevation",
    "light_pollution", "distance.to.arable.land", "distance.to.permanent.crops",
    "distance.to.pastures", "distance.to.wetlants",
    "distance.to.heterogeneous.agricultural.areas",
    "distance.to.forest.and.semi.natural.areas", 
    "distance.to.water.bodies")            
REMOVE_LAYERS_FROM_MASTER <- c("CLC", "distance to NODATA")
CORRELATED_LAYERS_IN_MASTER <- list(  # keep first name for each vector in list
    c("CHELSA_tas", "Elevation")
) 

if (authorise_overwrite(EXPLORE_PATH)) {
    unlink(EXPLORE_PATH, recursive = TRUE)
    dir.create(EXPLORE_PATH)
} else {
    stop("User refused to overwrite previous analysis.")
}

##### Local functions ##### ---------------------------------------------------
save_exploration <- function() {
    make_plots <- typeline(
        prompt = paste0(paste0(
            "Show exploration plots of data, save figures & occurrences?",
            "(overwrites previous files in 'outputs/data_exploration')? [Y/n]: ")))
    cleaned_answer <- tolower(trimws(make_plots))

    if (cleaned_answer %in% c("y", "yes")) {
        cli_alert_warning(paste0(
            "All figures are automatically saved as .pdf inside '", 
            EXPLORE_PATH, "' ."))
        cli_alert_info("Showing plots to the user.")

        explore_dataset(
            df = fix_names(vect(BIODICAPT_OBS_FULL)), 
            x_cols = EXPLORE_COLUMNS, 
            save_folder = EXPLORE_PATH, 
            save_name = "biodicapt",
            sp_cols = NULL)
        explore_dataset(
            df = fix_names(vect(ENI500_OBS_FULL)), 
            x_cols = EXPLORE_COLUMNS, 
            save_folder = EXPLORE_PATH, 
            save_name = "eni500",
            sp_cols = NULL)
        explore_dataset(
            df = fix_names(vect(STOC_OBS_FULL)), 
            x_cols = EXPLORE_COLUMNS, 
            save_folder = EXPLORE_PATH, 
            save_name = "stoc",
            sp_cols = fix_names(NAMES_SPECIES), 
            top = 10)
        
        # creating a single df with all datasets together
        combined_datasets <- rbind(
            fix_names(vect(BIODICAPT_OBS_FULL)), 
            fix_names(vect(ENI500_OBS_FULL)), 
            fix_names(vect(STOC_OBS_FULL))
        )
        explore_dataset(
            df = combined_datasets, 
            x_cols = EXPLORE_COLUMNS, 
            save_folder = EXPLORE_PATH, 
            save_name = "concat")
    } else {
        cli_alert_info("Skipping datasets exploration plots.")
    }
}

histogram_covariables <- function(df) {
    long_df <- df |>
        pivot_longer(
            cols = everything(), 
            names_to = "variable", 
            values_to = "value") 
    plot <- ggplot(data = long_df, aes(x = value, y = after_stat(density))) +
        geom_histogram(bins = 25) +
        facet_wrap(~ variable, scales = "free", nrow = 3) +
        theme_minimal()
    plot <- my_custom_ggplot_theme(plot, with_palette = FALSE) +
        labs(caption = "Total area under each histogram sums up to 1")
    return(plot)
}

histogram_covariables_plus_densities <- function(hist_plot, df, x_cols) {
    # Build one long df combining all species' presence rows
    long_species_df <- map_dfr(fix_names(NAMES_SPECIES), function(sp) {
        as.data.frame(df) |>
            filter(.data[[sp]] == 1) |>
            select(all_of(x_cols)) |>
            pivot_longer(cols = everything(),
                        names_to = "variable", values_to = "value") |>
            mutate(species = sp)
    })

    hist_plot <- hist_plot +
        geom_density(
            data = long_species_df,
            aes(x = value, y = after_stat(density), color = species),
            inherit.aes = FALSE,
            linewidth = 0.4,
            alpha = 0.7
        ) +
        scale_color_manual(values = rainbow(40))
    return(hist_plot)
}

ellipses_of_all_covariables <- function(old_pca, df) {
    ellipses_df <- data.frame()
    for (sp in fix_names(NAMES_SPECIES)) {
        local_results <- suppressMessages(project_on_existing_pca(
            old_pca$pca, old_pca$plot_ind, df, sp))
        local_ellipse <- compute_ellipse(local_results$sp_coords)
        local_ellipse$species <- sp
        ellipses_df <- rbind(ellipses_df, local_ellipse)
    }
    extended_pca_plot <- ggplot() +
        scale_color_manual(values = rainbow(40)) +
        scale_fill_manual(values = rainbow(40)) +
        geom_polygon(
            data = ellipses_df,
            aes(x = Dim.1, y = Dim.2, group = species, color = species, fill = species),
            alpha = 0.15
        )
    return( my_custom_ggplot_theme(extended_pca_plot, with_palette = FALSE))
}


##### Load datasets ##### -----------------------------------------------------
cli_alert_info("Loading datasets...")
# environmental data for virtual species
envs_raster <- get_master_raster(
    mode = "raw", 
    rm.lyr = REMOVE_LAYERS_FROM_MASTER, 
    rm.cor.lyr = CORRELATED_LAYERS_IN_MASTER)

# species + environmental dataset (a sort of biased "ground truth")
stoc_df <- fix_names(vect(STOC_OBS_FULL))
cli_alert_success("Dataset loaded!\n\n")


##### Full data summary ##### -------------------------------------------------
. <- save_exploration()


##### Analysis of env variables ##### -----------------------------------------
cli_alert_info("Analysis of env variables...")

env_df <- as.data.frame(envs_raster)

### PCA
result_pca <- make_pca_and_plots(env_df)
margin_pca <- ggMarginal(result_pca$plot_ind, type = "density")
print(margin_pca)
print(result_pca$plot_var)
standardised_ggplot_save(
    margin_pca, 
    file.path(EXPLORE_PATH, "pca_on_rasters.pdf"))
saveRDS(result_pca$pca, file = file.path(EXPLORE_PATH, "pca_data.rds"))

### Histograms
hist_grid_plot <- histogram_covariables(env_df)
print(hist_grid_plot)
standardised_ggplot_save(
    hist_grid_plot, 
    file.path(EXPLORE_PATH, "hist_on_rasters.pdf"))
cli_alert_success("Plots are ready!\n\n")


##### Analysis of real-world sightings ##### ----------------------------------
cli_alert_info("Analysis of sightings...")
### Plot all ellipsoids
ellipses_pca_plot <- ellipses_of_all_covariables(result_pca, stoc_df)
print(ellipses_pca_plot)
standardised_ggplot_save(
    ellipses_pca_plot, file.path(EXPLORE_PATH, "pca_ellipses_stoc.pdf"))

### Plot all densities
densities_plot <- histogram_covariables_plus_densities(
    hist_grid_plot, stoc_df, x_cols = names(env_df))
print(densities_plot)
standardised_ggplot_save(densities_plot, 
    file.path(EXPLORE_PATH, "hist_rasters_and_stoc.pdf"))

### Plot observed rank-occurence curve 
# (proxy to rank-abundance for presence/absence data)
rank_occurence_plot <- rank_occurence_curves(stoc_df, fix_names(NAMES_SPECIES))
print(rank_occurence_plot$roc_plot)
standardised_ggplot_save(
    rank_occurence_plot$roc_plot, 
    file.path(EXPLORE_PATH, "rank-occurence_stoc.pdf"))
write_csv(
    rank_occurence_plot$occ_table, 
    file.path(EXPLORE_PATH, "rank-occurence_stoc.csv"))

cli_alert_success("File '2-data_analysis.R' finished running!\n\n")