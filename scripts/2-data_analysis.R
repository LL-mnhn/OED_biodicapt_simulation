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
library(patchwork)
library(bestNormalize)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_analysis.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters
set.seed(42) # reproducibility 

if (authorise_overwrite(EXPLORE_PATH)) {
    unlink(EXPLORE_PATH, recursive = TRUE)
    dir.create(EXPLORE_PATH)
} else {
    stop("User refused to overwrite previous analysis.")
}

FULL_PATHS <- c(
    BIODICAPT_OBS_FULL, ENI500_OBS_FULL, STOC_OBS_FULL,
    FEATURES_RAW_FULL, FEATURES_RES_FULL 
)
TRF_PATHS <- c(
    BIODICAPT_OBS_TRF, ENI500_OBS_TRF, STOC_OBS_TRF,
    FEATURES_RAW_TRF, FEATURES_RES_TRF
)
 

##### Local functions ##### ---------------------------------------------------
fit_transforms <- function(train) {
    train <- fix_names(train)
    vals <- as.data.frame(train)   # SpatVector -> data.frame
    setNames(lapply(X_VARIABLES, function(v) {
        bestNormalize(
            vals[[v]],
            allow_orderNorm = FALSE,  # rank-based, can't be applied sensibly to rasters
            out_of_sample   = FALSE,  # no repeated CV, much faster
            quiet           = TRUE
        )
    }), X_VARIABLES)
}

transform_dataset <- function(dataset, bn_list = X_BN) {
    dataset <- fix_names(dataset)

    for (variable in X_VARIABLES) {
        bn <- bn_list[[variable]]

        if (inherits(dataset, "SpatRaster")) {
            layer <- terra::app(dataset[[variable]],
                                function(x) predict(bn, newdata = x))
            names(layer) <- variable
            dataset[[variable]] <- layer
        } else if (inherits(dataset, "SpatVector")) {
            dataset[[variable]] <- predict(bn, newdata = values(dataset)[[variable]])
        } else {
            stop("Unknown type of dataset.")
        }
    }
    dataset
}

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

        # Transformed data
        explore_dataset(
            df = fix_names(vect(BIODICAPT_OBS_TRF)), 
            x_cols = X_VARIABLES, 
            save_folder = EXPLORE_PATH, 
            save_name = "transformed_biodicapt",
            sp_cols = NULL)
        explore_dataset(
            df = fix_names(vect(ENI500_OBS_TRF)), 
            x_cols = X_VARIABLES, 
            save_folder = EXPLORE_PATH, 
            save_name = "transformed_eni500",
            sp_cols = NULL)
        explore_dataset(
            df = fix_names(vect(STOC_OBS_TRF)), 
            x_cols = X_VARIABLES, 
            save_folder = EXPLORE_PATH, 
            save_name = "transformed_stoc",
            sp_cols = fix_names(NAMES_SPECIES), 
            top = 10)
        
        # creating a single df with all datasets together
        combined_trf_datasets <- rbind(
            fix_names(vect(BIODICAPT_OBS_TRF)), 
            fix_names(vect(ENI500_OBS_TRF))
        )
        explore_dataset(
            df = combined_trf_datasets, 
            x_cols = X_VARIABLES, 
            save_folder = EXPLORE_PATH, 
            save_name = "transformed_biodicapt&eni500")
    } else {
        cli_alert_info("Skipping datasets exploration plots.")
    }
}

histogram_covariables <- function(df, n_threshold = 5000) {
    long_df <- df |>
        pivot_longer(
            cols = everything(),
            names_to = "variable",
            values_to = "value"
        )

    # Build the annotation label for one variable
    make_label <- function(x) {
        if (length(x) <= n_threshold) {
            # Small sample: Shapiro-Wilk (shapiro.test supports 3 <= n <= 5000)
            sw <- tryCatch(shapiro.test(x), error = function(e) NULL)
            if (is.null(sw)) return("Shapiro-Wilk: not computable")

            p_txt <- if (sw$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", sw$p.value)
            sprintf("Shapiro-Wilk W = %.3f\n%s", unname(sw$statistic), p_txt)
        } else {
            # Large sample: skewness and excess kurtosis
            m       <- mean(x)
            s2      <- mean((x - m)^2)  # population variance
            skew    <- mean((x - m)^3) / s2^1.5
            ex_kurt <- mean((x - m)^4) / s2^2 - 3
            sprintf("Skewness = %.2f\nExcess kurtosis = %.2f", skew, ex_kurt)
        }
    }

    stats_df <- long_df |>
        filter(!is.na(value)) |>
        group_by(variable) |>
        summarise(label = make_label(value), .groups = "drop")

    plot <- ggplot(data = long_df, aes(x = value, y = after_stat(density))) +
        geom_histogram(bins = 25) +
        geom_text(
            data = stats_df,
            aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.3,
            colour = "black", size = 3
        ) +
        facet_wrap(~ variable, scales = "free", nrow = 3) +
        theme_minimal()

    plot <- my_custom_ggplot_theme(plot, with_palette = FALSE) +
        labs(caption = paste0(
            "Total area under each histogram sums up to 1\n",
            "Shapiro-Wilk test shown when n <= ", n_threshold, " (H0: normality); ",
            "otherwise skewness and excess kurtosis (normal = 0 for both)"
        ))

    return(plot)
}

qqplot_covariables <- function(df, n_points = 10000) {
    probs <- ppoints(n_points)

    qq_df <- df |>
        pivot_longer(
        cols = everything(),
        names_to = "variable",
        values_to = "value"
        ) |>
        filter(!is.na(value)) |>
        group_by(variable) |>
        # standardize so the reference line is y = x
        mutate(z = (value - mean(value)) / sd(value)) |>
        reframe(
        theoretical = qnorm(probs),
        sample      = quantile(z, probs = probs, names = FALSE)
        )

    plot <- ggplot(qq_df, aes(x = theoretical, y = sample)) +
        geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.6) +
        geom_point(size = 0.6, alpha = 0.7) +
        facet_wrap(~ variable, scales = "free_y", nrow = 3) +
        labs(
        x = "Theoretical normal quantiles",
        y = "Standardized sample quantiles"
        ) +
        theme_minimal()

    plot <- my_custom_ggplot_theme(plot, with_palette = FALSE) +
        labs(caption = paste0(
        "Normal QQ plot on standardized variables (", n_points, " quantiles per variable)\n",
        "Points on the red line = normal; S-shapes = heavy/light tails; curves = skew"
        ))

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


##### Transform datasets ##### ------------------------------------------------
# Transformations are fitted on BIODICAPT + 500 ENI (since it is our training data)
cli_alert_info("Fitting transformations to normalise datasets...")
combined_training_data <- vect(list(
    fix_names(vect(BIODICAPT_OBS_FULL)),
    fix_names(vect(ENI500_OBS_FULL))
))
combined_training_data <- combined_training_data[X_VARIABLES]
X_BN <- fit_transforms(combined_training_data)
saveRDS(X_BN, file.path(EXPLORE_PATH, "bestnormalize_fits.rds"))

# Get new units
X_TRANSFORMED_UNITS <- X_UNITS
for (i in seq_along(X_VARIABLES)) {
    cls <- class(X_BN[[i]]$chosen_transform)[1]
    if (cls == "double_reverse_log") next
    u <- X_UNITS[[i]]
    X_TRANSFORMED_UNITS[[i]] <-
        if (is.expression(u)) {
            as.expression(as.call(list(as.name(cls), u[[1]])))
        } else {
            paste0(cls, "(", u, ")")
        }
}
saveRDS(X_TRANSFORMED_UNITS, file.path(EXPLORE_PATH, "tranformed_units.rds"))


cli_alert_info("Application of transformation on environmental variables...")
for (p in seq_along(FULL_PATHS)) {
    if (endsWith(FULL_PATHS[p], '.tif')) {
        data <- rast(FULL_PATHS[p])
        data <- transform_dataset(data)
        writeRaster(data, TRF_PATHS[p], overwrite = TRUE)
    } else if (endsWith(FULL_PATHS[p], '.gpkg')) {
        data <- vect(FULL_PATHS[p])
        data <- transform_dataset(data)
        writeVector(data, TRF_PATHS[p], overwrite = TRUE)
    } else {
        stop("Unknown type of dataset.")
    }
}
cli_alert_info("Transformed datasets are ready!")


##### Analysis of env variables ##### -----------------------------------------
cli_alert_info("Analysis of env variables...")

# The HMSC model will be trained on obs from BIODICAPT and ENI500
# So we use them to study the distribution of our environmental variables
combined_paths <- list(
    base = c(BIODICAPT_OBS_FULL, ENI500_OBS_FULL),
    transformed = c(BIODICAPT_OBS_TRF, ENI500_OBS_TRF)
)
for (comb in names(combined_paths)) {
    sample_pool_df <- as.data.frame(
        vect(list(
            fix_names(vect(combined_paths[[comb]][1])),
            fix_names(vect(combined_paths[[comb]][2]))
        ))
    )
    sample_pool_df <- sample_pool_df[X_VARIABLES]

    ### Histograms
    hist_grid_plot <- histogram_covariables(sample_pool_df)
    print(hist_grid_plot)
    standardised_ggplot_save(
        hist_grid_plot, 
        file.path(EXPLORE_PATH, paste0(comb, "hist_on_biodicap&eni500.pdf")))

    qqplot_grid_plot <- qqplot_covariables(sample_pool_df)
    print(qqplot_grid_plot)
    standardised_ggplot_save(
        qqplot_grid_plot, 
        file.path(EXPLORE_PATH, paste0(comb, "qqplot_on_biodicap&eni500.pdf")))
    
}
cli_alert_success("Plots are ready!\n\n")


##### Full data summary ##### -------------------------------------------------
. <- save_exploration()


##### Load datasets ##### -----------------------------------------------------
cli_alert_info("Showing raster layers...")
# environmental data for virtual species
t_envs_raster <- get_master_t_rasters(mode = "raw", keep.lyr = X_VARIABLES)
t_envs_units <- readRDS(file.path(EXPLORE_PATH, "tranformed_units.rds"))

# Save all layers to a figure
plots <- lapply(seq_along(names(t_envs_raster)), function(i) {
    n <- names(t_envs_raster)[i]

    p <- ggplot() +
        geom_spatraster(data = t_envs_raster[[i]]) +
        ggtitle(n) +
        scale_fill_viridis_c(
        option = "magma",
        name = t_envs_units[[i]]
        )

    my_custom_ggplot_theme(p, with_palette = FALSE) +
        theme(plot.title = element_text(size = 9))
})
standardised_ggplot_save(
    wrap_plots(plots), 
    file.path(EXPLORE_PATH, "env_raster_layers_vtransform.pdf"), 
    .width = 36, .height = 12)


##### Computing PCA (for `virtualspecies`) ##### ------------------------------
### Compute PCA on raster
t_env_df <- as.data.frame(t_envs_raster)
result_pca <- make_pca_and_plots(t_env_df)
margin_pca <- ggMarginal(result_pca$plot_ind, type = "density")
print(margin_pca)
print(result_pca$plot_var)
standardised_ggplot_save(
    margin_pca, 
    file.path(EXPLORE_PATH, "pca_on_t_rasters.pdf"))
saveRDS(result_pca$pca, file = file.path(EXPLORE_PATH, "pca_data.rds"))

### Visualise real-world sigthings on PCA
cli_alert_info("Analysis of sightings...")
# species + environmental dataset (a sort of biased "ground truth")
t_stoc_df <- fix_names(vect(STOC_OBS_TRF))

# Plot ellipsoids
ellipses_pca_plot <- ellipses_of_all_covariables(result_pca, t_stoc_df)
print(ellipses_pca_plot)
standardised_ggplot_save(
    ellipses_pca_plot, file.path(EXPLORE_PATH, "pca_ellipses_stoc.pdf"))
project_on_existing_pca(
    result_pca$pca, result_pca$plot_ind, t_stoc_df, group = NULL)

### Plot all densities
densities_plot <- histogram_covariables_plus_densities(
    hist_grid_plot, t_stoc_df, x_cols = names(t_env_df))
print(densities_plot)
standardised_ggplot_save(densities_plot, 
    file.path(EXPLORE_PATH, "hist_rasters_and_stoc.pdf"))


### Plot observed rank-occurence curve 
# (proxy to rank-abundance for presence/absence data)
rank_occurence_plot <- rank_occurence_curves(t_stoc_df, fix_names(NAMES_SPECIES))
standardised_ggplot_save(
    rank_occurence_plot$roc_plot, 
    file.path(EXPLORE_PATH, "rank-occurence_stoc.pdf"))
write_csv(
    rank_occurence_plot$occ_table, 
    file.path(EXPLORE_PATH, "rank-occurence_stoc.csv"))

cli_alert_success("File '2-data_analysis.R' finished running!\n\n")