# This script:
#   - Analyses the results of the I-optimisation runs

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_simul.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_model.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters

EFFECTS <- readRDS(file.path(RESULTS_PATH, "effects.rds"))
REFERENCE_EFFECTS <- EFFECTS[1][[1]]

if (authorise_overwrite(ANALYSES_PATH)) {
    unlink(ANALYSES_PATH, recursive = TRUE)
    dir.create(ANALYSES_PATH)
} else {
    if (file.exists(file.path(ANALYSES_PATH, "all_scores.csv"))) {
        cli_alert_warning("Auto-authorisation: re-using existing 'all_scores.csv'.")
    } else {
        stop("User refused to overwrite previous analyses of model's runs.")
    }
    
}

ORDER_STRATS <- c(
    "business-as-usual", "gap-filling", "simplified-I-optimality", 
    "2-part-simplified-I-optimality", "5-part-simplified-I-optimality",
    "10-part-simplified-I-optimality")
EXPLORE_COLUMNS <- c(
    "chelsa_hurs", "chelsa_pr", "chelsa_tas", "ndvi", "light_pollution", 
    "distance.to.arable.land", "distance.to.permanent.crops",
    "distance.to.pastures", "distance.to.wetlants",
    "distance.to.heterogeneous.agricultural.areas",
    "distance.to.forest.and.semi.natural.areas", 
    "distance.to.water.bodies")   


##### Local functions ##### ---------------------------------------------------
fetch_all_scores <- function() {
    cli_alert_info("Loading scores...")

    # load datasets
    datasets <- unwrap_simulations(readRDS(
        file.path(RESULTS_PATH, "datasets.rds")))
    all_scores <- tibble()
    run_paths <- list.dirs(RESULTS_PATH)[-1]

    for (i in 1:length(run_paths)) {
        cli_alert_info("[Treating run {i}/{length(run_paths)}]")
        # load parameters 
        local_parameters <- readRDS(
            file.path(run_paths[i], "local_parameters.rds"))
        
        # load model chains
        chains <- readRDS(file.path(run_paths[i], "chains.rds"))

        # select dataset
        k_dataset <- datasets[[local_parameters$k]]

        for (j in 1:length(k_dataset)) {
            cli_alert_info("    - computing scores on dataset {j}/{length(k_dataset)} ...")

            if (names(k_dataset)[j] == "new_pool") {
                next
            }

            if (nrow(k_dataset[[j]]) > 500) {
                k_set <- k_dataset[[j]] |> slice_sample(n = 500)
            } else {
                k_set <- k_dataset[[j]]
            }

            scores <- evaluate_hmsc_performances(
                hM = chains, data_set = k_set, 
                x_cols = X_VARIABLES, sp_cols = NAME_SP_SIMUL)   

            scores <- scores |>
                lapply(function(x) tibble(sp = names(x), score = unname(x))) |>
                bind_rows(.id = "name")
            
            for (p in 1:length(local_parameters)){
                scores[names(local_parameters)[p]] <- local_parameters[[names(local_parameters)[p]]]
            }
            scores["subset"] <- names(k_dataset)[j]

            all_scores <- rbind(all_scores, scores)
        } 
    }

    cli_alert_success("Scores loaded in data.frame!\n\n")
    return(all_scores)
}

fetch_all_new_samples_id <- function() {
    run_paths <- list.dirs(RESULTS_PATH)[-1]
    all_new_samples_df <- data.frame(
        sample_id = character(),
        r_effect = character(),
        new_samples_size = numeric(),
        strategies = character(),
        k_folds = integer()
    )

    for (i in 1:length(run_paths)) {
        # load parameters 
        local_parameters <- readRDS(
            file.path(run_paths[i], "local_parameters.rds"))
        
        if (local_parameters$NEW_SAMPLE_SIZES == 0) {
            next
        }
        
        # load model chains
        new_samples <- readRDS(file.path(run_paths[i], "ID_new_samples.rds"))
        new_samples_df <- data.frame(
            sample_id = new_samples,
            r_effect = local_parameters$R_EFFECTS,
            new_samples_size = local_parameters$NEW_SAMPLE_SIZES,
            strategies = local_parameters$STRATEGIES,
            k_folds = local_parameters$k
        )

        # add to dataframe
        all_new_samples_df <- rbind(all_new_samples_df, new_samples_df)
    }

    return(all_new_samples_df)
}


##### Fetch scores ##### ------------------------------------------------------
# Compute models scores
if (file.exists(file.path(ANALYSES_PATH, "all_scores.csv"))) {
    scores_df <- read_csv(
        file.path(ANALYSES_PATH, "all_scores.csv"),
        show_col_types = FALSE)
} else {
    scores_df <- fetch_all_scores()
    write_csv(scores_df, file.path(ANALYSES_PATH, "all_scores.csv"))
}


# Add MSE metric
scores_df <- scores_df |>
    bind_rows(
        scores_df |>
            filter(name == "RMSE") |>
            mutate(name = "MSE", score = score^2)
    )
scores_df <- scores_df |>
    mutate(name = if_else(name == "SD", "SD of predictions", name))


# Divide reference and new values
ref_df <- scores_df |>
    semi_join(as_tibble(REFERENCE_EFFECTS), by = names(REFERENCE_EFFECTS))
compare_df <- scores_df |>
    anti_join(as_tibble(REFERENCE_EFFECTS), by = names(REFERENCE_EFFECTS))


##### Analyses ##### ----------------------------------------------------------
### 1st objective: which strategy gives us the best model predictions
### Quality of prediction
compare_plot(
        ref_scores = ref_df,
        compare_scores = compare_df,
        metric = "MSE",
        x_var = "STRATEGIES",
        x_order = ORDER_STRATS,
        fill_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset",
        panel_order = c("train", "new_pool", "test"),
        group_species = TRUE,
        save_to = file.path(ANALYSES_PATH, "MSE_deltas.pdf"))

# Compute improvement in MSE
improvement_df <- count_improving(
    ref_df, compare_df,
    metric = "RMSE",
    threshold = 0.01,
    group_vars = c("STRATEGIES", "NEW_SAMPLE_SIZES"),
    panel_var = "subset")
improvement_df["name"] <- "# of species with\nincrease in ΔMSE >1%"
improvement_df <- improvement_df |> rename(diff_score = n_sp_above)
compare_plot(
        diffs_scores = improvement_df,
        metric = unique(improvement_df["name"])[[1]],
        x_var = "STRATEGIES",
        x_order = ORDER_STRATS,
        fill_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset",
        panel_order = c("train", "new_pool", "test"),
        group_species = TRUE,
        save_to = file.path(ANALYSES_PATH, "MSE_delta_n_sp.pdf"))


### Uncertainty of prediction
compare_plot(
        ref_scores = ref_df,
        compare_scores = compare_df,
        metric = "SD of predictions",
        x_var = "STRATEGIES",
        x_order = ORDER_STRATS,
        fill_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset",
        panel_order = c("train", "new_pool", "test"),
        group_species = TRUE,
        save_to = file.path(ANALYSES_PATH, "SD_deltas.pdf"))


##### MAPS ###### -------------------------------------------------------------
### Map locations sampled in new_pool
### TODO: turn it into a nice function
all_datasets <- unwrap_simulations(
    readRDS(file = file.path(RESULTS_PATH, "datasets.rds")))
# new pool of points
available_points <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map(), 
        df = as.data.frame(all_datasets[[1]]$new_pool), 
        lon_c = "lon",
        lat_c = "lat",
        col = "black", shape = 20, size = 1)
available_points <- available_points + 
    labs(caption = paste0(
            "Red squares: training points (biodicapt); ",
            "Black circles: new points available (500 ENI)."))
print(available_points)

# fetch points that were selected from new_pool during each run
selected_new_df <- fetch_all_new_samples_id()

# add lon and lat info
selected_new_df <- selected_new_df |>
  left_join(
    as.data.frame(all_datasets[[1]]$new_pool)[c("lon", "lat", "row_id")], 
    by = c("sample_id" = "row_id"))

# Order each category
selected_new_df <- selected_new_df |>
    mutate(strategies = factor(strategies, levels = ORDER_STRATS))

# Frequency across k_folds
new_samples_freq <- selected_new_df |>
    group_by(strategies, r_effect, new_samples_size, sample_id, lon, lat) |>
    summarise(n_folds = n_distinct(k_folds), .groups = "drop")

plot_grid <- ggplot_quantitative_df_on_background_map(
    background_map = available_points,
    df = new_samples_freq |> filter(r_effect == "none", new_samples_size == 50),
    cmap = viridisLite::plasma(256),
    lon_c = "lon", lat_c = "lat",
    column = "n_folds", unit = "# of selections",
    facet_formula = ~ strategies
)
print(plot_grid)
standardised_ggplot_save(
    plot_grid, 
    file.path(ANALYSES_PATH, "locations_freq_strategies_new-50.pdf"))

plot_grid <- ggplot_quantitative_df_on_background_map(
    background_map = available_points,
    df = new_samples_freq |> filter(r_effect == "none", strategies == "simplified-I-optimality"),
    cmap = viridisLite::plasma(256),
    lon_c = "lon", lat_c = "lat",
    column = "n_folds", unit = "# of selections",
    facet_formula = ~ new_samples_size
)
print(plot_grid)
standardised_ggplot_save(
    plot_grid, 
    file.path(ANALYSES_PATH, "locations_freq_new-sample-sizes_strategy-SIO.pdf"))



### Grid of histograms/densities
### TODO: turn it into a nice function
# fetch points that were selected from new_pool during each run
selected_new_features <- fetch_all_new_samples_id()

# add features
selected_new_features <- selected_new_features |>
  left_join(
    as.data.frame(all_datasets[[1]]$new_pool)[c(EXPLORE_COLUMNS, "row_id")], 
    by = c("sample_id" = "row_id"))


# Reference data in long format, labelled "reference"
to_long <- function(df, group_label = NULL, group_col = "k_folds") {
    out <- df |>
        select(any_of(group_col), all_of(EXPLORE_COLUMNS)) |>
        pivot_longer(
        cols = all_of(EXPLORE_COLUMNS),
        names_to = "variable",
        values_to = "value"
        ) |>
        mutate(variable = factor(variable, levels = EXPLORE_COLUMNS))

    if (!is.null(group_label)) out[[group_col]] <- group_label
    out[[group_col]] <- as.character(out[[group_col]])
    out
}

plot_with_references <- function(
        selected, ref_left, ref_right,
        left_label = "reference",
        right_label = "reference_bis") {
    folds_long <- to_long(selected)
    fold_levels <- as.character(sort(unique(as.numeric(folds_long$k_folds))))

    all_levels <- c(left_label, fold_levels, right_label)

    plot_all <- bind_rows(
        to_long(ref_left,  left_label),
        folds_long,
        to_long(ref_right, right_label)
    ) |>
        mutate(k_folds = factor(k_folds, levels = all_levels))

    # Colors: gray shades for references, default palette for folds
    fill_values <- c(
        setNames("gray60", left_label),
        setNames(scales::hue_pal()(length(fold_levels)), fold_levels),
        setNames("gray30", right_label)
    )

    # Separator positions: between left ref and first fold, and between last fold and right ref
    sep_left  <- 1.5
    sep_right <- length(all_levels) - 0.5

    ggplot(plot_all, aes(x = k_folds, y = value, fill = k_folds)) +
        geom_violin(alpha = 0.4, adjust = .5) +
        geom_boxplot(width = 0.15, outlier.size = 0) +
        geom_vline(xintercept = c(sep_left, sep_right),
                colour = "gray40", linetype = "dashed") +
        scale_fill_manual(values = fill_values) +
        facet_wrap(~ variable, scales = "free_y", ncol = 3) +
        labs(x = NULL, y = NULL) +
        theme_bw() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )
}

# --- Usage -----------------------------------------------------------------

selected <- selected_new_features |>
  filter(
    r_effect == "none",
    new_samples_size == 50,
    strategies == "simplified-I-optimality"
  )

p <- plot_with_references(
    selected  = selected,
    ref_left  = as.data.frame(envs_raster_subres)[EXPLORE_COLUMNS],
    ref_right = as.data.frame(all_datasets[[1]]$new_pool)[EXPLORE_COLUMNS],
    left_label  = "France",
    right_label = "500 ENI"
)

standardised_ggplot_save(
    my_custom_ggplot_theme(p, with_palette = FALSE)  +
        theme(
            aspect.ratio = 0.5,
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        ),
    file.path(ANALYSES_PATH, "distributions_in_selected_locations.pdf"),
    .width = 36, .height = 12
)
