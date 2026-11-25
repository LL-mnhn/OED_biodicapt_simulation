# This script:
#   - Analyses the results of the different runs

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(terra)
library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_simul.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_model.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_analysis.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data","config","config.R"))) # Global parameters
set.seed(42) # reproducibility 

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

EFFECTS <- readRDS(file.path(RESULTS_PATH, "effects.rds"))
REFERENCE_EFFECTS <- EFFECTS[1][[1]]
ORDER_STRATS <- c(
    "none", "random-sampling", "by-region", "gap-filling", "max-uncertainty", 
    "ds-optimality", "i-optimality")

# FIX Random effect (we are not interested in its effect on predictions)
R_EFFECT_PLOT <- "site"
REFERENCE_EFFECTS$R_EFFECTS <-  R_EFFECT_PLOT

# Plot parameters
SUBSETS_PLOT <- c("train", "test_local", "test_stoc") # c("train", "test_local", "test_stoc")


##### Local functions ##### ---------------------------------------------------
# Turns a list of named vectors (one value per species) into a long tibble
# and attaches the parameters of the current run and the subset name.
format_scores <- function(raw_scores, local_parameters, subset_name) {
    scores <- raw_scores |>
        lapply(function(x) tibble(sp = names(x), score = unname(x))) |>
        bind_rows(.id = "name")
    for (p in 1:length(local_parameters)) {
        scores[names(local_parameters)[p]] <- local_parameters[[names(local_parameters)[p]]]
    }
    scores["subset"] <- subset_name
    return(scores)
}

fetch_all_scores <- function() {
    cli_alert_info("Loading scores...")

    # load datasets
    datasets <- unwrap_simulations(readRDS(
        file.path(RESULTS_PATH, "datasets.rds")))
    all_scores <- tibble()
    run_paths <- list.dirs(RESULTS_PATH)[-1]

    # initialise reference model import
    ref_cache <- list(path = NULL, chains = NULL)

    for (i in 1:length(run_paths)) {
        cli_alert_info("[Treating run {i}/{length(run_paths)}]")
        # load parameters / re_trained model / new pool IDs
        local_parameters <- readRDS(
            file.path(run_paths[i], "local_parameters.rds"))
        chains <- readRDS(file.path(run_paths[i], "chains.rds"))
        id_new_samples <- readRDS(file.path(run_paths[i], "ID_new_samples.rds"))

        # Make train + test_local + test_stoc datasets
        k_dataset <- datasets[[local_parameters$k]]
        new_pool_set <- k_dataset$new_pool[
            k_dataset$new_pool$row_id %in% id_new_samples, ]
        unselected_new_pool_set <- k_dataset$new_pool[
            !(k_dataset$new_pool$row_id %in% id_new_samples), ]

        training_set <- vect(list(k_dataset$train, new_pool_set))
        test_local_set <- unselected_new_pool_set
        test_stoc_set <- k_dataset$test

        local_datasets <- list(
            train = training_set,
            test_local = test_local_set,
            test_stoc = test_stoc_set
        )

        ref_chains <- NULL # reset reference model
        ref_info <- NULL # reset information of reference model
        # if needed (aka new samples added) load reference model
        if (local_parameters$NEW_SAMPLE_SIZES > 0) {
            # Update random effects in reference model
            ref_params <- modifyList(
                REFERENCE_EFFECTS,
                list(R_EFFECTS = local_parameters$R_EFFECTS))
            ref_path <- file.path(
                make_run_path(
                    folder = RESULTS_PATH, 
                    params = ref_params, 
                    k_fold = local_parameters$k),
                "chains.rds")
            if (!identical(ref_cache$path, ref_path)) {
                ref_cache$chains <- readRDS(ref_path)
                ref_cache$path <- ref_path
            }
            ref_chains <- ref_cache$chains

            # Get information of the reference model on the selected new sites
            # (same X scaling as the realised metrics, so they are comparable)
            selected_df <- drop_na(as.data.frame(new_pool_set)[c(X_VARIABLES)])
            ref_info <- get_hmsc_information(
                hM = ref_chains, 
                newdata = selected_df, 
                x_cols = X_VARIABLES,
                scaling_data = as.data.frame(k_dataset$train))
        }        

        # Ratio of posterior SDs (retrained / reference), same X scaling
        ratio_scores <- NULL
        if (!is.null(ref_chains)) {
            ratio_scores <- get_beta_sd_ratio(
                hM_new = chains,
                hM_ref = ref_chains,
                x_cols = X_VARIABLES,
                actionable = ACTIONABLE_VARIABLES,
                scaling_data = as.data.frame(k_dataset$train))
        }

        for (j in 1:length(local_datasets)) {
            cli_alert_info("    - computing scores on dataset {j}/{length(local_datasets)} ...")

            if (nrow(local_datasets[[j]]) == 0) {
                next
            }

            if (nrow(local_datasets[[j]]) > 500) {
                k_set <- local_datasets[[j]] |> slice_sample(n = 500)
            } else {
                k_set <- local_datasets[[j]]
            }

            # Base metrics (RMSE, AUC, Tjur2, SD) and parameter-based metrics
            # (BetaSD, BetaSD_actionable, Ds, I).
            scores <- evaluate_hmsc_performances(
                hM = chains, data_set = k_set, 
                x_cols = X_VARIABLES, sp_cols = NAME_SP_SIMUL,
                actionable = ACTIONABLE_VARIABLES,
                scaling_data = as.data.frame(k_dataset$train))  
            scores_formatted <- format_scores(
                scores, local_parameters, names(local_datasets)[j])
            all_scores <- rbind(all_scores, scores_formatted)
            
            if (length(unique(scores_formatted$sp)) != 40) {
                print(unique(scores_formatted$sp))
                stop(paste("stopped at", i, j))
            }

            if (!is.null(ratio_scores)) {
                all_scores <- rbind(
                    all_scores,
                    format_scores(
                        ratio_scores, local_parameters, names(local_datasets)[j]))
            }
            
            # Expected (approximated) Ds and I from the reference model
            if (!is.null(ref_info)) {
                expected_scores <- get_expected_design_metrics(
                    info = ref_info,
                    actionable = ACTIONABLE_VARIABLES,
                    target = as.data.frame(k_set) |> 
                                drop_na() |> 
                                select(all_of(X_VARIABLES))
                    )
                all_scores <- rbind(
                    all_scores,
                    format_scores(
                        expected_scores, 
                        local_parameters, 
                        names(local_datasets)[j]))
            }            

            # Optional: to explore KL divergence
            # # KL divergence (ref || retrained)
            # # Only on "test_stoc" set because thats the only dataset with the 
            # # same number of sites in reference and retrained models.
            # if ((!is.null(ref_chains)) && (names(local_datasets)[j] == "test_stoc")) {
            #     cli_alert_info("    - computing KL divergence vs reference ...")
            #     kl_scores <- kl_between_models(
            #         hM_ref = ref_chains, 
            #         hM_new = chains, 
            #         df = tidyr::drop_na(as.data.frame(k_set)), 
            #         x_cols = X_VARIABLES,
            #         direction = "new_vs_ref")
            #     all_scores <- rbind(
            #         all_scores, 
            #         format_scores(kl_scores, local_parameters, names(local_datasets)[j]))
            # }
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

show_map_selected_new_pool <- function(
    complete_new_pool, keep_strats = NULL, keep_sizes = NULL, save_to = NULL) {
    # new pool of points
    available_points <- ggplot_categorical_df_on_background_map(
            background_map = ggplot_get_france_base_map(), 
            df = as.data.frame(complete_new_pool), 
            lon_c = "lon",
            lat_c = "lat",
            col = "black", shape = 20, size = 1)
    available_points <- available_points + 
        labs(caption = paste0(
                "Red squares: training points (biodicapt); ",
                "Black circles: new points available (500 ENI)."))
    
    # fetch points that were selected from new_pool during each run
    selected_new_df <- fetch_all_new_samples_id()

    # add lon and lat info
    selected_new_df <- selected_new_df |>
    inner_join(
        as.data.frame(complete_new_pool)[c("lon", "lat", "row_id")], 
        by = c("sample_id" = "row_id"))
    
    # Order each category
    selected_new_df <- selected_new_df |>
        mutate(strategies = factor(strategies, levels = ORDER_STRATS))

    # Frequency across k_folds
    new_samples_freq <- selected_new_df |>
        group_by(strategies, r_effect, new_samples_size, sample_id, lon, lat) |>
        summarise(n_folds = n_distinct(k_folds), .groups = "drop")

    if (!is.null(keep_strats)) {
        new_samples_freq <- new_samples_freq |> 
            filter(strategies %in% keep_strats)
    }

    if (!is.null(keep_sizes)) {
        new_samples_freq <- new_samples_freq |> 
            filter(new_samples_size %in% keep_sizes)
    }


    plot_grid <- ggplot_quantitative_df_on_background_map(
        background_map = available_points,
        df = new_samples_freq |> filter(
            r_effect == R_EFFECT_PLOT),
        lon_c = "lon", lat_c = "lat",
        column = "n_folds", unit = "# of selections",
        facet_formula = new_samples_size ~ strategies
    )
    print(plot_grid)
    if (!is.null(save_to)) {
        standardised_ggplot_save(plot_grid, save_to)
    }
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


##### Quick comparison of random effect influence ##### -----------------------
dummy_scores_df <- scores_df
dummy_scores_df <- dummy_scores_df |>
    bind_rows(
        dummy_scores_df |>
            filter(name == "RMSE") |>
            mutate(name = "MSE", score = score^2)
        ) |>
    filter(
        name == "MSE",
        NEW_SAMPLE_SIZES == 0,
        STRATEGIES == "none")

dummy_summary_df <- dummy_scores_df |>
    group_by(R_EFFECTS, subset) |>
    summarise(
        mean_mse = mean(score),
        sd_mse = sd(score),
        .groups = "drop_last")

effect_plot <- ggplot(
        dummy_summary_df,
        aes(x = subset, y = mean_mse, colour = R_EFFECTS)) +
    geom_pointrange(
        aes(ymin = mean_mse - sd_mse, ymax = mean_mse + sd_mse),
            position = position_dodge(width = 0.4)) +
    labs(
        x = NULL, y = "MSE (mean ± SD)", 
        colour = "Random effects",
        title = "Influence of random effects",
        caption = "Comparing base models (fitted on BIODICAPT sites only)") +
    theme_minimal()
effect_plot <- my_custom_ggplot_theme(effect_plot, light = TRUE)
standardised_ggplot_save(
    effect_plot, file.path(ANALYSES_PATH, "random_effects.pdf"))

# Conclusion : 
#   - shows small gains in MSE -> latent variables are doing their jobs
#   - keeping "site" random effect to improve prediction


##### Preprocessing scores ##### ----------------------------------------------
# Select fixed effects of interest
scores_df <- scores_df |> filter(R_EFFECTS == R_EFFECT_PLOT)
# Select strategies of interest
scores_df <- scores_df |> filter(STRATEGIES %in% ORDER_STRATS)

# Add MSE and rename metrics
scores_df <- scores_df |>
    bind_rows(
        scores_df |>
            filter(name == "RMSE") |>
            mutate(name = "MSE", score = score^2)
        ) |>
    mutate(name = recode(
        name,
        "SD" = "SD of predictions",
        "BetaSD" = "SD of parameters (all)",
        "BetaSD_actionable" = "SD of actionable parameters",
        "BetaSD_ratio" = "SD ratio of parameters (all)",
        "BetaSD_actionable_ratio" = "SD ratio of actionable parameters",
        "Ds" = "Ds criterion (realised)",
        "I" = "I criterion (realised)",
        "Ds_expected" = "Ds criterion (expected)",
        "I_expected" = "I criterion (expected)"))
# Add gain in MSE per number of added samples
scores_df <- scores_df |>
    bind_rows(
        scores_df |>
            filter(name == "MSE") |>
            mutate(
                name = "MSE per sample", 
                score = score/NEW_SAMPLE_SIZES,
                score = if_else(is.infinite(score), 0, score)
            )
        )

# Partition reference effects vs other effects
ref_df <- scores_df |>
    semi_join(as_tibble(REFERENCE_EFFECTS), by = names(REFERENCE_EFFECTS))
compare_df <- scores_df |>
    anti_join(
        as_tibble(REFERENCE_EFFECTS), 
        by = names(REFERENCE_EFFECTS)
)

# Add number of species improving by >1% in MSE 
# (to compare with doi.org/10.1111%2F2041-210X.14355)
improvement_df <- count_improving(
    ref_df, compare_df,
    metric = "RMSE",
    threshold = 0.01,
    group_vars = c("STRATEGIES", "NEW_SAMPLE_SIZES"),
    panel_var = "subset")
improvement_df["name"] <- "# of species with increase in ΔMSE >1%"
improvement_df <- improvement_df |> rename(diff_score = n_sp_above)


##### Check criteria approximations ##### -------------------------------------
# Rescaled because extremums in NEW_SAMPLE_SIZES 10 and 150 skew the plots
ds_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "Ds criterion (expected)",
    y_metric = "Ds criterion (realised)",
    title = "Expected vs realised Ds criterion",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    panel_order = c("train"),
    lower_is_better = FALSE,
    save_to = file.path(ANALYSES_PATH, "ds-criterion_expected_vs_realised.pdf"))
ds_corr_plot$plot <- ds_corr_plot$plot +
    coord_cartesian(ylim = c(6, 11))
print(ds_corr_plot$plot)
standardised_ggplot_save(
    ds_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_realised.pdf"))

dsparam_corr_plot <- corr_metrics(
    scores = compare_df,  
    x_metric = "Ds criterion (expected)",
    y_metric = "SD of actionable parameters",
    title = "Expected Ds criterion VS SD of actionable parameters",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    panel_order = c("train"),
    lower_is_better = FALSE,
    save_to = NULL)
dsparam_corr_plot$plot <- dsparam_corr_plot$plot +
    labs(subtitle = "(Ds: higher is better, SD: lower is better)") +
    coord_cartesian(ylim = c(0.067, 0.2))
print(dsparam_corr_plot$plot)
standardised_ggplot_save(
    dsparam_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_actionable_SD.pdf"))

dsparam_ratio_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "Ds criterion (expected)",
    y_metric = "SD ratio of actionable parameters",
    title = "Expected Ds criterion VS SD ratio of actionable parameters",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    panel_order = c("train"),
    lower_is_better = FALSE,
    save_to = NULL)
dsparam_ratio_corr_plot$plot <- dsparam_ratio_corr_plot$plot + 
    labs(subtitle = "(Ds: higher is better, SD ratio: lower is better)") +
    coord_cartesian(ylim = c(0.33, 1))
print(dsparam_ratio_corr_plot$plot)
standardised_ggplot_save(
    dsparam_ratio_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_actionable_SD_ratio.pdf"))

dsallparam_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "Ds criterion (expected)",
    y_metric = "SD of parameters (all)",
    title = "Expected Ds criterion VS SD of all parameters",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    panel_order = c("train"),
    lower_is_better = FALSE,
    save_to = NULL)
dsallparam_corr_plot$plot <- dsallparam_corr_plot$plot + 
    labs(subtitle = "(Ds: higher is better, SD: lower is better)") +
    coord_cartesian(ylim = c(0.095, 0.2))
print(dsallparam_corr_plot$plot)
standardised_ggplot_save(
    dsallparam_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_all_parameters_SD.pdf"))

dsallparam_ratio_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "Ds criterion (expected)",
    y_metric = "SD ratio of parameters (all)",
    title = "Expected Ds criterion VS SD ratio of all parameters",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    panel_order = c("train"),
    lower_is_better = FALSE,
    save_to = NULL)
dsallparam_ratio_corr_plot$plot <- dsallparam_ratio_corr_plot$plot + 
    labs(subtitle = "(Ds: higher is better, SD ratio: lower is better)") +
    coord_cartesian(ylim = c(0.45, 1))
print(dsallparam_ratio_corr_plot$plot)
standardised_ggplot_save(
    dsallparam_ratio_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_all_parameters_SD_ratio.pdf"))

i_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "I criterion (expected)",
    y_metric = "I criterion (realised)",
    title = "Expected vs realised I criterion",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    lower_is_better = FALSE,
    panel_order = c("test_local"),
    save_to = NULL)
i_corr_plot$plot <- i_corr_plot$plot +
    coord_cartesian(ylim = c(-0.045, -0.0033))
print(i_corr_plot$plot)
standardised_ggplot_save(
    i_corr_plot$plot, 
    file.path(ANALYSES_PATH, "i-criterion_expected_vs_realised.pdf"))

ipred_corr_plot <- corr_metrics(
    scores = compare_df, 
    x_metric = "I criterion (expected)",
    y_metric = "SD of predictions",
    title = "Expected I criterion vs SD of predictions",
    strat_var = "STRATEGIES", 
    strat_order = ORDER_STRATS,
    size_var = "NEW_SAMPLE_SIZES", 
    size_order = NULL,
    panel_var = "subset", 
    lower_is_better = TRUE,
    panel_order = c("test_local"),
    save_to = NULL)
ipred_corr_plot$plot <- ipred_corr_plot$plot + 
    labs(subtitle = "(I: higher is better, SD: lower is better)")
print(ipred_corr_plot$plot)
standardised_ggplot_save(
    ipred_corr_plot$plot, 
    file.path(ANALYSES_PATH, "i-criterion_expected_vs_predictions.pdf"))


##### Check strategy effect on estimations/predictions ##### ------------------
# metric name = list(pdf name, lower_is_better boolean)
dist_metrics <- list(
    "SD of predictions" = c("Delta_SD_distrib.pdf", TRUE),
    "SD of parameters (all)" = c("Delta_SD_parameters_all.pdf", TRUE),
    "SD of actionable parameters" = c("Delta_SD_parameters_actionable_distrib.pdf", TRUE),
    "MSE" = c("Delta_MSE_distrib.pdf", TRUE),
    "MSE per sample" = c("Delta_MSE_per_sample_distrib.pdf", FALSE),
    "TjurR2" = c("Delta_TjurR2_distrib.pdf", FALSE)
)

for (m in names(dist_metrics)) {
    if (grepl("param", m)) {
        subsets_plot <- SUBSETS_PLOT[1]
    } else {
        subsets_plot <- SUBSETS_PLOT
    }
    metric_lineplot <- compare_lineplot(
        ref_scores = ref_df,
        compare_scores = compare_df,
        metric = m,
        title = NULL, 
        x_var = "NEW_SAMPLE_SIZES",
        x_order = NULL,
        strat_var = "STRATEGIES", 
        strat_order = ORDER_STRATS,
        panel_var = "subset", 
        panel_order = subsets_plot,
        lower_is_better = dist_metrics[[m]][2],
        save_to = file.path(ANALYSES_PATH, dist_metrics[[m]][1]))
    print(metric_lineplot$plot)
}


##### Compare with doi.org/10.1111%2F2041-210X.14355 ##### --------------------
lineplot_mse <- compare_lineplot(
    ref_scores = ref_df,
    diffs_scores = compare_df |> 
        rename(diff_score = score),
    metric = "MSE",
    x_var = "NEW_SAMPLE_SIZES", 
    x_order = NULL,
    strat_var = "STRATEGIES",
    strat_order = ORDER_STRATS,
    panel_var = "subset", 
    panel_order = SUBSETS_PLOT,
    random_strategy = "random-sampling",
    group_species = TRUE,
    summary_fun = mean,
    lower_is_better = TRUE,
    save_to = file.path(ANALYSES_PATH, "MSE.pdf"))
print(lineplot_mse$plot)
# remove error bars for readability
lineplot_mse$plot$layers <- Filter(
    function(l) !inherits(l$geom, "GeomErrorbar"), lineplot_mse$plot$layers)
print(lineplot_mse$plot)
standardised_ggplot_save(
    lineplot_mse$plot, 
    file.path(ANALYSES_PATH, "MSE_no_bars.pdf"))

lineplot_n_mse <- compare_lineplot(
    diffs_scores = improvement_df,
    metric = unique(improvement_df["name"])[[1]],
    x_var = "NEW_SAMPLE_SIZES", 
    x_order = NULL,
    strat_var = "STRATEGIES",
    strat_order = ORDER_STRATS,
    panel_var = "subset", 
    panel_order = SUBSETS_PLOT,
    random_strategy = "random-sampling",
    group_species = TRUE,
    summary_fun = mean,
    lower_is_better = FALSE,
    save_to = file.path(ANALYSES_PATH, "MSE_delta_n_sp.pdf"))
print(lineplot_n_mse$plot)


##### MAPS ###### -------------------------------------------------------------
### Map frequency of sampled locations in new_pool
all_datasets <- unwrap_simulations(
    readRDS(file = file.path(RESULTS_PATH, "datasets.rds")))
show_map_selected_new_pool(
    complete_new_pool = all_datasets[[1]]$new_pool,
    keep_strats = c(
        "by-region", "max-uncertainty", "ds-optimality", "i-optimality"),
    keep_sizes = c(10, 50, 150),
    save_to = file.path(ANALYSES_PATH, "freq-new-locations.pdf"))