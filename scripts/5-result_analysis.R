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
    "random-sampling", "gap-filling", "max-uncertainty", 
    "2-part-max-uncertainty", "5-part-max-uncertainty",
    "10-part-max-uncertainty", "ds-optimality", "i-optimality")

# FIX Random effect (we are not interested in its effect on predictions)
R_EFFECT_PLOT <- "site"
REFERENCE_EFFECTS$R_EFFECTS <-  R_EFFECT_PLOT

# Plot parameters
SUBSETS_PLOT <- c("train", "test_local", "test_stoc") # c("train", "test_local", "test_stoc")
NEW_SAMPLE_SIZE_MAP <- 50
STRATEGY_MAP <- "max-uncertainty"


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
            all_scores <- rbind(
                all_scores, 
                format_scores(
                    scores, local_parameters, names(local_datasets)[j]))
            
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

show_map_selected_new_pool <- function(complete_new_pool, save_to = NULL) {
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

show_new_pool_distributions <- function(
        complete_new_pool, variable,
        strat_var = "strategies", strat_order = NULL,
        size_var  = "new_samples_size", size_order = NULL,
        fold_var  = "k_folds",
        random_strategy = "random-sampling",
        aggregate = c("fold_mean", "pool"),
        summary_fun = mean,
        max_ref_n = 5000, seed = 1,
        save_to = NULL) {

    var_levels <- variable   # keep the argument safe from column masking
    aggregate <- match.arg(aggregate)
    set.seed(seed)

    # 1. Selected points + features
    new_points_df <- fetch_all_new_samples_id() |>
        left_join(
            as.data.frame(complete_new_pool)[c(variable, "row_id")],
            by = c("sample_id" = "row_id")) |>
        pivot_longer(all_of(variable), names_to = "variable", values_to = "value") |>
        filter(!is.na(value))

    # 2. Aggregation over folds
    if (aggregate == "fold_mean") {
        new_points_df <- new_points_df |>
            group_by(across(all_of(c(fold_var, size_var, strat_var, "variable")))) |>
            summarise(value = summary_fun(value, na.rm = TRUE), .groups = "drop")
    }

    # 3. References: France + pool of candidates (subsampled)
    france_df <- clip_raster_to_france_shape(rast(FEATURES_RES_TRF), buffer = 0) |>
        as.data.frame() |> fix_names()
    to_long <- function(df, label) {
        df[variable] |>
            pivot_longer(everything(), names_to = "variable", values_to = "value") |>
            filter(!is.na(value)) |>
            mutate(ref = label)
    }
    ref_long <- bind_rows(
        to_long(france_df, "France"),
        to_long(as.data.frame(complete_new_pool), "New points (500 ENI)"))

    # 4. Factor levels and colours
    lvls <- .strategy_levels(new_points_df[[strat_var]], strat_order, random_strategy)
    colors <- .strategy_colours(lvls, random_strategy, random_colour = "grey50")
    ref_lvls <- c("France", "New points (500 ENI)")
    ref_cols <- c("France" = "grey15", "New points (500 ENI)" = "grey85")
    all_lvls <- c(ref_lvls, lvls)
    all_cols <- c(ref_cols, colors)

    new_points_df <- new_points_df |>
        .set_factor(strat_var, all_lvls) |>
        .set_factor(size_var, size_order) |>
        mutate(variable = factor(variable, levels = var_levels),
               is_ref = FALSE)
    size_levels <- levels(new_points_df[[size_var]])

    # 5. References as extra "strategies", repeated in every size column
    ref_rows <- ref_long |>
        group_by(variable, ref) |>
        filter(row_number() %in% sample.int(n(), min(max_ref_n, n()))) |>
        ungroup() |>
        rename(!!strat_var := ref) |>
        tidyr::crossing(!!size_var := factor(size_levels, levels = size_levels)) |>
        mutate(!!strat_var := factor(.data[[strat_var]], levels = all_lvls),
               variable = factor(variable, levels = var_levels),
               is_ref = TRUE)
    plot_df <- bind_rows(new_points_df, ref_rows)

    # 6. Plot: strategies on x, values on y
    p <- ggplot(plot_df, aes(x = .data[[strat_var]], y = value,
                             fill = .data[[strat_var]])) +
        geom_violin(scale = "width", trim = TRUE, adjust = 0.5,
                    colour = "grey30", alpha = 0.7, linewidth = 0.4)
    
    if (aggregate == "pool") {
        p <- p + geom_boxplot(width = 0.05, outlier.shape = NA,
                              show.legend = FALSE)
    } else {
        p <- p +
            geom_jitter(data = \(d) filter(d, !is_ref),
                        width = 0.12, height = 0, size = 0.8,
                        alpha = 0.5, colour = "grey20",
                        show.legend = FALSE) +
            stat_summary(fun = summary_fun, geom = "point",
                         shape = 23, size = 2.5, fill = "white", colour = "black",
                         show.legend = FALSE)
    }

    p <- p +
        facet_grid(rows = vars(variable), cols = vars(.data[[size_var]]),
                   scales = "free_y") +   # y shared within a row; x identical everywhere
        labs(title = NULL, x = NULL, y = NULL,
             caption = if (aggregate == "pool")
                 "Selected points pooled over k folds; references: all pixels / candidates" else
                     "Strategies: one value per k fold; references: all pixels / candidates")

    p <- my_custom_ggplot_theme(p, with_palette = FALSE) +
        scale_fill_manual(values = all_cols, breaks = all_lvls, name = NULL) +
        scale_x_discrete(limits = all_lvls) +   # references first, then strategies
        theme(legend.position = "bottom",
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())

    print(p)
    if (!is.null(save_to)) {
        standardised_ggplot_save(p, save_to)
    }
    invisible(list(data = plot_df, plot = p))
}

LEGACY_plot_ridgeline <- function(
        dfs, labels, group_folds,
        zero_inflated = NULL,       # names of the variables to transform
        zero_transform = sqrt,      # function applied to those variables
        zero_tol = 1e-5,
        bandwidth_adjust = 1, clip_quantiles = c(0.01, 0.99)) {
    library(ggridges)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    # 1. Security check
    n <- length(dfs)
    if (!all(c(length(labels), length(group_folds)) == n)) {
        stop("`labels` and `group_folds` must have the same length as `dfs`.")
    }

    labels <- as.list(labels)
    group_folds <- as.list(group_folds)


    # 2. Transform list of dataframes into a list of df with a "source" column
    prepare_source <- function(df, label, pool) {
        if ("k_folds" %in% names(df)) {
            long <- to_long(df) |> 
                mutate(k = as.character(k_folds))
            if (isTRUE(pool)) {
                long$source <- label
                lvls <- label
            } else {
                ks   <- as.character(sort(unique(as.numeric(long$k))))
                lvls <- paste(label, ks)
                long$source <- paste(label, long$k)
            }
        } else {
            long <- to_long(df, label)
            long$source <- label
            lvls <- label
        }

        # Colours: filled in later
        n_lv <- length(lvls)
        cols <- rep(NA_character_, n_lv)

        list(
            data = dplyr::select(long, source, variable, value),
            levels = lvls,
            colors = setNames(cols, lvls)
        )
    }
    prepared <- Map(prepare_source, dfs, labels, group_folds)

    # 3. Combine dataframes together
    all_levels <- unlist(
        lapply(prepared, `[[`, "levels"), 
        use.names = FALSE)
    if (anyDuplicated(all_levels)) {
        stop("Labels must be unique across `dfs`.")
    }

    all_dfs <- bind_rows(lapply(prepared, `[[`, "data")) |>
        mutate(source = factor(source, levels = all_levels))

    # 4. Transform zero-inflated variables (before rescaling)
    if (!is.null(zero_inflated)) {
        unknown <- setdiff(zero_inflated, unique(as.character(all_dfs$variable)))
        if (length(unknown) > 0) {
            warning("Unknown variables in `zero_inflated`: ",
                    paste(unknown, collapse = ", "))
        }

        # 4a. Flag zeros (the reference level is shared by all sources)
        all_dfs <- all_dfs |>
            group_by(variable) |>
            mutate(
                is_zi    = as.character(variable) %in% zero_inflated,
                zero_lvl = min(value, na.rm = TRUE),
                is_zero  = is_zi & value <= zero_lvl + zero_tol
            ) |>
            ungroup()

        # 4b. One label per variable, with the share of zeros per source
        zero_labels <- all_dfs |>
            filter(is_zi) |>
            group_by(variable, source) |>
            summarise(
                pct_zero = mean(is_zero, na.rm = TRUE),
                .groups  = "drop"
            ) |>
            group_by(variable) |>
            summarise(
                zero_label = paste0(
                    source, " (", round(100 * pct_zero), "%)",
                    collapse = ", "
                ),
                .groups = "drop"
            )

        # 4c. Drop the spike, transform the rest, attach the labels
        all_dfs <- all_dfs |>
            filter(!is_zero) |>
            group_by(variable) |>
            mutate(
                value = if (first(is_zi)) zero_transform(value - zero_lvl) else value
            ) |>
            ungroup() |>
            left_join(zero_labels, by = "variable") |>
            mutate(
                variable = ifelse(
                    is_zi,
                    paste0(variable, " *\nshare of zeros: ", zero_label),
                    as.character(variable))
            ) |>
            dplyr::select(-is_zi, -zero_lvl, -is_zero, -zero_label)
    }
    
    # 5. Rescale per variables
    all_dfs <- all_dfs |>
        group_by(variable) |>
        mutate(
            lo = quantile(value, clip_quantiles[1], na.rm = TRUE),
            hi = quantile(value, clip_quantiles[2], na.rm = TRUE),
            # fall back to the full range if the quantiles collapse
            hi = ifelse(hi == lo, max(value, na.rm = TRUE), hi),
            value_scaled = pmin(pmax((value - lo) / (hi - lo), 0), 1)) |>
        ungroup() |>
        mutate(variable = factor(variable, levels = rev(unique(variable))))

    # 6. Make plot
    p <- ggplot(
            all_dfs, 
            aes(x = value_scaled, 
                y = variable, 
                fill = source, 
                colour = source)
            ) +
        geom_density_ridges(
            stat   = "density",
            aes(height = after_stat(scaled)),
            alpha = 0.35,
            scale = 0.9,
            adjust = bandwidth_adjust,
            trim = FALSE) +
        scale_x_continuous(
            expand = c(0, 0), 
            breaks = c(0, 1),
            labels = c("low", "high")) +
        labs(
            x = "Value (rescaled per variable)", y = NULL,
            caption = if (!is.null(zero_inflated)) {
                "*: zeros excluded from the distribution (in brackets: number of points kept brackets)."
            } else {
                NULL
            }
        )

    my_custom_ggplot_theme(p, with_palette = TRUE) +
        guides(
            fill = guide_legend(title = "Data source:"),
            colour = guide_legend(title = "Data source:"))
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

# Select fixed effects of interest
scores_df <- scores_df |> filter(R_EFFECTS == R_EFFECT_PLOT)

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
        "Ds" = "Ds criterion (realised)",
        "I" = "I criterion (realised)",
        "Ds_expected" = "Ds criterion (expected)",
        "I_expected" = "I criterion (expected)"))

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
print(ds_corr_plot$plot)

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
    labs(subtitle = "(Ds: higher is better, SD: lower is better)")
print(dsparam_corr_plot$plot)
standardised_ggplot_save(
    dsparam_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_actionable_SD.pdf"))

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
    labs(subtitle = "(Ds: higher is better, SD: lower is better)")
print(dsallparam_corr_plot$plot)
standardised_ggplot_save(
    dsallparam_corr_plot$plot, 
    file.path(ANALYSES_PATH, "ds-criterion_expected_vs_all_parameters_SD.pdf"))

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
    lower_is_better = TRUE,
    panel_order = c("test_local", "test_stoc"),
    save_to = file.path(ANALYSES_PATH, "i-criterion_expected_vs_realised.pdf"))
print(i_corr_plot$plot)

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
    panel_order = c("test_local", "test_stoc"),
    save_to = file.path(ANALYSES_PATH, "i-criterion_expected_vs_predictions.pdf"))
print(ipred_corr_plot$plot)


##### Check strategy effect on estimations/predictions ##### ------------------
dist_metrics <- c(
    "SD of predictions" = "SD_distrib.pdf",
    "SD of parameters (all)" = "SD_parameters_all.pdf",
    "SD of actionable parameters" = "SD_parameters_actionable_distrib.pdf",
    "MSE" = "MSE_distrib.pdf")

for (m in names(dist_metrics)) {
    if (grepl("param", m)) {
        subsets_plot <- SUBSETS_PLOT[1]
    } else {
        subsets_plot <- SUBSETS_PLOT
    }
    strat_violin_plot <- distribution_plot(
        scores = compare_df,
        ref_scores = ref_df,
        metric = m,
        title = NULL, 
        strat_var = "STRATEGIES", 
        strat_order = ORDER_STRATS,
        size_var = "NEW_SAMPLE_SIZES", 
        size_order = NULL,
        panel_var = "subset", 
        panel_order = subsets_plot,
        lower_is_better = TRUE,
        ref_label = "Reference model",
        save_to = file.path(
            ANALYSES_PATH, paste0("violinplots_", dist_metrics[[m]])))
    print(strat_violin_plot$plot)
}


##### Compare with doi.org/10.1111%2F2041-210X.14355 ##### --------------------
violin_mse <- compare_plot(
    ref_scores = ref_df,
    compare_scores = compare_df,
    metric = "MSE",
    x_var = "NEW_SAMPLE_SIZES",
    x_order = NULL,
    strat_var = "STRATEGIES",
    strat_order = ORDER_STRATS,
    panel_var = "subset",
    panel_order = SUBSETS_PLOT,
    group_species = TRUE,
    lower_is_better = TRUE,
    save_to = file.path(ANALYSES_PATH, "boxplot_MSE.pdf"))
print(violin_mse$plot)

violin_D_mse <- compare_plot(
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
    show_folds = TRUE,
    summary_fun = mean,
    lower_is_better = FALSE,
    save_to = file.path(ANALYSES_PATH, "MSE_delta_n_sp.pdf"))
print(violin_D_mse$plot)


##### MAPS ###### -------------------------------------------------------------
### Map frequency of sampled locations in new_pool
all_datasets <- unwrap_simulations(
    readRDS(file = file.path(RESULTS_PATH, "datasets.rds")))
show_map_selected_new_pool(
    complete_new_pool = all_datasets[[1]]$new_pool,
    save_to = file.path(ANALYSES_PATH, "freq-new-locations.pdf"))


#### VAR DISTRIBUTIONS ##### --------------------------------------------------
show_new_pool_distributions(
        complete_new_pool = all_datasets[[1]]$new_pool, 
        variable = c("light.pollution", "swf"),
        strat_var = "strategies", 
        strat_order = ORDER_STRATS,
        size_var = "new_samples_size", 
        size_order = NULL,
        aggregate   = "pool",
        random_strategy = "random-sampling",
        save_to = file.path(ANALYSES_PATH, "distributions_grid.pdf")) 

cli_alert_success("Analysis of HMSC is complete!\n\n")