# This script:
#   - Analyses the results of the I-optimisation runs

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(readr)
library(dplyr)
library(purrr)

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
    stop("User refused to overwrite previous analyses of model's runs.")
}


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


##### Fetch scores ##### ------------------------------------------------------
# Compute models scores
scores_df <- fetch_all_scores()
# Save them
write_csv(scores_df, file.path(ANALYSES_PATH, "all_scores.csv"))

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
        fill_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset",
        panel_order = c("train", "new_pool", "test"),
        group_species = TRUE,
        save_to = file.path(ANALYSES_PATH, "SD_deltas.pdf"))


##### MAPS ###### -------------------------------------------------------------
# TODO :
#   - A map showing the frequency of selection of each point in new_pool
#   - Grid of histograms showing the features of the most selected points