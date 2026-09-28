# Set of utility functions used to analyse data.

##### Liraries #####
library(cli)
library(ade4)
library(uwot)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(GGally)
library(ggplot2)
library(parallel)
library(tidyverse)
library(factoextra)
library(ggcorrplot) 
library(RColorBrewer)


##### Parameters #####
source(here::here("data/config/config.R")) # all parameters are grouped together
source(here::here("R/utils_figures.R"))  # needs my_custom_ggplot_theme, standardised_ggplot_save


##### Global functions #####
# A function to summarise the significance of a p-value
# ARGS:
#   - p: a numeric, the p-value to evaluate.
# RETURNS:
#   - a string with stars to show statistical significance
sig_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1)  "·"
    else ""
}

# A function to choose which correlation to apply to a association of variables
# when using GGally::ggpairs
# ARGS: what's expected by GGpairs
assoc_fun <- function(data, mapping, ...) {
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y) 

    if (is.numeric(x) && is.numeric(y)) {
        # both continuous -> Pearson correlation
        test <- cor.test(x, y, method = "pearson")
        r <- test$estimate
        stars <- sig_stars(test$p.value)
        lbl <- paste0("Corr:\nr = ", round(r, 2), stars)

    } else if (is.numeric(x) != is.numeric(y)) {
        # one continuous, one categorical -> correlation ratio (eta)
        if (is.numeric(x)) { 
            num_var <- x; cat_var <- factor(y) 
        } else { 
            num_var <- y; cat_var <- factor(x) 
        }

        fit <- aov(num_var ~ cat_var)
        fit_summary <- summary(fit)[[1]]
        eta2 <- fit_summary[1, "Sum Sq"] / sum(fit_summary[, "Sum Sq"])
        p_val <- fit_summary[1, "Pr(>F)"]
        stars <- sig_stars(p_val)
        lbl <- paste0("Corr ratio:\nη = ", round(sqrt(eta2), 2), stars)

    } else {
        # both categorical -> Cramer's V
        tbl <- table(x, y)
        chi <- suppressWarnings(chisq.test(tbl))
        n   <- sum(tbl)
        V   <- sqrt((chi$statistic / n) / min(nrow(tbl) - 1, ncol(tbl) - 1))
        stars <- sig_stars(chi$p.value)
        lbl <- paste0("Cramer:\nV = ", round(V, 2), stars)
    }

    ggally_text(label = lbl, mapping = aes(), color = "black", ...) +
    theme_void()
}

# A function to add a ellipse and loess to plotted points in ggpair diag
# ARGS: what's expected by GGpairs
lower_cont_fun <- function(data, mapping, ..., max_n = 2000) {
    # if too many rows, sample down for faster computation
    if (nrow(data) > max_n) {
        data <- data[sample(nrow(data), max_n), ]
    }
    ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.4, color = "grey20") +
        stat_ellipse(color = PALETTE[2], type = "norm", level = 0.95, linewidth = 0.6) +
        geom_smooth(method = "loess", color = PALETTE[1], se = FALSE,
                    linewidth = 0.6, n = 50, ...)
}

# A function that create's a draftman's plot.
# ARGS:
#   - df: a dataframe.
#   - columns: a vector. 
#       Contains the names of the columns to plot against each other.
# RETURNS:
#   - a ggplot figure
ggplot_custom_draftman <- function(df, columns) {
    plot <- ggpairs(df,
        columns = columns,
        progress = TRUE,
        upper = list(continuous = assoc_fun, combo = assoc_fun, discrete = assoc_fun),
        lower = list(
            continuous = lower_cont_fun, 
            combo = "box_no_facet",
            discrete = "count"
        )
    )

    plot <- my_custom_ggplot_theme(plot) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
    return(plot)
}

# A function that computes a generalized "correlation" matrix 
# (Pearson r / eta / Cramer's V)
# ARGS:
#   - df: a dataframe.
#   - columns: a vector. 
#       Contains the names of the columns to plot against each other.
# RETURNS:
#   - a matrix
get_assoc_matrix <- function(df, columns) { 
    n <- length(columns)
    mat <- matrix(NA_real_, n, n, dimnames = list(columns, columns))

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
        if (i == j) {
            mat[i, j] <- 1
            next
        }
        if (j < i) next  # fill lower triangle by symmetry at the end

        x <- df[[columns[i]]]
        y <- df[[columns[j]]]

        if (is.numeric(x) && is.numeric(y)) {
            val <- suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs"))

        } else if (is.numeric(x) != is.numeric(y)) {
            if (is.numeric(x)) { num_var <- x; cat_var <- factor(y) } 
            else               { num_var <- y; cat_var <- factor(x) }
            fit <- aov(num_var ~ cat_var)
            fit_summary <- summary(fit)[[1]]
            eta2 <- fit_summary[1, "Sum Sq"] / sum(fit_summary[, "Sum Sq"])
            val <- sqrt(eta2)

        } else {
            tbl <- table(x, y)
            chi <- suppressWarnings(chisq.test(tbl))
            val <- sqrt((chi$statistic / sum(tbl)) / min(nrow(tbl) - 1, ncol(tbl) - 1))
        }

        mat[i, j] <- val
        mat[j, i] <- val  # symmetric
        }
    }
    mat
}

# A function that summarises the contents of a dataset with X variables and 
# Y species. Saves a dataframe of occurences and a draftman's plot in PDF.
# ARGS:
#   - df: a data.frame 
#   - x_cols: a list of strings. The columns containing explanatory variables.
#   - save_folder: a string. Path to a folder where "draftman_plot.pdf" will be saved.
#   - sp_cols: a list of strings (optional, default is NULL). 
#       The columns containing species occurrences.
#   - top: a numeric (default is 5). Controls the number of species to show as most and least represented in dataset
# RETURNS:
#   - NULL, prints results in console and silently saves .CSV and .PDF
explore_dataset <- function(
    df, x_cols, save_folder, save_name, sp_cols = NULL, top = 5) {
    ### PLOT
    cli_alert_info("Exploration of dataset.")

    # Make sure df is a dataframe
    if (inherits(df, "SpatVector")) {
        df <- as.data.frame(df)
    } else if (!inherits(df, "data.frame")) {
        stop(paste0("df must be a dataframe, got ", class(df), "."))
    }
   
    # For many columns, draftman is way too slow, use a simple correlation plot
    if (length(x_cols) > 7) {
        corr_plot <- ggcorrplot(
                get_assoc_matrix(df[x_cols], columns = x_cols), 
                lab = TRUE, 
                type = "lower", 
                colors = c("#5793cf", "white", "#EE6677"))
        corr_plot <- my_custom_ggplot_theme(corr_plot, with_palette = FALSE) + 
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        suppressMessages(print(corr_plot))
        standardised_ggplot_save(
            figure = corr_plot, 
            save_path = file.path(
                save_folder, paste0(save_name, "_corr_plot.pdf")),
            .width = 36, .height = 12)
        cli_alert_success("Saved correlation plot.")

    } else {
        draft_plot <- ggplot_custom_draftman(df, columns = x_cols)
        suppressMessages(print(draft_plot))
        standardised_ggplot_save(
            figure = draft_plot, 
            save_path = file.path(
                save_folder, paste0(save_name, "_draftman_plot.pdf")),
            .width = 36, .height = 12)
        cli_alert_success("Saved Draftman's plot.")
    }

    ### FACTORS
    factor_df <- df |>
        select(all_of(x_cols)) |>
        select(where(~ !is.numeric(.x)))

    for (col in names(factor_df)) {
        # for each column, compute frequence and proportion of each variable
        freq <- table(factor_df[[col]])
        prop <- round(freq/sum(freq), digits=3)
        combined_table <- rbind(freq = as.vector(freq), 
                            prop = as.vector(prop))

        colnames(combined_table) <- names(freq)

        cli_alert_info(paste0("Distribution of '",  col ,"' in dataset:"))
        print(t(combined_table))
        cli_alert_info(" ")
    }

    ### PREDICTED VARIABLES
    if (!is.null(sp_cols)) {
        # Make a table with the number of occurences
        top_sp <- df[if (!is.null(sp_cols)) sp_cols else NAMES_SPECIES] |>
            summarise(across(where(is.numeric), sum)) |>
            pivot_longer(everything(), names_to = "column", values_to = "sum")  |>
            mutate(frequency = sum / nrow(df))

        # print top-x species and bottom-x species
        cli_alert_info(paste0(
            "Top-", top, " MOST sighted species (from occurences):"))
        print(as.data.frame(top_sp |> slice_max(sum, n = top)))
        cli_alert_info(" ")
        cli_alert_info(paste0(
            "Top-", top," LEAST sighted species (from occurences):"))
        print(as.data.frame(top_sp |> slice_min(sum, n = top)))
        write_csv(
                as.data.frame(top_sp), 
                file.path(
                    save_folder, 
                    paste0(save_name, "_occurences_in_full_dataset.csv")))
        cli_alert_success("Saved table of occurrences in dataset.")
    }

    cat("\n")
    return(NULL)
}

# A function that automatically selects numerical columns in a dataframe, 
# standardizes them and computes a PCA
# ARGS:
#   - df: a dataframe (or SpatVector object)
#   - color_var: the name of a column of categorical values in df
# RETURNS: a list:
#   - pca: the PCA result values
#   - p_ind: the PCA plot showing individual observations in 2D
#   - p_var: the contribution plot for all variables, to each axis in 2D
make_pca_and_plots <- function(df, color_var = NULL) {
    # Make sure df is a dataframe
    if (inherits(df, "SpatVector")) {
        df <- as.data.frame(df)
    } else if (!inherits(df, "data.frame")) {
        stop(paste0("df must be a dataframe, got ", class(df), "."))
    }

    # Make sure NAs are excluded
    df_clean <- df |> drop_na()

    # Keep only numeric columns for PCA
    df_num <- df_clean |> select(where(is.numeric))

    # compute PCA (always standardize data)
    pca_res <- dudi.pca(
        df_num, center = TRUE, scale = TRUE, scannf = FALSE, nf = 5)

    # Individuals plot (samples), colored by a grouping variable if provided
    p_ind <- fviz_pca_ind(
            pca_res,
            geom.ind = "point",
            col.ind = if (!is.null(color_var)) df_clean[[color_var]] else "grey20",
            alpha.ind = 0.4,
            palette = "jco",
            addEllipses = !is.null(color_var),
            legend.title = color_var %||% "Group",
            repel = TRUE) +
        theme_minimal(base_size = 13) +
        labs(title = "PCA - Individuals")
    p_ind <- my_custom_ggplot_theme(p_ind)

    # Variables plot (loadings / contributions)
    p_var <- fviz_pca_var(
            pca_res,
            col.var = "contrib",
            gradient.cols = rev(brewer.pal(11, "RdYlBu")),
            repel = TRUE
        ) +
        theme_minimal(base_size = 13) +
        labs(title = "PCA - Variables")
    p_var <- my_custom_ggplot_theme(p_var, with_palette = FALSE)

    return(list(pca = pca_res, plot_ind = p_ind, plot_var = p_var))
}

# A function to compute ellipses in the same way as "stat_ellipse" in 
# make_pca_and_plots
# ARGS:
#   - df: a dataframe (or SpatVector object)
#   - x, y: strings. Names of the first and second dimension 
#   - level: a numeric between 0 and 1. The confidence level for the ellipse
#   - segments: an integer. The number of segments used to draw the ellipse
# RETURNS:
#   - a dataframe of points forming an ellipse
compute_ellipse <- function(df, x = "Dim.1", y = "Dim.2", level = 0.95, segments = 51) {
    data <- df[, c(x, y)]
    dfn <- 2
    dfd <- nrow(data) - 1

    if (dfd < 3) {
        warning("Too few points to calculate an ellipse")
        return(data.frame(setNames(list(NA_real_, NA_real_), c(x, y))))
    }

    v <- stats::cov.wt(data)  # "norm" method
    shape <- v$cov
    center <- v$center
    chol_decomp <- chol(shape)
    radius <- sqrt(dfn * stats::qf(level, dfn, dfd))

    angles <- (0:segments) * 2 * pi / segments
    unit_circle <- cbind(cos(angles), sin(angles))
    ellipse <- t(center + radius * t(unit_circle %*% chol_decomp))

    ellipse <- as.data.frame(ellipse)
    colnames(ellipse) <- c(x, y)
    ellipse
}

# A function to add new points on an already existing PCA
# ARGS:
#   - pca_object: a PCA object computed with dudi.pca
#   - plot_obs: a ggplot containing the output of a previous PCA
#   - new_df: a data.frame containing the data to add to the existing PCA
#   - group: a string. If given, filters out all rows for which df$group == 0
# RETURNS: a list:
#   - sp_pca: the coordinates of the old PCA + new data
#   - sp_plot_ind: the PCA plot showing individual observations in 2D
#   - sp_coords: the coordinates of the new data
project_on_existing_pca <- function(pca_object, plot_obs, new_df, group = NULL) {
    if (inherits(new_df, "SpatVector")) {
        new_df <- as.data.frame(new_df)
    } else if (!inherits(new_df, "data.frame")) {
        stop(paste0("df must be a dataframe, got ", class(new_df), "."))
    }

    # select group
    if (!is.null(group)) {
        new_df_filtered <- new_df |> filter(.data[[group]] != 0)
    } else {
        group <- "single species"
        new_df_filtered <- new_df
    }

    # clean the dataset in the same way that make_pca_and_plots does it
    new_df_clean <- new_df_filtered |> drop_na()
    pca_vars <- names(pca_object$tab)
    new_df_num <- new_df_clean |> select(all_of(pca_vars))
    if (length(names(new_df_num)) != length(pca_vars)) {
        stop("New dataframe given has different columns than those of pca.")
    }

    # predict projection for new dataset
    sup <- suprow(pca_object, new_df_num)
    new_pca_coords <- sup$lisup
    coords_df <- as.data.frame(new_pca_coords)[, 1:2]
    colnames(coords_df) <- c("Dim.1", "Dim.2")

    # Make plot manually
    p_ind <- plot_obs +
        geom_point(
            data = coords_df,
            aes(x = Dim.1, y = Dim.2),
            color = PALETTE[1]) +
        stat_ellipse(
            data = coords_df,
            aes(x = Dim.1, y = Dim.2, color = PALETTE[1], fill = PALETTE[1]),
            geom = "polygon", alpha = 0.15, type = "norm") +
        labs(caption = paste0("Colored group: ", group, ""))
    p_ind <- my_custom_ggplot_theme(p_ind, light = TRUE) +
        theme(legend.position = "none")

    return(list(
        sp_pca = new_pca_coords, 
        sp_coords = coords_df,
        sp_plot_ind = p_ind)) 
}

# A function that automatically selects numerical columns in a dataframe, 
# standardizes them and computes a UMAP
# ARGS:
#   - df: a dataframe (or SpatVector object)
#   - color_var: the name of a column of categorical values in df
#   - max_sample: maximum size of df (if nrow(df) > max_sample, runs umap on 
#       a subset of size=max_sample, then transforms the remaining data).
# RETURNS:
#   - umap: the UMPA result values
#   - p_ind: the UMAP plot showing individual observations in 2D
make_umap_and_plots <- function(
        df, color_var = NULL, max_sample = 500000, seed = 42) {
    # Make sure df is a dataframe
    if (inherits(df, "SpatVector")) {
        df <- as.data.frame(df)
    } else if (!inherits(df, "data.frame")) {
        stop(paste0("df must be a dataframe, got ", class(df), "."))
    }

    # Make sure NAs are excluded
    df_clean <- df |> drop_na()

    # Keep only numeric columns for PCA
    df_num <- df_clean |> select(where(is.numeric))

    # Standardize (always) - fit on full data so scaling is consistent
    df_scaled <- scale(df_num)

    n_total <- nrow(df_scaled)
    do_subsample <- n_total > max_sample

    if (do_subsample) {
        set.seed(seed)
        fit_idx <- sample(seq_len(n_total), size = max_sample)
    } else {
        fit_idx <- seq_len(n_total)
    }

    df_fit <- df_scaled[fit_idx, , drop = FALSE]

    # compute UMAP on the (sub)sample - always standardize data
    umap_res <- umap(
        df_fit,
        n_neighbors = 15,
        n_components = 2,
        metric = "euclidean",
        n_epochs = 200,
        min_dist = 0.1,
        spread = 1,
        init = "random",
        n_threads = parallel::detectCores() - 1,
        n_sgd_threads = parallel::detectCores() - 1,
        nn_method = "annoy",  # faster than NNDescent
        seed = seed,          # for reproducible result
        verbose = TRUE,
        ret_model = TRUE      # keeps the model in case you want to project new data later
    )

    if (do_subsample) {
        remaining_idx <- setdiff(seq_len(n_total), fit_idx)
        df_remaining <- df_scaled[remaining_idx, , drop = FALSE]

        message(sprintf(
            "Fit UMAP on %d sampled rows; projecting remaining %d rows...",
            length(fit_idx), length(remaining_idx)
        ))

        emb_remaining <- umap_transform(
            df_remaining,
            umap_res,
            n_threads = parallel::detectCores() - 1,
            n_epochs = 0,      # not applicable, ignored for transform
            verbose = TRUE
        )

        # Stitch fit + projected embeddings back into original row order
        emb_full <- matrix(NA_real_, nrow = n_total, ncol = 2)
        emb_full[fit_idx, ]       <- umap_res$embedding
        emb_full[remaining_idx, ] <- emb_remaining
        emb <- as.data.frame(emb_full)
    } else {
        emb <- as.data.frame(umap_res$embedding)
    }

    colnames(emb) <- c("UMAP1", "UMAP2")

    if (!is.null(color_var)) {
        emb$.color <- df_clean[[color_var]]
        p_ind <- ggplot(emb, aes(x = UMAP1, y = UMAP2)) +
            geom_point(aes(color = .color), alpha = 0.4, size = 0.3) +
            stat_ellipse(aes(color = .color), linewidth = 0.6) +
            scale_color_manual(values = ggsci::pal_jco()(length(unique(emb$.color))))
    } else {
        p_ind <- ggplot(emb, aes(x = UMAP1, y = UMAP2)) +
            geom_point(color = "grey20", alpha = 0.4, size = 0.3)
    }

    # Individuals plot (samples), colored by a grouping variable if provided
    p_ind <- p_ind +
        theme_minimal(base_size = 13) +
        labs(title = "UMAP - Individuals", color = color_var %||% "Group")
    p_ind <- my_custom_ggplot_theme(p_ind, with_palette = TRUE, light = TRUE)

    return(list(umap = umap_res, plot_ind = p_ind, embedding_full = emb))
}

##### Ecology functions #####
# A function that computes occurence rank curves from a dataframe of presence
# absence data.
# ARGS:
#   - ref_df: a data.frame containing presence-absence data (0s and 1s).
#   - ref_sp_names: the name of the columns of ref_df to use.
#   - list_dfs: a list containing data.frames (with presence-absence data).
#   - list_sp_names: the name of the columns of the data.frames in list_dfs to
#       use. All data.frames must have the same column names.
# RETURNS: list
#   - roc_plot: the ggplot.
#   - occ_freq: the occurence dataframe used to make the plot.
rank_occurence_curves <- function(
        ref_df, ref_sp_names, 
        list_dfs = NULL, list_sp_names = NULL) {
    
    # filter out NAs
    clean_ref_df <- ref_df |> drop_na()
    
    # Compute occurences
    ref_occ_freq <- data.frame(
            species = names(clean_ref_df[[ref_sp_names]]),
            n_detections = colSums(clean_ref_df[[ref_sp_names]]),
            n_total = nrow(clean_ref_df[[ref_sp_names]])) |>
        mutate(freq_occurrence = n_detections / n_total) |>
        arrange(desc(freq_occurrence)) |>
        mutate(rank = row_number())

    # Empty base plot
    roc_plot <- ggplot(ref_occ_freq, aes(x = rank, y = freq_occurrence)) +
        labs(x = "Species rank",
             y = "Frequency of occurrence (Sightings / All observations)",
             title = "Rank-Occurrence Curve")

    # Add simulation curves (if there was some)
    if (!is.null(list_dfs) && !is.null(list_sp_names)) {
        if (length(list_dfs) > length(PALETTE)) {
            # quick check because i know myself, I won't do it right each time
            stop(paste0(
                "PALETTE has only ", length(PALETTE), 
                " colors, list_dfs needs ", length(list_dfs),"."))
        }

        list_occ_freq_tables <- list()

        for (local_df in list_dfs) {
            clean_local_df <- local_df |> drop_na()

            # compute the occurence table for each dataframe in list_dfs
            list_occ_freq_tables[[length(list_occ_freq_tables)+1]] <- data.frame(
                    species = names(clean_local_df[[list_sp_names]]),
                    n_detections = colSums(clean_local_df[[list_sp_names]]),
                    n_total = nrow(clean_local_df[[list_sp_names]])) |>
                mutate(freq_occurrence = n_detections / n_total) |>
                arrange(desc(freq_occurrence)) |>
                mutate(rank = row_number())

            # add line to original plot
            roc_plot <- roc_plot + 
                geom_line(
                    data = list_occ_freq_tables[[length(list_occ_freq_tables)]],
                    aes(x = rank, y = freq_occurrence),
                    color = PALETTE[length(list_occ_freq_tables)])
        }
        roc_plot <- roc_plot + 
            labs(caption = "Colours: black = reference, rainbow = simulations.")

    } else if (!is.null(list_dfs)) {
        stop("'list_dfs' was given but 'list_sp_names' is missing.")
    } else if (!is.null(list_sp_names)) {
        stop("'list_sp_names' was given but 'list_dfs' is missing.")
    }
    
    # add black reference last (foreground)
    roc_plot <- roc_plot +
            geom_line() +
            geom_point(size = 2)
    roc_plot <- my_custom_ggplot_theme(roc_plot, with_palette = FALSE)
    
    return(list(roc_plot = roc_plot, occ_table = ref_occ_freq))
}
