# Set of functions used to plot harmonised figures

##### Liraries ##### ----------------------------------------------------------
library(colorspace)
library(dotwhisker)
library(tidyterra)
library(patchwork)
library(reshape2)
library(lmerTest)
library(MCMCvis)
library(ggplot2)
library(GGally)
library(readr)
library(dplyr)
library(rlang)
library(terra)
library(Hmsc)
library(coda)
library(lme4)
library(cli)
library(sf)

options(bitmapType = "cairo")
source(here::here("R/utils_data.R")) 

##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # all parameters are grouped together
source(here::here("R/utils_data.R"))  # needs get_metropolitan_france_shapefile

##### Global functions ##### --------------------------------------------------
# A wrapper to create a custom ggplot theme 
# ARGS:
#   - figure: a ggplot object.
#   - with_palette: a boolean. 
#       If TRUE (default) uses custom color/fill/shape/sizes.
#   - light: a boolean.
#       If TRUE (default is FALSE) when with_palette is TRUE,
#       uses the same color palette for color and fill values.
# RETURNS:
#   - figure, with customized theme and colors.
my_custom_ggplot_theme <- function(figure, with_palette=TRUE, light=FALSE){
    customised_fig <- figure +
        theme_linedraw(
            base_family = FONT
        ) +
        theme(
            # title and subtitle styling
            plot.title.position = "plot",
            plot.title = element_text(
                size = 18,
                face = "bold",
                color = "#000000",  
                margin = margin(b = 10)
            ),
            plot.subtitle = element_text(
                size = 14,
                color = "#777777", 
                margin = margin(b = 10)
            ),
            
            # plot styling
            plot.caption.position = "plot",
            plot.caption = element_text(
                size = 9,
                color = "#999999", 
                margin = margin(t = 15),
                hjust = 0
            ),
            axis.text = element_text(
                size = 11,
                color = "#000000"
            ),
            aspect.ratio = 1,
            
            # external grid
            axis.ticks = element_line(
                linetype = "solid",
                linewidth = 0.50,
                color = "#000000"
            ),
            panel.border = element_rect(
                colour = "#000000",
                linewidth = 1,
                fill = NA
            ),
            
            # internal grid
            panel.grid.major = element_line(
                linetype = "solid",
                linewidth = 0.15,
                color = "#999999"
            ),
            panel.grid.minor = element_blank(),
        )

        if (with_palette){
            if (light) {
                return(customised_fig + LIGHT_CUSTOM_SCALES)
            } else {
                return(customised_fig + CUSTOM_SCALES)
            }
        } else {
            return(customised_fig)
        }
}

# A function to save a ggplot figure to pdf
# ARGS:
#   - figure: a ggplot object.
#   - save_path: a filepath to create a pdf file.
# RETURNS:
#   - NULL, silently saves figure to a local file.
standardised_ggplot_save <- function(figure, save_path, .width = 18, .height = 6){
    # check if string ends with ".pdf"
    if (!endsWith(save_path, ".pdf")){
        stop(paste("Provided save_path must end with '.pdf', got", save_path))
    }
  
    ggsave(
        filename = save_path,
        plot = figure,
        dpi = 300,
        width = .width,             # large width to account for plots with very wide legends
        height = .height,
        device = cairo_pdf)
}

# A function that creates basic histograms with ggplots
# ARGS:
#   - df: a data.frame with columns x and category.
#   - x: a string. The name of a column to compute histograms on.
#   - category: a string. The name of a column to color the histograms with (default is NULL, no color added).
#   - bins: A numeric. Controls the number of bins (default is 1).
#   - breaks: A numeric. Default is NULL (not toggled) if given, overwrites bins.
# RETURNS:
#   - a ggplot object
ggplot_bars <- function(df, x, category = NULL, bins = 10, breaks = NULL, underlayers = NULL) {
  if (is.null(category)) {
    graph <- ggplot(data = df, aes(x = .data[[x]])) +
      underlayers +                          # drawn first, behind bars
      geom_histogram(
        color = darken(PALETTE[4], amount = 0.5),
        fill = PALETTE[4],
        bins = if (is.null(breaks)) bins else NULL,
        breaks = breaks)
  } else {
    graph <- ggplot(data = df, aes(x = .data[[x]], color = .data[[category]], fill = .data[[category]])) +
      underlayers +                          # drawn first, behind bars
      geom_histogram(binwidth = binwidth, position = "identity", alpha = 0.5)
  }
  return(my_custom_ggplot_theme(graph, with_palette = TRUE))
}


##### Maps functions ##### ----------------------------------------------------
# A function that creates a simple background map of france in ggplot2
# ARGS:
#   - borders_type: a string. Either "national" (default) or "regional". 
#       If "regional", draws highest level inner borders ofthe country.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in EPSG:4326).
# RETURNS:
#   - a ggplot figure
ggplot_get_france_base_map <- function(
        borders_type = "national",
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX
    ){
    europe_shp = rnaturalearth::ne_countries(
        continent = "Europe", scale = "large", returnclass = "sf")
    france_shp = get_metropolitan_france_shapefile(borders_type)
    
    base_map <- ggplot(europe_shp) +
        geom_sf(fill = "grey80", color = "white") +                     # color of countries
        geom_sf(data = france_shp, fill = "white", color = "black") +   # color of France
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "lightcyan1", color = NA),
            panel.grid.major = element_line(color = "lightcyan1"),
            panel.grid.minor = element_line(color = "lightcyan1")
        ) +
        labs(
            x = "longitude",
            y = "latitude"
        ) +
        coord_sf(
            xlim = c(lon_min, lon_max), 
            ylim = c(lat_min, lat_max))
    
    return(my_custom_ggplot_theme(base_map, with_palette=FALSE))
}

# A function to create a ggplot that shows a shapefile of categorical values 
# on a given map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - color_df: a dataframe. Contains columns "Value", "hex" 
# RETURNS:
#   - a ggplot figure
ggplot_categorical_shapefile_on_background_map <- function(
        background_map,
        shapefile,
        layer_name,
        color_df,
        label_layer_name = "NEW_LABEL3"
    ) {

    shapefile[[layer_name]] <- as.factor(shapefile[[layer_name]])
    color_key <- setNames(color_df$hex, color_df$Value)
    label_key <- setNames(color_df[[label_layer_name]], color_df$Value)
    
    # make ggplot
    map_category_grid <- suppressMessages(background_map +
        geom_sf( # color field to hide the anti-aliasing between polygons
            data = shapefile, 
            alpha = 1,
            aes(fill = .data[[layer_name]], color = after_scale(fill)),
            linewidth = 0.1) +
        scale_fill_manual(
            values = color_key,
            labels = label_key,
            na.value = "transparent") +
        scale_color_manual(
            values = color_key,
            labels = label_key,
            na.value = "transparent") +
        labs(
            x = "longitude", 
            y = "latitude", 
            fill = "Land Cover", color = "Land Cover") +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )
    
    return(my_custom_ggplot_theme(map_category_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows a shapefile of continuous values on a map
# ARGS:
#   - background_map: a ggplot object. 
#       The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer to show.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). 
#       Sets hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale 
#       (values are rounded to closest precision_auto_limits)
# RETURNS:
#   - a ggplot figure
ggplot_quantitative_shapefile_on_background_map <- function(
        background_map,
        shapefile,
        layer_name,
        unit="°C",
        limits=NULL,
        precision_auto_limits = 1e-5
    ) {

    if (is.vector(limits) && length(limits) == 2){
        low_limit <- limits[1]
        high_limit <- limits[2]
    } else if (is.null(limits)) {
        # round palette scale to the bottom and top nearest multiple of 5
        low_limit <- floor(min(shapefile[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
        high_limit <- ceiling(max(shapefile[[layer_name]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
    } else {
        stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got", limits))
    }

    map_quantity_grid <- suppressMessages(background_map +
        geom_sf( # color field to hide the anti-aliasing between polygons
            data = shapefile, 
            alpha = 1,
            aes(fill = .data[[layer_name]], color = after_scale(fill)),
            linewidth = 0.1) +
        scale_fill_continuous(
            na.value = "transparent", 
            palette = "turbo",
            limits = c(low = low_limit, high = high_limit)) +
        scale_color_continuous(
            na.value = "transparent", 
            palette = "turbo",
            limits = c(low = low_limit, high = high_limit)) +
        labs(x = "longitude", y = "latitude", fill = unit, color = unit) +
        coord_sf(
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )

    return(my_custom_ggplot_theme(map_quantity_grid, with_palette = FALSE))
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - lon_c: a string. The name of the column with longitude values.
#   - lat_c: a string. The name of the column with latitude values.
#   - column: a string (optional). The name of a column with categorical values.
#   - col, size, shape, stroke : default markers for ggplot.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in EPSG:4326).
# RETURNS:
#   - a ggplot figure
ggplot_categorical_df_on_background_map <- function(
        background_map, 
        df, 
        lon_c = "LON",
        lat_c = "LAT",
        column = NULL,
        legend_title = "Sampling location",
        col = PALETTE[1], size = SIZES[1], shape = SHAPES[1], stroke = STROKES[1],
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
       
    # make sure that coordinates are in the right coordinates system
    data <- st_as_sf(df, coords = c(lon_c, lat_c), crs = 4326)
    
    # shuffle to avoid biased overlaps
    data <- slice_sample(data, prop = 1)
    
    if (!is.null(column)) {
        uniques_vals <- unique(data[[column]])
        if (length(uniques_vals) > length(palette)){
            stop(paste(c(
                "Given palette can handle up to ", length(palette), 
                " unique values. Found ", length(uniques_vals),
                " in column '", column,"'."
            )))
        } else {
            # make plot
            map_obs <- suppressMessages(background_map +
                geom_sf(
                    data = data,
                    stroke = 0.8,
                    aes(
                        color = .data[[column]],
                        fill = .data[[column]],
                        shape = .data[[column]],
                        size = .data[[column]],
                    )
                ) +
                labs(
                    x = "longitude",
                    y = "latitude",
                    color = legend_title,
                    fill = legend_title,
                    shape = legend_title,
                    size = legend_title
                ) +
                coord_sf(
                    xlim = c(LON_MIN, LON_MAX),
                    ylim = c(LAT_MIN, LAT_MAX))
            )
        }
        return(my_custom_ggplot_theme(map_obs, with_palette = TRUE))
        
    } else {
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                size = size,
                shape = shape,        
                fill = col, 
                color = darken(col, amount = 0.66), 
                stroke = stroke
            ) +
            labs(
                x = "longitude",
                y = "latitude"
            ) +
            coord_sf(
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
    }
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - lon_c: a string. The name of the column with longitude values.
#   - lat_c: a string. The name of the column with latitude values.
#   - column: a string (optional). The name of a column with categorical values.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). 
#       Sets a hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale 
#       (values are rounded to closest precision_auto_limits)
#   - col, size, shape, stroke : default markers for ggplot.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in EPSG:4326).
# RETURNS:
#   - a ggplot figure
ggplot_quantitative_df_on_background_map <- function(
        background_map, 
        df, 
        lon_c = "LON",
        lat_c = "LAT",
        column = NULL,
        unit = NULL,
        limits = NULL,
        precision_auto_limits = 1e-5,
        col = PALETTE[1], size = SIZES[1], shape = SHAPES[1], stroke = STROKES[1],
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
       
    # make sure that coordinates are in the right coordinates system
    data <- st_as_sf(df, coords = c(lon_c, lat_c), crs = 4326)
    
    # shuffle to avoid biased overlaps
    data <- slice_sample(data, prop = 1)

   
    if (!is.null(column)) {
        if (is.vector(limits) && length(limits) == 2){
            low_limit <- limits[1]
            high_limit <- limits[2]
        } else if (is.null(limits)) {
            # round palette scale to the bottom and top nearest multiple of 5
            low_limit <- floor(min(df[[column]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
            high_limit <- ceiling(max(df[[column]], na.rm = TRUE) / precision_auto_limits) * precision_auto_limits
        } else {
            stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got", limits))
        }

        
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                stroke = 0.8,
                aes( 
                    fill = .data[[column]], 
                    color = .data[[column]]
                )
            ) +
            scale_fill_continuous(
                na.value = "transparent", 
                palette = "turbo",
                limits = c(low = low_limit, high = high_limit)) +
            scale_color_continuous(
                na.value = "transparent", 
                palette = "turbo",
                limits = c(low = low_limit, high = high_limit)) +
            labs(
                x = "longitude",
                y = "latitude",
                fill = unit,
                color = unit) +
            coord_sf(
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
        
    } else {
        # make plot
        map_obs <- suppressMessages(background_map +
            geom_sf(
                data = data,
                size = size,
                shape = shape,        
                fill = col, 
                color = darken(col, amount = 0.66), 
                stroke = stroke
            ) +
            labs(
                x = "longitude",
                y = "latitude"
            ) +
            coord_sf(
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        return(my_custom_ggplot_theme(map_obs, with_palette = FALSE))
    }
}


##### HMSC interpretation ##### -----------------------------------------------
# A function that mimicks Hmsc::plotBeta
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - post: post posterior summary of Beta parameters obtained from getPostEstimate().
#   - supportLevel: a numeric threshold for plotting, values between 0.5 and 1 (default 0.95).
# RETURNS:
#   - a ggplot object.
ggplot_custom_plotBeta <- function(hM, post, supportLevel = 0.95) {
  
    # Reproduce the Support calculation from source
    betaP  <- post$support
    toPlot <- 2 * betaP - 1
    toPlot <- toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
    betaMat <- matrix(toPlot, nrow = hM$nc, ncol = ncol(hM$Y))

    rownames(betaMat) <- hM$covNames
    colnames(betaMat) <- hM$spNames

    # Long format for ggplot
    df <- as.data.frame(as.table(betaMat))
    colnames(df) <- c("Covariate", "Species", "value")

    plot <- ggplot(df, aes(x = Covariate, y = Species, fill = value)) +
        geom_tile(color = "grey90") +
        scale_fill_gradient2(
            low     = PALETTE[3],
            mid     = "white",
            high    = PALETTE[1],
            midpoint = 0,
            limits  = c(-1, 1),
            name    = "Support"
        ) +
        theme_minimal() +
        theme(
            axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
            axis.text.y  = element_text(face = "italic"),
            panel.grid   = element_blank()
        ) +
        labs(
            x = NULL, y = NULL, 
            subtitle = paste0("Showing support levels >=", supportLevel, ".")) 

    return(my_custom_ggplot_theme(plot, with_palette = FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1)))
}

# A function that mimicks Hmsc::computeAssociations + Corplot
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - supportLevel: a numeric threshold for plotting, values between 0.5 and 1 (default 0.95).
# RETURNS:
#   - a ggplot object.
ggplot_custom_random_corr_associations <- function(hM, supportLevel = 0.95) {
    OmegaCor = computeAssociations(hM)
    toPlot = ((OmegaCor[[1]]$support>supportLevel)
        + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean

    # Convert matrix to long format for ggplot
    toPlot_df <- melt(toPlot)
    colnames(toPlot_df) <- c("Var1", "Var2", "value")

    plot <- ggplot(toPlot_df, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
            low  = "blue",
            mid  = "white",
            high = "red",
            midpoint = 0,
            limits = c(-1, 1),
            name = "Correlation"
        ) +
        scale_y_discrete(limits = rev(levels(factor(toPlot_df$Var2)))) + 
        labs(subtitle = paste("random effect level:", fitted_model$rLNames[1])) +
        theme_minimal() +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y  = element_text(size = 8),
            axis.title   = element_blank(),
            plot.title   = element_text(hjust = 0.5),
            panel.grid   = element_blank()
        ) +
        coord_fixed()

    return(my_custom_ggplot_theme(plot, with_palette = FALSE) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1)))
}

# A function that mimicks Hmsc::plotVariancePartitioning
# ARGS:
#   - hM: a fitted Hmsc model object.
#   - VP: a matrix obtained from Hmsc::computeVariancePartitioning.
# RETURNS:
#   - a ggplot object.
ggplot_custom_plotVariancePartitioning <- function(hM, VP) {

    # Build labels with means
    if (!(length(fitted_model$rLNames) == 0)) {
        labels <- c(VP$groupnames, paste0("Random: ", hM$rLNames))
    } else {
        labels <- c(VP$groupnames)
    }
    
    means <- round(100 * rowMeans(VP$vals), 1)
    labels <- paste0(labels, " (mean = ", means, ")")

    # Long format
    df <- as.data.frame(VP$vals)
    df$Group <- factor(labels, levels = rev(labels))
    df_long <- pivot_longer(df, -Group, names_to = "Species", values_to = "Proportion")

    plot <- ggplot(df_long, aes(x = Species, y = Proportion, fill = Group)) +
        geom_col() +
        labs(x = "Species", y = "Variance proportion") 

    return(my_custom_ggplot_theme(plot, with_palette = TRUE) +
        theme(axis.text.x  = element_text(angle = 45, hjust = 1)))
}

# A function to compare scores between k_folds, subset and model type.
# ARGS:
#   - parent_folder: a string. 
#       Path to parent of subfolders containing `{subset_name}_scores.csv`.
#   - reference_model_combination: a list of single parameters. 
#   - loop_model_combination: a list of parameters (single, 
#       except one parameter, that is a vector of several elements).
#   - loop_on: a string. The name of a parameter in a combination. 
#       The parameter which has several values.
#   - k_fold: a numeric. The number of cross-validation subsets to make. 
#   - metric: a string. Metric to extract from file (MSE, RMSE, AUC or TjuR2).
#   - subset_names: a list string. Usually c("train", "val", "test").
#   - xlabel: a string. The xlabel for the plot (default is "Effect").
#   - group_species: whether to take the mean 
#       accross all k_folds and species (TRUE, default) or 
#       only accross k_folds (FALSE).
#   - species_names: names of species in CSV 
#       (rownames are not available from csvs).
#   - save_to: a string. Path which should end with .pdf. 
#       If NULL, does not save pdf.
boxplot_compare_scores <- function(
        parent_folder,
        reference_model_combination,
        loop_model_combination,
        loop_on,
        k_fold,
        metric = "MSE",
        subset_names = c("train", "val", "test"),
        xlabel = "Effect",
        group_species = TRUE,
        species_names = NULL,
        save_to = NULL) {

    # auto compute differences between reference scores and the other scores
    diffs <- compute_score_diffs(
        parent_folder, reference_model_combination, loop_model_combination,
        loop_on, k_fold, metric, subset_names, group_species, species_names)
    scores_df <- diffs$summary
    raw_diffs_df <- diffs$raw_diffs

    # load reference values for display in captions
    refs_list <- list()
    for (subset_name in subset_names) {
        for (k in 1:k_fold) {
            run_path <- make_run_path(
                parent_folder, reference_model_combination, k)
            refs <- load_metric_scores(run_path, subset_name, metric)

            refs_list[[length(refs_list) + 1]] <- data.frame(
                subset = subset_name,
                k_fold = k,
                scores = refs,
                species = species_names,
                stringsAsFactors = FALSE
            )
        }
    }
    refs_df <- bind_rows(refs_list)
    if (group_species) {
        mean_refs_df <- refs_df |>
            group_by(subset) |>
            summarise(mean_score = mean(scores, na.rm = TRUE), .groups = "drop")
    } else {
        mean_refs_df <- refs_df |>
            group_by(subset, species) |>
            summarise(mean_score = mean(scores, na.rm = TRUE), .groups = "drop")
    }

    # just to be sure
    if (!group_species) {
        scores_df <- scores_df |> 
            mutate(species = factor(species, levels = species_names))
        raw_diffs_df <- raw_diffs_df |> 
            mutate(species = factor(species, levels = species_names))
    }

    if (group_species) {
        bottom_caption <- paste(
            "Distribution of per-fold differences",
            "(species-averaged within each fold, across k_fold).")
    } else {
        bottom_caption <- paste(
            "Distribution of per-fold differences across k_fold, per species.")
    }

    # format formula to number of variables
    temp_loop <- loop_model_combination
    temp_loop$HMSC_XFORMULAS <- lapply(lapply(
        temp_loop$HMSC_XFORMULAS, all.vars), length)
    temp_ref <- reference_model_combination
    temp_ref$HMSC_XFORMULAS <- lapply(lapply(
        temp_ref$HMSC_XFORMULAS, all.vars), length)
    
    # add model type to caption
    ref_type <- ""
    model_type <- ""
    for (param in names(temp_loop)) {
        if (length(temp_loop[[param]]) > 1) {
            model_type <- paste0(
                model_type, tolower(param), "=see x-axis, ")
        } else {
            model_type <- paste0(
                model_type, tolower(param), "=", temp_loop[[param]], ", ")
        }
        ref_type <- paste0(
                ref_type, tolower(param), "=", temp_ref[[param]], ", ")
    }
    model_type <- paste0(substr(model_type, 1, nchar(model_type)-2), ".")
    ref_type <- paste0(substr(ref_type, 1, nchar(ref_type)-2), ".")
    bottom_caption <- paste0(
        bottom_caption,
        "\nReference: ", ref_type, 
        ".\nCompared with: ", model_type)
    
    p <- ggplot(
            raw_diffs_df, 
            aes(y = diff_value, x = loop_element, fill = subset)) +
        geom_boxplot(
            position = position_dodge(width = 0.75), 
            width = 0.6, 
            outlier.shape = 16) +
        labs(caption = bottom_caption, fill = "Subset") +
        ylab(paste("Delta in average", metric)) +
        xlab(xlabel) +
        geom_hline(yintercept = 0, linetype = "dashed")

    # add text for mean reference
    if (group_species) {
        ref_values <- paste0(mean_refs_df |> 
            mutate(subset = factor(subset, levels = subset_names)) |>
            arrange(subset) |>
            mutate(label = paste0(subset, ": ", round(mean_score, 3))) |>
            pull(label), collapse=", ")

        p <- p + annotation_custom(
            grob = grid::textGrob(
                paste0("Reference scores: ", ref_values, "."),
                x = unit(0, "npc"), y = unit(0, "npc"),
                hjust = -0.02, vjust = -0.75,
                gp = grid::gpar(
                    fontsize = 9, fontface = "italic", lineheight = 0.8)
            )
        )
    }

    if (!group_species) p <- p + facet_wrap(~species, nrow = 1)

    p <- my_custom_ggplot_theme(p)  +
        scale_fill_manual(values = c(PALETTE[2], PALETTE[3], PALETTE[1])) +
        scale_x_discrete(labels = abbreviate_loop_labels)

    if (grepl("MSE", metric)) {
        p <- p + scale_y_reverse()
    }

    finalize_plot(p, save_to, what = "data")
    return(list(scores = scores_df, diffs = raw_diffs_df, plot = p))
}