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
#   - is_map: a boolean.
#       If TRUE, drops the fixed aspect ratio (default is FALSE).
# RETURNS:
#   - figure, with customized theme and colors.
my_custom_ggplot_theme <- function(
    figure, with_palette=TRUE, light=FALSE, is_map=FALSE){
    customised_fig <- figure +
        theme_linedraw(
            base_family = FONT
        ) +
        theme(
            aspect.ratio = if (is_map) NULL else 1,
            
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

# A function that abbreviates hyphenated labels of a ggplot. A label containing 
# a single word is left as-is, a multi-word (e.g. "foo-bar-baz") element 
# becomes initials (e.g. "FBB").
# ARGS:
#   - x: a vector of string.
# RETURNS:
#   - The transformed vector of strings
abbreviate_labels <- function(x) {
  sapply(strsplit(x, "-"), function(words) {
    if (length(words) == 1) {
      words
    } else {
      is_num <- grepl("^[0-9]+$", words)
      paste0(ifelse(is_num, words, toupper(substr(words, 1, 1))), collapse = "")
    }
  })
}

##### Maps functions ##### ----------------------------------------------------
# A function that creates a simple background map of france in ggplot2
# ARGS:
#   - borders_type: a string. Either "national" (default) or "regional". 
#       If "regional", draws highest level inner borders ofthe country.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
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
            crs = PROJECT_CRS,              # CRS the map is drawn in
            xlim = c(lon_min, lon_max), 
            ylim = c(lat_min, lat_max))
    
    return(my_custom_ggplot_theme(base_map, with_palette=FALSE, is_map = TRUE))
}

# A function to create a ggplot that shows a categorical SpatRaster
# on a given map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - raster: a SpatRaster. Categorical values (integer codes or factor).
#   - layer_name: a string or integer. The layer of the raster to show.
#       Defaults to the first layer.
#   - color_df: a dataframe. Contains columns "Value", "hex" and the label column.
#   - label_layer_name: a string. Column of color_df used for legend labels.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
#   - max_cells: a numeric, the maximum number of cells before downsampling.
# RETURNS:
#   - a ggplot figure
ggplot_categorical_raster_on_background_map <- function(
        background_map,
        raster,
        layer_name = "NEW_LABEL3",
        max_cells = 5e5,
        lon_min = LON_MIN, lon_max = LON_MAX,
        lat_min = LAT_MIN, lat_max = LAT_MAX) {

    # keep only the requested layer (avoids the "mixed layer classes" warning)
    raster <- raster[[layer_name]]

    # build label -> hex color from the embedded color table
    ct <- terra::coltab(raster)[[1]]
    lv <- terra::levels(raster)[[1]]
    m  <- merge(lv, ct, by.x = names(lv)[1], by.y = "value")
    color_key <- setNames(
        grDevices::rgb(m$red, m$green, m$blue, m$alpha, maxColorValue = 255),
        m[[layer_name]])

    # make ggplot
    map_category_grid <- suppressMessages(
        background_map +
            tidyterra::geom_spatraster(
                data = raster,
                aes(fill = .data[[layer_name]]),
                maxcell = max_cells) +
            scale_fill_manual(
                values = color_key,
                na.value = "transparent",
                na.translate = FALSE) +
            labs(
                x = "longitude",
                y = "latitude",
                fill = "Land Cover") +
            coord_sf(
                crs = PROJECT_CRS,
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
    )

    my_custom_ggplot_theme(
        map_category_grid, with_palette = FALSE, is_map = TRUE)
}

# A function to create a ggplot that shows a shapefile of categorical values 
# on a given map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer containing the values to show.
#   - color_df: a dataframe. Contains columns "Value", "hex" 
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
# RETURNS:
#   - a ggplot figure
ggplot_categorical_shapefile_on_background_map <- function(
        background_map,
        shapefile,
        layer_name,
        color_df,
        label_layer_name = "NEW_LABEL3",
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX
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
            crs = PROJECT_CRS,              # CRS the map is drawn in
            xlim = c(lon_min, lon_max), 
            ylim = c(lat_min, lat_max))
    )
    
    return(my_custom_ggplot_theme(
        map_category_grid, with_palette = FALSE, is_map = TRUE))
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - lon_c: a string. The name of the column with longitude values (EPSG:4326).
#   - lat_c: a string. The name of the column with latitude values (EPSG:4326).
#   - column: a string (optional). The name of a column with categorical values.
#   - col, size, shape, stroke : default markers for ggplot.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
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
    
    # Check that the dataframe is in EPSG:4326 coordinates (approx.)
    if (!(all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[lon_c])) &
        all(as.data.frame(df[lon_c]) <= FRANCE_BBOX_4326["xmax"])   &
        all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[lat_c]))   &
        all(as.data.frame(df[lat_c]) <= FRANCE_BBOX_4326["ymax"])    
    )) {
        stop(paste0(
            "ggplot_categorical_df_on_background_map() assumes that your data ",
            "is in EPSG:4326 over France's extent. It looks like 'df' is in ",
            "another CRS, or on a different area. Please check your data."
        ))
    }
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
                    crs = PROJECT_CRS,              # CRS the map is drawn in
                    xlim = c(lon_min, lon_max),
                    ylim = c(lat_min, lat_max))
            )
        }
        return(my_custom_ggplot_theme(
            map_obs, with_palette = TRUE, is_map = TRUE))
        
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
                crs = PROJECT_CRS,              # CRS the map is drawn in
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        return(my_custom_ggplot_theme(
            map_obs, with_palette = FALSE, is_map = TRUE))
    }
}

# A function to create a ggplot that shows a SpatRaster of continuous values on a map
# ARGS:
#   - background_map: a ggplot object.
#       The background map that will be used.
#   - raster: a SpatRaster.
#   - layer_name: a string or integer. The layer to show (defaults to the first).
#   - unit: a string. A label that will be shown along the palette displayed.
#   - transform_scale: a string. Transforms only the color scale.
#   - limits: a vector of 2 values (optional).
#       Sets hard limits on the values considered by the palette.
#   - cols: a list of custom colors.
#   - precision_auto_limits: when limits is NULL, precision of color scale
#       (values are rounded to closest precision_auto_limits)
#   - max_cells: a numeric, the maximum number of cells before downsampling.
# RETURNS:
#   - a ggplot figure
ggplot_quantitative_raster_on_background_map <- function(
        background_map,
        raster,
        layer_name = 1,
        unit = "°C",
        transform_scale = "sqrt",
        limits = NULL,
        cols = CONTRASTED_PLASMA_COLOURS,
        precision_auto_limits = 1e-5,
        max_cells = 5e5) {

    r <- raster[[layer_name]]
    names(r) <- "value"

    # limits from the full-resolution data, before aggregating,
    # so that aggregation doesn't shrink the extremes
    if (is.vector(limits) && length(limits) == 2) {
        low_limit  <- limits[1]
        high_limit <- limits[2]
    } else if (is.null(limits)) {
        mm <- terra::global(r, c("min", "max"), na.rm = TRUE)
        low_limit  <- floor(mm$min / precision_auto_limits) * precision_auto_limits
        high_limit <- ceiling(mm$max / precision_auto_limits) * precision_auto_limits
    } else {
        stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got",
                   paste(limits, collapse = ", ")))
    }

    # 2. downsample if there are more cells than needed
    if (terra::ncell(r) > max_cells) {
        fact <- ceiling(sqrt(terra::ncell(r) / max_cells))
        r <- terra::aggregate(r, fact = fact, fun = "mean", na.rm = TRUE)
    }  

    # raster -> data.frame of cell centers (drops NA cells)
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)

    if (is.vector(limits) && length(limits) == 2) {
        low_limit <- limits[1]
        high_limit <- limits[2]
    } else if (is.null(limits)) {
        low_limit <- floor(min(df$value) / precision_auto_limits) * precision_auto_limits
        high_limit <- ceiling(max(df$value) / precision_auto_limits) * precision_auto_limits
    } else {
        stop(paste("'limits' is not recognised. Expected vector of length 2 or NULL, got",
                   paste(limits, collapse = ", ")))
    }

    map_quantity_grid <- suppressMessages(
        background_map +
            geom_raster(data = df, aes(x = x, y = y, fill = value)) +
            scale_fill_gradientn(
                colours = cols,
                trans = transform_scale,
                limits = c(low_limit, high_limit),
                name = unit,
                na.value = "transparent"
            ) +
            labs(x = "longitude", y = "latitude", fill = unit) +
            coord_sf(
                crs = PROJECT_CRS,
                xlim = c(LON_MIN, LON_MAX),
                ylim = c(LAT_MIN, LAT_MAX)
            )
    )

    return(my_custom_ggplot_theme(
        map_quantity_grid, with_palette = FALSE, is_map = TRUE))
}

# A function to create a ggplot that shows a shapefile of continuous values on a map
# ARGS:
#   - background_map: a ggplot object. 
#       The background map that will be used.
#   - shapefile: a shapefile.
#   - layer_name: a string. The name of the layer to show.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - transform_scale: a string. Transforms only the color scale.
#   - limits: a vector of 2 values (optional). 
#       Sets hard limits on the values considered by the palette.
#   - cols: a list of custom colors.
#   - precision_auto_limits: when limits is NULL, precision of color scale 
#       (values are rounded to closest precision_auto_limits)
# RETURNS:
#   - a ggplot figure
ggplot_quantitative_shapefile_on_background_map <- function(
        background_map,
        shapefile,
        layer_name,
        unit = "°C",
        transform_scale = "sqrt",
        limits = NULL,
        cols = CONTRASTED_PLASMA_COLOURS,
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
        scale_fill_gradientn(
            colours = cols,
            trans = transform_scale,
            limits = c(low_limit, high_limit),
            name = unit,
            na.value = "transparent") +
        labs(x = "longitude", y = "latitude", fill = unit, color = unit) +
        coord_sf(
            crs = PROJECT_CRS,
            xlim = c(LON_MIN, LON_MAX), 
            ylim = c(LAT_MIN, LAT_MAX))
    )

    return(my_custom_ggplot_theme(
        map_quantity_grid, with_palette = FALSE, is_map = TRUE))
}

# A function to create a ggplot that shows scattered locations on a map
# ARGS:
#   - background_map: a ggplot object. The background map that will be used.
#   - df: a dataframe.
#   - lon_c: a string. The name of the column with longitude values (EPSG:4326).
#   - lat_c: a string. The name of the column with latitude values (EPSG:4326).
#   - column: a string (optional). The name of a column with categorical values.
#   - unit: a string. A label that will be shown along the palette displayed.
#   - limits: a vector of 2 values (optional). 
#       Sets a hard limits on the values considered by the palette.
#   - precision_auto_limits: when limits is NULL, precision of color scale 
#       (values are rounded to closest precision_auto_limits)
#   - cols: a list of custom colors.
#   - transform_scale: a string. Transforms only the color scale.
#   - col, size, shape, stroke : default markers for ggplot.
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS6).
# RETURNS:
#   - a ggplot figure
ggplot_quantitative_df_on_background_map <- function(
        background_map, 
        df, 
        lon_c = "LON",
        lat_c = "LAT",
        column = NULL,
        facet_formula = NULL,
        unit = NULL,
        limits = NULL,
        precision_auto_limits = 1e-5,
        cols = CONTRASTED_PLASMA_COLOURS,
        transform_scale = "sqrt",
        col = PALETTE[1], size = SIZES[1], shape = SHAPES[1], stroke = STROKES[1],
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {

    # Check that the dataframe is in EPSG:4326 coordinates (approx.)
    if (!(all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[lon_c])) &
        all(as.data.frame(df[lon_c]) <= FRANCE_BBOX_4326["xmax"])   &
        all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[lat_c]))   &
        all(as.data.frame(df[lat_c]) <= FRANCE_BBOX_4326["ymax"])    
    )) {
        stop(paste0(
            "ggplot_categorical_df_on_background_map() assumes that your data ",
            "is in EPSG:4326 over France's extent. It looks like 'df' is in ",
            "another CRS, or on a different area. Please check your data."
        ))
    }
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
            scale_fill_gradientn(
                colours = cols,
                trans = transform_scale,
                limits = c(low_limit, high_limit),
                name = unit,
                na.value = "transparent") +
            scale_color_gradientn(
                colours = cols,
                trans = transform_scale,
                limits = c(low_limit, high_limit),
                name = unit,
                na.value = "transparent") +
            labs(
                x = "longitude",
                y = "latitude",
                fill = unit,
                color = unit) +
            coord_sf(
                crs = PROJECT_CRS,              # CRS the map is drawn in
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        
        if (!is.null(facet_formula)) {
            map_obs <- map_obs + facet_grid(facet_formula)
        }

        return(my_custom_ggplot_theme(
            map_obs, with_palette = FALSE, is_map = TRUE))
        
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
                crs = PROJECT_CRS,              # CRS the map is drawn in
                xlim = c(lon_min, lon_max),
                ylim = c(lat_min, lat_max))
        )
        return(my_custom_ggplot_theme(
            map_obs, with_palette = FALSE, is_map = TRUE))
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
#   - ref_scores: a data.frame or tibble. 
#       Must contain columns: 
#           name (with metrics' names), 
#           score (with metrics' values),
#           sp (a list of species' names),
#           k (for each k_fold).
#   - compare_scores: a data.frame or tibble.
#       Must contain columns: 
#           name (with metrics' names), 
#           score (with metrics' values),
#           sp (a list of species' names),
#           k (for each k_fold),
#           x_var (see details below),
#           panel_var (see details below).
#   - diffs_scores: a data.frame or tibble.
#       Default is NULL (computed from ref_scores and compare_scores).
#       When given, skips computation and uses it directly.
#   - metric: a string. The name of a category in column "name".
#   - x_var: a string. 
#       The name of column in compare_scores, will be displayed on the x-axis.
#   - x_order: a vector of string. Must contains all unique values in 
#       compare_scores[[x_order]]. The order of the panels displayed.
#   - panel_var: a string. 
#       The name of column in compare_scores, will be displayed in different 
#       panels on a 1-row grid of plots.
#   - panel_order: a vector of string. Must contains all unique values in 
#       compare_scores[[panel_order]]. The order of the panels displayed.
#   - group_species: a boolean. When TRUE (default) averages scores over all
#       species. When FALSE, computes each species average over k_folds and
#       displays a grid of panels with species in column and panel_var in row.
#   - save_to: a string. Filepath to save the plot to. Must end with ".pdf".
compare_plot <- function(
        ref_scores = NULL, compare_scores = NULL, metric, diffs_scores = NULL,
        x_var = NULL, x_order = NULL,
        fill_var = NULL, fill_order = NULL,
        color_order = NULL, random_strategy = "random-sampling",
        panel_var = NULL, panel_order = NULL,
        group_species = TRUE,
        save_to = NULL) {

    if (is.null(diffs_scores)) {
        # species is kept as a key only when not averaging over it
        sp_key    <- if (group_species) NULL else "sp"
        join_keys <- c("k", panel_var, sp_key)

        # 1. Reference: mean score per fold (and panel, and species if kept)
        ref_summary <- ref_scores |>
            filter(name == metric) |>
            group_by(across(all_of(join_keys))) |>
            summarise(ref_score = mean(score, na.rm = TRUE), .groups = "drop")

        # 2. Compared scores: mean per fold, x, fill, panel (and species)
        compare_summary <- compare_scores |>
            filter(name == metric) |>
            group_by(across(all_of(c(join_keys, x_var, fill_var)))) |>
            summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop")

        # 3. Per-fold difference with the matching reference
        diffs_df <- compare_summary |>
            left_join(ref_summary, by = join_keys) |>
            mutate(diff_score = avg_score - ref_score)
    } else {
        diffs_df <- diffs_scores
    }

    # dummy x so the violin still has a discrete position
    if (is.null(x_var)) {
        diffs_df$.x <- factor("all")
        x_col <- ".x"
    } else {
        diffs_df[[x_var]] <- if (is.null(x_order)) {
            factor(diffs_df[[x_var]])
        } else {
            factor(diffs_df[[x_var]], levels = x_order)
        }
        x_col <- x_var
    }

    if (!is.null(panel_var) && !is.null(panel_order)) {
        diffs_df <- diffs_df |> filter(.data[[panel_var]] %in% panel_order)
    }

    # dummy fill so the violin still has a single group to colour
    if (is.null(fill_var)) {
        diffs_df$.fill <- factor("all")
        fill_col <- ".fill"
    } else {
        diffs_df[[fill_var]] <- if (is.null(fill_order)) {
            factor(diffs_df[[fill_var]])
        } else {
            factor(diffs_df[[fill_var]], levels = fill_order)
        }
        fill_col <- fill_var
    }

    # ordered panel
    if (!is.null(panel_var)) {
        diffs_df[[panel_var]] <- if (is.null(panel_order)) {
            factor(diffs_df[[panel_var]])
        } else {
            factor(diffs_df[[panel_var]], levels = panel_order)
        }
    }

    # 4. Caption, with reference means per panel
    if (is.null(diffs_scores)) {
        ref_means <- ref_summary |>
            group_by(across(all_of(panel_var))) |>
            summarise(m = mean(ref_score, na.rm = TRUE), .groups = "drop")

        if (!is.null(panel_var)) {
            if (!is.null(panel_order)) {
                ref_means <- ref_means |>
                    mutate(across(all_of(panel_var), ~ factor(.x, levels = panel_order))) |>
                    arrange(across(all_of(panel_var)))
            }
            ref_text  <- paste0(ref_means[[panel_var]], ": ", round(ref_means$m, 3),
                                collapse = ", ")
            ref_label <- paste0(" (", panel_var, ")")
        } else {
            ref_text  <- as.character(round(ref_means$m, 3))
            ref_label <- ""
        }
    }

    if (group_species) {
        bottom_caption <- "Distribution of per-fold differences (species-averaged within each fold, across k_fold)."
    } else {
        bottom_caption <- "Distribution of per-fold differences across k_fold, per species."
    }
    if (is.null(diffs_scores)) {
        bottom_caption <- paste0(
            bottom_caption,
            "\nReference mean ", metric, ref_label, ": ", ref_text, "."
        )
        y_label <- "Delta in average "
    } else {
        y_label <- ""
    }

    # one summary value per violin: x value x fill value (x facet cells)
    summary_keys <- unique(c(x_col, fill_col, panel_var,
                             if (!group_species) "sp"))
    violin_summary <- diffs_df |>
        group_by(across(all_of(summary_keys))) |>
        summarise(m = mean(diff_score, na.rm = TRUE), .groups = "drop")

    # 5. Colours: tied to values (via color_order), not to plotting position
    lvls_now   <- levels(diffs_df[[fill_col]])
    color_lvls <- if (is.null(color_order)) {
        sort(lvls_now)
    } else {
        c(intersect(color_order, lvls_now), setdiff(lvls_now, color_order))
    }

    if (!is.null(random_strategy) && random_strategy %in% color_lvls) {
        # grey for the random strategy, palette for the others (as in distribution_plot)
        rest        <- setdiff(color_lvls, random_strategy)
        fill_colors <- c(setNames("grey70", random_strategy),
                         setNames(rep_len(PALETTE, length(rest)), rest))
    } else {
        fill_colors <- setNames(rep_len(PALETTE, length(color_lvls)), color_lvls)
    }

    fill_title <- if (is.null(fill_var)) {
        NULL
    } else {
        fill_var |> str_replace_all("_", " ") |> str_to_sentence()
    }

    # 6. Plot
    dodge <- position_dodge(width = 0.75)

    p <- ggplot(
            diffs_df,
            aes(y = diff_score, x = .data[[x_col]], fill = .data[[fill_col]])) +
        geom_violin(position = dodge, width = 0.6,
                    scale = "width", trim = TRUE,
                    colour = "grey30", alpha = 0.7, linewidth = 0.4) +
        geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
                   size = 0.8, alpha = 0.5, colour = "grey20", show.legend = FALSE) +
        geom_point(data = violin_summary,
                   aes(x = .data[[x_col]], y = m, group = .data[[fill_col]]),
                   inherit.aes = FALSE, position = dodge,
                   shape = 23, size = 2.5, fill = "white", colour = "black") +
        labs(caption = bottom_caption, fill = fill_title) +
        ylab(paste0(y_label, metric)) +
        xlab(if (is.null(x_var)) {
                NULL
            } else {
                x_var |> str_replace_all("_", " ") |> str_to_sentence()
            }
        ) +
        geom_hline(yintercept = 0, linetype = "dashed")

    # faceting
    if (group_species) {
        if (!is.null(panel_var)) {
            p <- p + facet_wrap(vars(.data[[panel_var]]), nrow = 1)
        }
    } else {
        p <- p + if (!is.null(panel_var)) {
            facet_grid(rows = vars(.data[[panel_var]]), cols = vars(sp))
        } else {
            facet_grid(cols = vars(sp))
        }
    }

    p <- my_custom_ggplot_theme(p) +
        scale_fill_manual(values = fill_colors)

    if (is.null(fill_var)) p <- p + guides(fill = "none")

    if (is.null(x_var)) {
        p <- p + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank())
    } else {
        p <- p + scale_x_discrete(labels = abbreviate_labels)
    }

    # lower is better for MSE/RMSE, so flip the axis
    if ((grepl("MSE", metric) || grepl("SD", metric)) && (!grepl("Δ", metric))) {
        p <- p + scale_y_reverse()
    }

    # Save and return results
    if (!is.null(save_to)) standardised_ggplot_save(p, save_to)
    list(diffs = diffs_df, plot = p)
}

# A function to compare scores between k_folds, subset and model type.
# ARGS:
#   - ref_scores: a data.frame or tibble. 
#       Must contain columns in group_vars, in panel_var, sp and metric.
#   - compare_scores: a data.frame or tibble.
#       Must contain columns in group_vars, in panel_var, sp and metric.
#   - metric: a string. The name of a category in column "name".
#   - diff_order: a booelan. Default is TRUE: lower values are better.
#   - threshold: a numeric between 0 and 1. The relative increase in metric
#       to qualify as improvement.
#   - group_vars: a vector of string. The name of columns that are constants 
#       in ref_df.
#   - panel_var: a vector of string. The name of columns whose values appear 
#       in both tibbles/dataframes.
count_improving <- function(
        ref_scores, compare_scores, metric,
        diff_order = 1,
        threshold = 0.01,
        group_vars = NULL,
        panel_var = NULL) {

    # Species always kept: improvement is assessed per species
    join_keys <- c("k", panel_var, "sp")

    # 1. Reference: mean score per fold (and panel) and species
    ref_summary <- ref_scores |>
        filter(name == metric) |>
        group_by(across(all_of(join_keys))) |>
        summarise(ref_score = mean(score, na.rm = TRUE), .groups = "drop")

    # 2. Compared: mean per fold, panel, species and grouping columns
    compare_summary <- compare_scores |>
        filter(name == metric) |>
        group_by(across(all_of(c(join_keys, group_vars)))) |>
        summarise(cmp_score = mean(score, na.rm = TRUE), .groups = "drop")

    # 3. Per-fold, per-species relative improvement (lower is better)
    if (!diff_order) {
        order_dir = -1
    } else {
        order_dir = 1
    }
    per_sp <- compare_summary |>
        left_join(ref_summary, by = join_keys) |>
        mutate(rel_improvement = order_dir * (ref_score - cmp_score) / ref_score,
               improved = rel_improvement >= threshold)

    # 4. Count species above threshold
    improved_sp <- per_sp |>
        group_by(across(all_of(c(group_vars, panel_var, "k")))) |>
        summarise(
            n_sp_above = sum(rel_improvement > threshold, na.rm = TRUE),
            n_sp_total = sum(!is.na(rel_improvement)),
            .groups = "drop"
        )
    # # 5. Average across k
    # improved_sp |>
    #     group_by(across(all_of(c(group_vars, panel_var)))) |>
    #     summarise(
    #         mean_n_above = mean(n_sp_above),
    #         sd_n_above   = sd(n_sp_above),
    #         mean_prop    = mean(n_sp_above / n_sp_total),
    #         n_k          = n(),
    #         .groups = "drop"
    #     )

}

# A function to compare strategies against the distribution of random designs.
# Mirrors the "metric_hist" figures of the optimality example: the histogram
# shows the metric over random-sampling runs, vertical lines show the other
# strategies, and (optionally) a dashed black line shows the reference model.
# ARGS:
#   - scores: a data.frame or tibble (e.g. compare_df).
#       Must contain columns: name, score, k, strat_var, size_var, panel_var.
#   - metric: a string. The name of a category in column "name".
#   - ref_scores: a data.frame or tibble (e.g. ref_df). Default NULL.
#       When given, its mean (per panel_var) is drawn as a dashed black line.
#   - random_strategy: a string. The value of strat_var forming the histogram.
#   - strat_var: a string. Column holding the strategies.
#   - strat_order: a vector of strings. Order of the strategies in the legend.
#   - size_var: a string. Column giving the number of new samples (rows of the grid).
#   - panel_var: a string. Column giving the subset (columns of the grid).
#   - panel_order: a vector of strings. Order of the panels.
#   - summary_fun: a function used to summarise a non-random strategy across
#       folds into one line (default: mean).
#   - lower_is_better: a boolean. Default NULL (FALSE for Ds, TRUE otherwise).
#   - bins: number of histogram bins.
#   - save_to: a string. Filepath to save the plot to. Must end with ".pdf".
distribution_plot <- function(
        scores, metric, ref_scores = NULL,
        random_strategy = "random-sampling",
        strat_var = "STRATEGIES", strat_order = NULL,
        size_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset", panel_order = NULL,
        summary_fun = mean,
        lower_is_better = FALSE,
        show_points = TRUE,
        ref_label = "Reference model",
        save_to = NULL) {

    # 1. Species-averaged score for each fold, strategy, size and panel
    keys <- c(panel_var, size_var, strat_var)
    fold_df <- scores |>
        filter(name == metric) |>
        group_by(across(all_of(c("k", keys)))) |>
        summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

    if (!is.null(panel_order)) {
        fold_df <- fold_df |> filter(.data[[panel_var]] %in% panel_order)
    }

    if (!any(fold_df[[strat_var]] == random_strategy)) {
        stop("No rows for random_strategy '", random_strategy, "' in `scores`.")
    }

    # 2. Reference line (one value per panel, repeated across sizes)
    ref_line <- NULL
    if (!is.null(ref_scores)) {
        ref_line <- ref_scores |>
            filter(name == metric) |>
            group_by(across(all_of(c("k", panel_var)))) |>
            summarise(score = mean(score, na.rm = TRUE), .groups = "drop") |>
            group_by(across(all_of(panel_var))) |>
            summarise(score = summary_fun(score, na.rm = TRUE), .groups = "drop")
        if (!is.null(panel_order)) {
            ref_line <- ref_line |> filter(.data[[panel_var]] %in% panel_order)
        }
    }

    # 3. Factor ordering: random sampling first, then the other strategies
    others <- setdiff(unique(fold_df[[strat_var]]), random_strategy)
    other_lvls <- if (is.null(strat_order)) {
        sort(others)
    } else {
        intersect(setdiff(strat_order, random_strategy), others)
    }
    lvls <- c(random_strategy, other_lvls)
    fold_df <- fold_df |> filter(.data[[strat_var]] %in% lvls)
    fold_df[[strat_var]] <- factor(fold_df[[strat_var]], levels = lvls)

    set_panels <- function(df) {
        if (is.null(df)) return(df)
        df[[panel_var]] <- if (is.null(panel_order)) {
            factor(df[[panel_var]])
        } else {
            factor(df[[panel_var]], levels = panel_order)
        }
        df
    }
    fold_df  <- set_panels(fold_df)
    ref_line <- set_panels(ref_line)

    # Fill colours: grey for random sampling, palette for the others
    fills <- c("grey70", PALETTE[seq_along(other_lvls)])
    names(fills) <- lvls

    # 4. Plot (built in the vertical orientation, then flipped)
    p <- ggplot(fold_df, aes(x = score, y = .data[[strat_var]],
                            fill = .data[[strat_var]])) +
        geom_violin(
            orientation = "y",
            scale = "width", trim = TRUE, colour = "grey30",
            alpha = 0.7, linewidth = 0.4)

    if (show_points) {
        p <- p + geom_jitter(
            width = 0, height = 0.12, size = 0.8,
            alpha = 0.5, colour = "grey20")
    }

    p <- p + stat_summary(
        fun = summary_fun, geom = "point", orientation = "y",
        shape = 23, size = 2.5, fill = "white", colour = "black")

    if (!is.null(ref_line)) {
        p <- p +
            geom_vline(
                data = ref_line,
                aes(xintercept = score, linetype = ref_label),
                colour = "black", linewidth = 0.7) +
            scale_linetype_manual(
                name = NULL,
                values = setNames("dashed", ref_label))
    }

    p <- p +
        facet_grid(
            rows = vars(.data[[size_var]]),
            cols = vars(.data[[panel_var]]),
            scales = "free_x",   # x = score, shared within each column
            labeller = labeller(.rows = label_both)) +
        labs(
            title = metric,
            subtitle = paste0(
                "(", if (lower_is_better) "lower" else "higher", " is better)"),
            x = paste0("Mean ", metric, " per fold"),
            y = NULL,
            caption = paste0(
                "Violins: distribution over folds ",
                "(species-averaged, one value per k_fold); ",
                "points: individual folds; diamonds: summary."))

    p <- my_custom_ggplot_theme(p) +
        scale_fill_manual(values = fills, guide = "none") +
        # first strategy at the top: reverse the discrete limits
        scale_y_discrete(limits = rev(lvls)) +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.99, 0.01),
            legend.justification = c(1, 0),
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.8, "lines"),
            legend.key.width = unit(1.5, "lines"),
            legend.background = element_rect(
                fill = scales::alpha("white", 0.85), colour = "grey70"),
            legend.margin = margin(2, 4, 2, 4))

    if (!is.null(save_to)) standardised_ggplot_save(p, save_to)
    list(data = list(folds = fold_df, reference = ref_line),
         plot = p)
}

# Spearman rank correlation. Returns c(rho, p-value), or c(NA, NA) when it 
# cannot be computed (fewer than 3 points, constant vector).
# ARGS:
#   - x, y: two vectors of the same size
safe_spearman <- function(x, y) {
    if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
        return(c(rho = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
}

# A function that plots pre-computed (Fisher approx) against realised metrics.
# Each point is one strategy, averaged over the k-folds. It reports the Spearman 
# rank correlation between them: does the approximation rank the strategies 
# like the retrained model does?
# ARGS:
#   - scores: a data.frame or tibble (e.g. compare_df).
#       Must contain columns: name, score, k, strat_var, size_var, panel_var.
#   - x_metric: a string. Name (in column "name") of the metric on the x-axis.
#   - y_metric: a string. Name (in column "name") of the metric on the y-axis.
#   - title, subtitle: strings. Plot title and subtitle.
#   - x_label, y_label: strings. Axis labels (default: the metric names).
#   - strat_var: a string. Column holding the strategies.
#   - strat_order: a vector of strings. Order of the strategies in the legend.
#   - size_var: a string. Column with the number of new samples.
#   - panel_var: a string. Column with the subset.
#   - panel_order: a vector of strings. Subsets to keep, and their order.
#   - random_strategy: a string. Strategy drawn in grey (small points).
#   - show_folds: a boolean. If TRUE, the value of each individual fold is 
#       drawn as a small, transparent point behind the mean over folds.
#   - save_to: a string. Filepath to save the plot to. Must end with ".pdf".
# RETURNS:
#   - a list: 
#       - data: one row per strategy (mean and SE across folds), 
#       - data_folds: one row per strategy and fold,
#       - stats: Spearman per panel, 
#       - plot.
#   One panel per (panel_var, size_var) combination, so that designs of
#   different sizes are never ranked against each other.
corr_metrics <- function(
        scores, x_metric, y_metric, title = NULL, subtitle = NULL,
        x_label = NULL, y_label = NULL,
        strat_var = "STRATEGIES", strat_order = NULL,
        size_var = "NEW_SAMPLE_SIZES",
        panel_var = "subset", panel_order = NULL,
        random_strategy = "random-sampling",
        show_folds = TRUE,
        save_to = NULL) {

    # 1. One value per design and metric (aggregated over species)
    keys <- unique(c("k", panel_var, size_var, strat_var))
    per_design <- function(metric, new_name) {
        scores |>
            filter(name == metric) |>
            group_by(across(all_of(keys))) |>
            summarise(value = mean(score, na.rm = TRUE), .groups = "drop") |>
            rename(!!new_name := value)
    }
    designs <- inner_join(
        per_design(x_metric, "x"),
        per_design(y_metric, "y"),
        by = keys) |>
        filter(is.finite(x), is.finite(y))

    # Keep only folds where all strategies are present in the panel, 
    # so that every strategy is averaged over the same folds
    designs <- designs |>
        group_by(across(all_of(c(panel_var, size_var)))) |>
        mutate(.n_strat = n_distinct(.data[[strat_var]])) |>
        group_by(across(all_of(c("k", panel_var, size_var)))) |>
        filter(n_distinct(.data[[strat_var]]) == first(.n_strat)) |>
        ungroup() |>
        select(-.n_strat)

    # Per-fold values, kept to draw the individual points
    designs_folds <- designs

    # Mean over folds, with the standard error across folds
    designs <- designs |>
        group_by(across(all_of(setdiff(keys, "k")))) |>
        summarise(
            n_folds = n(),
            x_se = sd(x) / sqrt(n()),
            y_se = sd(y) / sqrt(n()),
            x = mean(x),
            y = mean(y),
            .groups = "drop")

    # Stop if x_metric x y_metric does not exist
    if (nrow(designs) == 0) {
        stop("No design has both '", x_metric, "' and '", y_metric, "'.")
    }

    # Arange the designs per panel_order 
    if (!is.null(panel_order)) {
        designs <- designs |> filter(.data[[panel_var]] %in% panel_order)
        designs_folds <- designs_folds |> 
            filter(.data[[panel_var]] %in% panel_order)
    }

    # 2. Strategies: assign grey to random sampling, then get the others.
    others <- setdiff(unique(designs[[strat_var]]), random_strategy)
    other_lvls <- if (is.null(strat_order)) {
        sort(others)
    } else {
        c(intersect(strat_order, others), setdiff(others, strat_order))
    }
    lvls <- c(intersect(random_strategy, unique(designs[[strat_var]])), other_lvls)
    designs[[strat_var]] <- factor(designs[[strat_var]], levels = lvls)
    designs_folds[[strat_var]] <- factor(designs_folds[[strat_var]], levels = lvls)
    designs <- designs |> arrange(.data[[strat_var]])
    designs_folds <- designs_folds |> arrange(.data[[strat_var]])

    is_random <- lvls == random_strategy
    colors <- setNames(
        c(rep("grey40", sum(is_random)), rep_len(PALETTE, length(other_lvls))),
        lvls)
    point_sizes <- setNames(ifelse(is_random, 1.5, 3), lvls)

    # 3. Spearman correlation per panel (size_var x panel_var)
    # computed on the means over folds
    panel_stats <- designs_folds |>
        group_by(across(all_of(c(panel_var, size_var)))) |>
        summarise(
            n_designs = n(),
            rho = safe_spearman(x, y)[["rho"]],
            p = safe_spearman(x, y)[["p"]],
            .groups = "drop") |>
        mutate(.panel = sprintf(
            "%s: %s | %s: %s\nSpearman = %.2f (p = %.2f)",
            panel_var, .data[[panel_var]], size_var, .data[[size_var]], rho, p))
    if (!is.null(panel_order)) {
        panel_stats <- panel_stats |>
            mutate(across(all_of(panel_var), ~ factor(.x, levels = panel_order)))
    }
    panel_stats <- panel_stats |>
        arrange(across(all_of(c(panel_var, size_var)))) |>
        mutate(.panel = factor(.panel, levels = unique(.panel)))

    # Attach the panel labels to both tables
    panel_labels <- panel_stats |>
        select(all_of(c(panel_var, size_var)), .panel) |>
        mutate(across(all_of(panel_var), as.character))
    designs <- left_join(designs, panel_labels, by = c(panel_var, size_var))
    designs_folds <- left_join(
        designs_folds, panel_labels, by = c(panel_var, size_var))

    # 4. Plot
    # Individual folds: small, transparent, drawn first (= underneath), 
    # no legend entry. NULL (ignored by ggplot) when hidden.
    if (show_folds) {
        fold_layer <-  geom_point(
            data = designs_folds,
            aes(colour = .data[[strat_var]]),
            size = 0.8, alpha = 0.5, show.legend = FALSE)
    } else {
        fold_layer <- NULL
    }

    if (show_folds) {
        caption <- paste0(
            "Large points: mean over k-folds (bars: SE across folds).",
            "Small points: individual folds.",
            "Dotted line: linear fit on the means. ",
            "\nSpearman's rank correlation computed with each k_fold (not the means).")
    } else {
        caption <- paste0(
            "Large points: mean over k-folds (bars: SE across folds).",
            "Dotted line: linear fit on the means. ",
            "\nSpearman's rank correlation computed with each k_fold (not the means).")
    }

    p <- ggplot(designs, aes(x = x, y = y)) +
        geom_abline(
            intercept = 0, slope = 1,
            colour = "grey60", linetype = "dotted", linewidth = 0.6,
            inherit.aes = FALSE) +
        fold_layer +
        geom_linerange(
            aes(ymin = y - y_se, ymax = y + y_se, colour = .data[[strat_var]]),
            alpha = 0.4) +
        geom_linerange(
            aes(xmin = x - x_se, xmax = x + x_se, colour = .data[[strat_var]]),
            orientation = "y", alpha = 0.4) +
        geom_point(
            aes(colour = .data[[strat_var]], size = .data[[strat_var]]),
            alpha = 0.9) +
        scale_colour_manual(values = colors, name = "Strategy") +
        scale_size_manual(values = point_sizes, name = "Strategy") +
        facet_wrap(vars(.panel), scales = "free") +
        labs(
            title = title,
            subtitle = subtitle,
            x = if (is.null(x_label)) x_metric else x_label,
            y = if (is.null(y_label)) y_metric else y_label,
            caption = caption)

    p <- my_custom_ggplot_theme(p, with_palette = FALSE)

    if (!is.null(save_to)) standardised_ggplot_save(p, save_to)
    list(
        data = designs, data_folds = designs_folds, 
        stats = panel_stats, plot = p)
}

# A custom helper for transform_scale parameters (transforms a linear color
# scale into a mirrored log one, starting at offset_high).
# ARGS:
#   - offset_high: a numeric. The starting point of the log scale (reversed).
# RETURNS:
#   - a scales::new_transform() object.
mirrored_log_trans <- function(offset_high) {
  scales::new_transform(
    name      = paste0("mirrored_log_", offset_high),
    transform = function(x) -log(offset_high - x),
    inverse   = function(x) offset_high - exp(-x),
    domain    = c(-Inf, offset_high)  # strictly below offset_high
  )
}