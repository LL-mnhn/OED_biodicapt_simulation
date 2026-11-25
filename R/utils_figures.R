# Set of functions used to plot harmonised figures

##### Liraries ##### ----------------------------------------------------------
library(colorspace)
library(dotwhisker)
library(tidyterra)
library(patchwork)
library(ggspatial)
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
    ) {
    europe_shp = rnaturalearth::ne_countries(
        continent = "Europe", scale = "large", returnclass = "sf")
    france_shp = get_metropolitan_france_shapefile(borders_type)
    
    base_map <- ggplot(europe_shp) +
        geom_sf(fill = "grey80", color = "white") +                     # color of countries
        geom_sf(data = france_shp, fill = "white", color = "black") +   # color of France
        coord_sf(
            crs = PROJECT_CRS,              # CRS the map is drawn in
            xlim = c(lon_min, lon_max), 
            ylim = c(lat_min, lat_max)) +
        annotation_scale(
            location = "bl", width_hint = 0.25, 
            style = "ticks", text_family = FONT) +     
        annotation_north_arrow(
            location = "tr", which_north = "true",
            style = north_arrow_fancy_orienteering(text_family = FONT)) +
        theme_minimal() +
        theme(
            panel.background = element_rect(fill = "lightcyan1", color = NA),
            panel.grid.major = element_line(color = "lightcyan1"),
            panel.grid.minor = element_line(color = "lightcyan1")) +
        labs(
            x = "longitude",
            y = "latitude")
    
    return(my_custom_ggplot_theme(base_map, with_palette=FALSE, is_map = TRUE))
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

# A custom helper for transform_scale parameters (transforms a linear color
# scale into a mirrored scale one, starting at offset_high).
# ARGS:
#   - offset_high: a numeric. The starting point of the log scale (reversed).
# RETURNS:
#   - a scales::new_transform() object.
mirrored_sqrt_trans <- function(offset_high) {
  scales::new_transform(
    name      = paste0("mirrored_sqrt_", offset_high),
    transform = function(x) -sqrt(offset_high - x),
    inverse   = function(x) offset_high - ((-x)^2),
    domain    = c(-Inf, offset_high)  # strictly below offset_high
  )
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
        # uniques_vals <- unique(data[[column]])
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
        transform_scale = "identity",
        limits = NULL,
        cols = CUSTOM_CMAP,
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
        transform_scale = "identity",
        limits = NULL,
        cols = CUSTOM_CMAP,
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
        cols = CUSTOM_POINT_CMAP,
        transform_scale = "identity",
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
            "ggplot_quantitative_df_on_background_map() assumes that your data ",
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


##### HMSC training analysis ##### --------------------------------------------
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


##### HMSC result analysis ##### ----------------------------------------------
# A function to compute number of improving species's scores between k_folds, 
# subset and model type.
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
# RETURNS:
#   - A data.frame
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
    return(improved_sp)
}

# Computes spearman rank correlation. Returns c(rho, p-value), or c(NA, NA) 
# when it cannot be computed (fewer than 3 points, constant vector).
# ARGS:
#   - x, y: two vectors of the same size
safe_spearman <- function(x, y) {
    if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
        return(c(rho = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
}

### A series of helpers for 
###     corr_metrics(), distribution_plot() and compare_plot().
# A helper that creates a text based on lower_is_better and reversed booleans.
.direction_note <- function(lower_is_better, reversed = FALSE) {
    msg <- "("
    if (lower_is_better) {
        msg <- paste0(msg, "lower")
    } else {
        msg <- paste0(msg, "higher")
    } 
    msg <- paste0(msg, " is better")

    if (reversed) {
        msg <- paste0(msg, ", axis reversed")
    }
    msg <- paste0(msg, ")")
    return(msg)
}
# A helper to organize strategies: 
#   random strategy first, then strat_order, then anything that is left
.strategy_levels <- function(x, strat_order, random_strategy) {
    present <- unique(as.character(x))
    others <- setdiff(present, random_strategy)
    others <- if (is.null(strat_order)) {
        sort(others)
    } else {
        c(intersect(strat_order, others), setdiff(others, strat_order))
    }
    c(intersect(random_strategy, present), others)
}
# A helper to organize colours tied ti strategies: 
#   One colour per level, tied to the level (not to its position in a plot)
.strategy_colours <- function(lvls, random_strategy, random_colour = "grey60") {
    is_random <- lvls %in% random_strategy
    cols <- character(length(lvls))
    cols[is_random] <- random_colour
    cols[!is_random] <- rep_len(PALETTE, sum(!is_random))
    setNames(cols, lvls)
}
# A helper to select panels in a dataframe : 
#   Keeps only the requested panels (no-op if panel_var / panel_order is NULL)
.keep_panels <- function(df, panel_var, panel_order) {
    if (is.null(df) || is.null(panel_var) || is.null(panel_order)) return(df)
    df[df[[panel_var]] %in% panel_order, , drop = FALSE]
}
# A helper to turn a column into a factor (ordered if `order` is given)
.set_factor <- function(df, var, order = NULL) {
    if (is.null(df) || is.null(var) || !var %in% names(df)) return(df)
    df[[var]] <- if (is.null(order)) {
        factor(df[[var]])
    } else {
        factor(df[[var]], levels = order)
    }
    df
}
# A helper to have a common facet layout: 
#   Rows get "var: value" labels, columns only the value.
.facets <- function(rows = NULL, cols = NULL, scales = "free_x") {
    if (length(rows) + length(cols) == 0) return(NULL)
    facet_grid(
        rows = vars(!!!rlang::syms(rows)),
        cols = vars(!!!rlang::syms(cols)),
        scales = scales,
        labeller = labeller(.rows = label_both))
}
# A helper to replace "_" by spaces
.pretty_label <- function(x) {
    x |> 
        str_replace_all("_", " ") |> 
        str_to_sentence()
}
# A helper to explain violin plots
.violin_caption <- function(unit, show_folds) {
    msg <- paste0("Violins: distribution over folds (", unit, ").")

    if  (show_folds) {
        msg <- paste0(msg, " Points: individual folds.")
    }
    msg <- paste0(msg, ". Diamonds: summary.")
    return(msg)
}

# Plots a pre-computed (Fisher approx) metric against a realised one. Each
# diamond is one strategy, averaged over the k-folds. 
# Plots also report the Spearman rank correlation per panel.
# ARGS:
#   - scores: a data.frame or tibble with columns name, score, k, strat_var,
#       size_var, panel_var.
#   - x_metric, y_metric: metric names (in column "name") for the x / y axes.
#   - title: string (default title: "<y> vs <x>").
#   - x_label, y_label: axis labels (default: the metric names).
#   - strat_var / strat_order: strategy column and its legend order.
#   - size_var / size_order: column with the number of new samples (facet rows).
#   - panel_var / panel_order: subset column (facet columns), kept and ordered.
#   - random_strategy: strategy drawn in grey (small points).
#   - show_folds: draw each individual fold as a small transparent point.
#   - save_to: filepath ending in ".pdf".
#   - summary_fun: a function to summarise k_folds (mean or median usually)
# RETURNS: 
#   a list containing:
#       - data: summary_fun and SE across folds, 
#       - data_folds:one row per strategy and fold, 
#       - stats: Spearman per panel,
#       - plot: the resulting ggplot
corr_metrics <- function(
        scores, x_metric, y_metric,
        title = NULL, x_label = NULL, y_label = NULL,
        strat_var = "STRATEGIES", strat_order = NULL,
        size_var = "NEW_SAMPLE_SIZES", size_order = NULL,
        panel_var = "subset", panel_order = NULL,
        random_strategy = "random-sampling", summary_fun = mean,
        show_folds = TRUE, lower_is_better = FALSE,
        save_to = NULL) {

    scores <- .keep_panels(scores, panel_var, panel_order)

    # 1. One value per design and metric (aggregated over species)
    keys <- unique(c("k", panel_var, size_var, strat_var))
    per_design <- function(metric, new_name) {
        scores |>
            filter(name == metric) |>
            group_by(across(all_of(keys))) |>
            summarise(value = mean(score, na.rm = TRUE), .groups = "drop") |>
            rename(!!new_name := value)
    }
    folds_df <- inner_join(
        per_design(x_metric, "x"),
        per_design(y_metric, "y"),
        by = keys) |>
        filter(is.finite(x), is.finite(y))

    # Keep only folds where all strategies are present in the panel, so that
    # every strategy is averaged over the same folds
    folds_df <- folds_df |>
        group_by(across(all_of(c(panel_var, size_var)))) |>
        mutate(.n_strat = n_distinct(.data[[strat_var]])) |>
        group_by(across(all_of(c("k", panel_var, size_var)))) |>
        filter(n_distinct(.data[[strat_var]]) == first(.n_strat)) |>
        ungroup() |>
        select(-.n_strat)

    if (nrow(folds_df) == 0) {
        stop("No design has both '", x_metric, "' and '", y_metric, "'.")
    }

    # summary_fun over folds, with the standard error across folds
    summary_df <- folds_df |>
        group_by(across(all_of(setdiff(keys, "k")))) |>
        summarise(
            n_folds = n(),
            x_se = sd(x) / sqrt(n()),
            y_se = sd(y) / sqrt(n()),
            x = summary_fun(x),
            y = summary_fun(y),
            .groups = "drop")

    # 2. Strategies: grey for random sampling, palette for the others
    lvls <- .strategy_levels(summary_df[[strat_var]], strat_order, random_strategy)
    colors <- .strategy_colours(lvls, random_strategy)
    point_sizes <- setNames(ifelse(lvls %in% random_strategy, 1.5, 3), lvls)

    # 3. Spearman correlation per panel, computed on the individual folds
    stats <- folds_df |>
        group_by(across(all_of(c(panel_var, size_var)))) |>
        summarise(
            n_designs = n(),
            .sp = list(safe_spearman(x, y)),
            rho = .sp[[1]][["rho"]],
            p = .sp[[1]][["p"]],
            .groups = "drop") |>
        select(-.sp) |>
        mutate(label = sprintf(
            "Spearman = %.2f (%s)", rho,
            ifelse(is.na(p), "p = NA",
                   ifelse(p < 0.01, "p < 0.01", sprintf("p = %.2f", p)))))

    # Same factor levels everywhere, so facets line up across the three tables
    prep <- function(df) {
        df |>
            .set_factor(strat_var, lvls) |>
            .set_factor(size_var, size_order) |>
            .set_factor(panel_var, panel_order)
    }
    summary_df  <- prep(summary_df) |> arrange(.data[[strat_var]])
    folds_df <- prep(folds_df) |> arrange(.data[[strat_var]])
    stats    <- prep(stats)

    # 4. Plot
    # Should points for each individual k_fold should be shown?
    name_fn_summary <- deparse1(substitute(summary_fun))
    caption <- paste0(
        "Diamonds: ", name_fn_summary, 
        " over k-folds (bars: SE across folds)."
    )
    if (show_folds) {
        fold_layer <- geom_point(
            data = folds_df,
            aes(colour = .data[[strat_var]]),
            size = 0.8, alpha = 0.5, show.legend = FALSE)
        caption <- paste0(caption, " Small points: individual folds.")
    } else {
        fold_layer <- NULL
    }
    caption <- paste0(caption, 
        " Dotted line: y = x.",
        "\nSpearman's rank correlation computed across individual folds ",
        "(not the ", name_fn_summary, "s).")

    if (is.null(x_label)) {
        x_lab <- x_metric 
    } else {
        x_lab <- x_label
    } 
    if (is.null(y_label)) {
        y_lab <- y_metric 
    } else {
        y_lab <- y_label
    } 

    p <- ggplot(summary_df, aes(x = x, y = y)) +
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
            aes(fill = .data[[strat_var]], size = .data[[strat_var]]),
            shape = 23, colour = "black", alpha = 0.9) +
        # Spearman annotation in the top-left corner of each panel
        geom_text(
            data = stats, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3,
            inherit.aes = FALSE) +
        .facets(rows = panel_var, cols = size_var, scales = "free") +
        labs(
            title = if (is.null(title)) paste(y_lab, "vs", x_lab) else title,
            subtitle = .direction_note(lower_is_better),
            x = x_lab, y = y_lab,
            caption = caption)

    p <- my_custom_ggplot_theme(p, with_palette = FALSE) +
        scale_colour_manual(values = colors, guide = "none") +
        scale_fill_manual(values = colors, name = "Strategy") +
        scale_size_manual(values = point_sizes, name = "Strategy")

    if (!is.null(save_to)) {
        standardised_ggplot_save(p, save_to)
    }
    list(data = summary_df, data_folds = folds_df, stats = stats, plot = p)
}

# Distribution over k-folds of a (species-averaged) metric, one violin per
# strategy, with an optional reference-model line.
# ARGS:
#   - scores: data.frame with columns name, score, k, strat_var, size_var,
#       panel_var.
#   - metric: metric name (in column "name").
#   - ref_scores: optional data.frame (name, score, k, panel_var) of the
#       reference model, drawn as a dashed vertical line per panel.
#   - title: string (default: the metric).
#   - strat_var / strat_order, size_var / size_order, panel_var / panel_order,
#       random_strategy, show_folds, save_to: as in corr_metrics().
#   - summary_fun: summary drawn as a diamond (and used for the reference line).
#   - lower_is_better: TRUE/FALSE, or NULL to guess from the metric name.
#   - ref_label: legend label of the reference line.
# RETURNS:  
#   a list containing:
#       - data: list(folds, reference), 
#       - plot: the resulting ggplot
distribution_plot <- function(
        scores, metric, ref_scores = NULL, title = NULL,
        strat_var = "STRATEGIES", strat_order = NULL,
        size_var = "NEW_SAMPLE_SIZES", size_order = NULL,
        panel_var = "subset", panel_order = NULL,
        random_strategy = "random-sampling",
        show_folds = TRUE,
        summary_fun = mean,
        lower_is_better = FALSE,
        ref_label = "Reference model",
        save_to = NULL) {

    # 1. Species-averaged score for each fold, strategy, size and panel
    keys <- c(panel_var, size_var, strat_var)
    fold_df <- scores |>
        .keep_panels(panel_var, panel_order) |>
        filter(name == metric) |>
        group_by(across(all_of(c("k", keys)))) |>
        summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

    # 2. Reference line (one value per panel, repeated across sizes)
    ref_line <- NULL
    if (!is.null(ref_scores)) {
        ref_line <- ref_scores |>
            .keep_panels(panel_var, panel_order) |>
            filter(name == metric) |>
            group_by(across(all_of(c("k", panel_var)))) |>
            summarise(score = mean(score, na.rm = TRUE), .groups = "drop") |>
            group_by(across(all_of(panel_var))) |>
            summarise(score = summary_fun(score, na.rm = T), .groups = "drop")
    }

    # 3. Factor ordering and colours (shared with the other two plots)
    lvls <- .strategy_levels(
        fold_df[[strat_var]], strat_order, random_strategy)
    colors <- .strategy_colours(
        lvls, random_strategy, random_colour = "grey70")

    fold_df <- fold_df |>
        .set_factor(strat_var, lvls) |>
        .set_factor(size_var, size_order) |>
        .set_factor(panel_var, panel_order)
    ref_line <- .set_factor(ref_line, panel_var, panel_order)

    # 4. Plot (horizontal violins: score on x, strategies on y)
    p <- ggplot(fold_df, aes(y = score, x = .data[[strat_var]],
                             fill = .data[[strat_var]])) +
        geom_violin(
            orientation = "x",
            scale = "width", trim = TRUE, colour = "grey30",
            alpha = 0.7, linewidth = 0.4)

    if (show_folds) {
        p <- p + geom_jitter(
            width = 0.12, height = 0, size = 0.8,
            alpha = 0.5, colour = "grey20")
    }

    p <- p + stat_summary(
        fun = summary_fun, geom = "point", orientation = "x",
        shape = 23, size = 2.5, fill = "white", colour = "black")

    if (!is.null(ref_line)) {
        p <- p +
            geom_hline(
                data = ref_line,
                aes(yintercept = score, linetype = ref_label),
                colour = "black", linewidth = 0.7) +
            scale_linetype_manual(
                name = NULL,
                values = setNames("dashed", ref_label))
    }

    p <- p +
        .facets(rows = panel_var, cols = size_var, scales = "free_y") +
        labs(
            title = if (is.null(title)) metric else title,
            subtitle = .direction_note(lower_is_better),
            y = paste0("Mean ", metric, " per fold"),
            x = NULL,
            caption = .violin_caption(
                "species-averaged, one value per k_fold", show_folds))

    p <- my_custom_ggplot_theme(p, with_palette = FALSE) +
        scale_fill_manual(values = colors, guide = "none") +
        # first strategy at the top: reverse the discrete limits
        scale_x_discrete(limits = rev(lvls)) +
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
            legend.margin = margin(2, 4, 2, 4),
            axis.text.x = element_text(angle = 45, hjust = 1))

    if (!is.null(save_to)) standardised_ggplot_save(p, save_to)
    list(data = list(folds = fold_df, reference = ref_line), plot = p)
}

# Mean +/- SD over k-folds of the difference between compared scores and a
# reference model (line plot: difference on y, x_var on x, one line per
# strat_var level).
# ARGS:
#   - ref_scores, compare_scores: data.frames (name, score, k, ...) used to
#       compute the per-fold differences. Alternatively give `diffs_scores`,
#       a table of pre-computed differences with a `diff_score` column.
#   - metric: metric name (in column "name").
#   - title: string (default: the metric).
#   - x_var / x_order: column drawn on the x-axis (NULL = a single position).
#   - strat_var / strat_order: column mapped to the colour (one line per
#       level). NULL = a single line.
#   - size_var / size_order, panel_var / panel_order, random_strategy,
#       summary_fun, lower_is_better, save_to: as in distribution_plot().
#   - group_species: TRUE = average over species inside each fold; FALSE =
#       one column of facets per species.
# RETURNS:
#   a list containing
#       - data: per-fold differences,
#       - summary: mean (summary_fun) and SD per line point,
#       - plot: the resulting ggplot
compare_lineplot <- function(
        ref_scores = NULL, compare_scores = NULL, diffs_scores = NULL,
        metric,
        title = NULL, x_var = NULL, x_order = NULL,
        strat_var = NULL, strat_order = NULL,
        size_var = NULL, size_order = NULL,
        panel_var = NULL, panel_order = NULL,
        random_strategy = "random-sampling",
        group_species = TRUE,
        summary_fun = mean,
        lower_is_better = NULL,
        save_to = NULL) {

    if (!is.null(diffs_scores) && !is.null(compare_scores)) {
        stop("Cannot have diff_scores and compare_scores together")
    }
    if (!is.null(ref_scores) && is.null(compare_scores) && is.null(diffs_scores)) {
        stop("Cannot have ref_scores without compare_scores or diffs_scores")
    }
    
    if (group_species) {
        sp_key <- NULL
    } else {
        sp_key <- "sp"
    }
    join_keys <- c("k", panel_var, sp_key)

    ref_summary <- NULL
    if (!is.null(ref_scores)) {
        # Reference: mean score per fold (and panel, and species if kept)
        ref_summary <- ref_scores |>
            .keep_panels(panel_var, panel_order) |>
            filter(name == metric) |>
            group_by(across(all_of(join_keys))) |>
            summarise(ref_score = mean(score, na.rm = TRUE), .groups = "drop")
    }
    if (!is.null(compare_scores)) {
        # Compared scores: mean per fold, x, strat, size, panel (and species)
        compare_summary <- compare_scores |>
            .keep_panels(panel_var, panel_order) |>
            filter(name == metric) |>
            group_by(across(all_of(unique(c(join_keys, size_var, x_var, strat_var))))) |>
            summarise(avg_score = mean(score, na.rm = TRUE), .groups = "drop")

        # Per-fold difference with the matching reference
        diffs_df <- compare_summary |>
            left_join(ref_summary, by = join_keys) |>
            mutate(diff_score = avg_score - ref_score)
    } 
    if (!is.null(diffs_scores)) {
        diffs_df <- diffs_scores |>
            filter(name == metric) |>
            .keep_panels(panel_var, panel_order)
    }

    # dummy x / strat so the lines still have a discrete position and a colour
    x_is_numeric <- FALSE
    if (is.null(x_var)) {
        diffs_df$.x <- factor("all")
        x_col <- ".x"
    } else {
        x_col <- x_var
        # numeric, and no explicit order requested -> keep as continuous
        x_is_numeric <- is.numeric(diffs_df[[x_var]]) && is.null(x_order)
        if (!x_is_numeric) {
            diffs_df <- .set_factor(diffs_df, x_var, x_order)
        }
    }
    if (is.null(strat_var)) {
        diffs_df$.strat <- factor("all")
        strat_col <- ".strat"
    } else {
        lvls_strat <- .strategy_levels(
            diffs_df[[strat_var]], strat_order, random_strategy)
        diffs_df <- .set_factor(diffs_df, strat_var, lvls_strat)
        strat_col <- strat_var
    }
    diffs_df <- diffs_df |>
        .set_factor(size_var, size_order) |>
        .set_factor(panel_var, panel_order)
    
    # One reference value per facet cell (panel, and species if not grouped)
    ref_lines <- NULL
    if (!is.null(ref_summary)) {
        label_cols <- c(panel_var, sp_key)
        ref_lines <- ref_summary |>
            group_by(across(all_of(label_cols))) |>
            summarise(ref_mean = mean(ref_score, na.rm = TRUE),
                      .groups = "drop") |>
            .set_factor(panel_var, panel_order) |>
            arrange(across(all_of(label_cols)))
    }

    # 4. Caption, with reference means per panel
    bottom_caption <- paste0(
        if (group_species) {
            "species-averaged within each fold"
        } else {
            "computed per species"
        },
        "; points = mean, bars/ribbon = \u00b1 1 SD across k_folds")

    y_label <- if (is.null(diffs_scores)) "Delta in average " else ""

    if (!is.null(ref_lines)) {
        label_cols <- c(panel_var, sp_key)
        ref_labels <- if (length(label_cols) > 0) {
            do.call(paste, c(as.list(ref_lines[label_cols]), sep = " / "))
        } else {
            NULL
        }
        ref_text <- if (is.null(ref_labels)) {
            format(round(ref_lines$ref_mean, 3), nsmall = 3)
        } else {
            paste0(ref_labels, ": ", round(ref_lines$ref_mean, 3))
        }
        bottom_caption <- paste0(
            bottom_caption,
            "\nReference ", metric, " (dashed line",
            if (length(label_cols) > 0) paste0(", by ", paste(label_cols, collapse = " / ")),
            "): ", paste(ref_text, collapse = ", "), ".")
    }

    # one mean / SD per (x, strat, panel, size, species) across the k folds
    summary_keys <- unique(c(x_col, strat_col, panel_var, size_var, sp_key))
    line_summary <- diffs_df |>
        group_by(across(all_of(summary_keys))) |>
        summarise(
            m = summary_fun(diff_score, na.rm = TRUE),
            s = sd(diff_score, na.rm = TRUE),
            .groups = "drop") |>
        mutate(ymin = m - s, ymax = m + s)

    # 5. Colours: tied to values, grey for the random strategy
    lvls_now <- levels(diffs_df[[strat_col]])
    line_colors <- .strategy_colours(lvls_now, random_strategy, random_colour = "grey50")

    # 6. Plot
    dodge_width <- if (x_is_numeric) {
        xs <- sort(unique(line_summary[[x_col]]))
        if (length(xs) > 1) min(diff(xs)) * 0.3 else 0.3
    } else {
        0.3
    }
    dodge <- position_dodge(width = dodge_width)
    eb_width <- if (x_is_numeric) dodge_width * 0.6 else 0.2

    p <- ggplot(
            line_summary,
            aes(x = .data[[x_col]], y = m,
                colour = .data[[strat_col]], fill = .data[[strat_col]],
                group = .data[[strat_col]]))

    p <- p +
        geom_line(position = dodge, linewidth = 0.8) +
        geom_errorbar(
            aes(ymin = ymin, ymax = ymax),
            position = dodge, width = eb_width, linewidth = 0.5) +
        geom_point(position = dodge, size = 2)
    
    if (!is.null(ref_lines)) {
        if (!is.null(diffs_scores)) {
            p <- p + geom_hline(
                data = ref_lines, aes(yintercept = ref_mean),
                linetype = "dashed", colour = "black", inherit.aes = FALSE)
        } else {
            p <- p + geom_hline(yintercept = 0, linetype = "dashed")
        }
    }
        
    p <- p + labs(
        title = if (is.null(title)) metric else title,
        subtitle = .direction_note(lower_is_better),
        caption = bottom_caption,
        colour = if (is.null(strat_var)) NULL else .pretty_label(strat_var),
        fill = if (is.null(strat_var)) NULL else .pretty_label(strat_var),
        x = if (is.null(x_var)) NULL else .pretty_label(x_var),
        y = paste0(y_label, metric))

    # faceting: rows = size (and panel if species take the columns)
    p <- p + if (group_species) {
        .facets(cols = panel_var, scales = "free_y")
    } else {
        .facets(cols = panel_var, rows = "sp", scales = "free_y")
    }

    guide_type <- if (is.null(strat_var)) "none" else "legend"
    p <- my_custom_ggplot_theme(p, with_palette = FALSE) +
        scale_colour_manual(values = line_colors, guide = guide_type) +
        scale_fill_manual(values = line_colors, guide = "none")

    if (is.null(x_var)) {
        p <- p + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())
    } else if (x_is_numeric) {
        p <- p + scale_x_continuous(
            breaks = sort(unique(line_summary[[x_col]])))
    } else {
        p <- p + scale_x_discrete(labels = abbreviate_labels)
    }

    if (!is.null(save_to)) standardised_ggplot_save(p, save_to)
    list(
        data = diffs_df, 
        ref_data = ref_lines,
        summary = line_summary, 
        plot = p)
}
