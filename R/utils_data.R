# Set of utility functions used accross scripts.

##### Liraries #####
library(sf)
library(cli)
library(readr)
library(dplyr)
library(gstat)
library(terra)
library(purrr)
library(tools)
library(readxl)
library(stringr)
library(dggridR)
library(tidyterra)
library(rnaturalearth)
library(exactextractr)


##### Parameters #####
source(here::here("data/config/config.R")) # all parameters are grouped together


##### Global functions #####
# A function to ask for a inputs by a user, which works in both 
# interactive (with user inputs) and batch mode (on remote cluster)
# ARGS:
#   - prompt: the message to display before asking for input
#   - remote_username : the name of the user for which there should not be an 
#       interactive session (check with Sys.info()["user"] in remote terminal).
# RETURNS:
#   - a string with the user answer or auto "Yes" when on a remote machine.
typeline <- function(prompt, remote_username = "", default_remote_behaviour = "No") {
    if (Sys.info()["user"] == remote_username) {
        # auto-accept on cluster runs
        cli_alert_warning(
            paste0("Cluster detected: automatic answer: '", 
                   default_remote_behaviour, "'."))
        txt <- default_remote_behaviour
    } else if (interactive() ) {
        txt <- readline(prompt)
    } else {
        cat(prompt)
        txt <- readLines("stdin", n=1)
    }
    return(txt)
}

# A function to ask the user if he wants to overwrite a file.
# ARGS:
#   - path: a path (folder or file) to check.
#   - remote_username : the name of the user for which there should not be an 
#       interactive session (check with Sys.info()["user"] in remote terminal).
# RETURNS:
#   - a boolean (TRUE to authorize, FALSE to refuse).
authorise_overwrite <- function(
    path, remote_username = REMOTE_USERNAME, default_remote_behaviour = "No") {
    # if file does not exist, no overwrite needed, return TRUE
    if (!file.exists(path)) {
        return(TRUE)
    } else if (dir.exists(path) && (length(list.files(path, all.files = TRUE, no.. = TRUE)) == 0)) {
        # is an empty folder (then no worries, we overwrite nothing)
        return(TRUE)
    }else { 
        # if file exist, ask user for what needs to be done
        user_input <- typeline(
            prompt = paste0(" Overwrite `", path, "`? [Y/n]: "),
            remote_username = remote_username,
            default_remote_behaviour = default_remote_behaviour)
        cleaned_input <- tolower(trimws(user_input))

        if (cleaned_input %in% c("y", "yes")) {
            cli_alert_info("User authorized process to overwrite file(s).")
            return(TRUE)

        } else if (cleaned_input %in% c("n", "no")) {
            cli_alert_info("User refused to allow process to overwrite file(s).")
            return(FALSE)

        } else {
            # if answer cannot be identified, return FALSE
            cli_alert_warning(paste0("Answer was '", user_input,  
                "', expected Y(es)/N(o)."))
            cli_alert_warning("\t↳ Defaults to `FALSE` (no overwrite).")
            return(FALSE)
        }
    }
}

# A function that fixes name formatting for dataframes columns, character 
# vector or SpatRasters/SpatVector layers
# ARGS:
#   - x: data.frame, SpatRaster, SpatVector or character
# RETURNS:
#   - The same object, with modified names
fix_names <- function(x) {
    if (    inherits(x, "SpatRaster") || 
            inherits(x, "SpatVector") || 
            inherits(x, "data.frame")) {
        names(x) <- tolower(make.names(names(x)))
    } else if( inherits(x, "character")) {
        x <- tolower(make.names(x))
    } else {
        stop("Unsupported class: ", paste(class(x), collapse = "/"),
            ". Expected character, SpatRaster, SpatVector, or data.frame.")
    }
    return(x)
}

# A function that splits a number into as-equal-as-possible integer parts
# ARGS:
#   - number: the number (integer) to divide into integers
#   - divisions: the number (integer) of divisions to make
# RETURNS:
#   - a vector of number, each item is the size of a division
split_evenly <- function(number, divisions) {
    base <- number %/% divisions        # integer division
    remainder <- number %% divisions    # remainder

    parts <- c(rep(base + 1, remainder), rep(base, divisions - remainder))
    return(parts)
}


##### Datasets #####
# A function that "blurs" coordinate within a dataframe by randomly shifting
#   each longitude and latitude coordinates.
# ARGS:
#   - df: a dataframe with gps coordinates.
#   - x_lon: a string. The name of the column with longitude values.
#   - y_lat: a string. The name of the colum with latitude values.
#   - blur_km: an integer/a float. 
#       The range within which the blurring occurs (in km).
#   - seed: an integer. When given, sets a seed at function level.
#   - lon_min, lon_max, lat_min, lat_max: coordinate extent of df. Used to 
#       define mean latitude and approximate res_km.
# RETURNS:
#   - orginal df with shifted coordinates.
blur_coordinates <- function(
        df, 
        x_lon, y_lat, 
        blur_km, seed = NULL, 
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
    # for reproducible results  
    if (!is.null(seed)){
        set.seed(seed) 
    }

    # if sf object, remove geometry field (we're gonna make a new one)
    is_sf <- inherits(df, "sf")
    if (is_sf) {
        crs_orig <- sf::st_crs(df)
        geom_col <- attr(df, "sf_column")
        df <- sf::st_drop_geometry(df) 
    }
  
    # Manual check of coordinate system (CRS 4326)
    if (    !between(min(df[[x_lon]]), lon_min, lon_max) ||
            !between(max(df[[x_lon]]), lon_min, lon_max) ||
            !between(min(df[[y_lat]]), lat_min, lat_max) ||
            !between(max(df[[y_lat]]), lat_min, lat_max)) {
        cli_alert_danger("Dataset failed basic extent test")
        stop(paste0(
            "Wrong coordinate system and/or borders are outside", 
            " France's metropolitan area."))
    }

    # get mean latitude and longitude resolutions
    lat_mean  <- (lat_min + lat_max) / 2
    res_lat   <- blur_km / 111.0
    res_lon   <- blur_km / (111.0 * cos(lat_mean * pi / 180))

    # Dividing by 1.96 so ~95% of points stay within a circle of blur_km 
    df <- df |>
        mutate(
            "{x_lon}" := .data[[x_lon]] + rnorm(dplyr::n(), mean = 0, sd = res_lon / 1.96),
            "{y_lat}" := .data[[y_lat]] + rnorm(dplyr::n(), mean = 0, sd = res_lat / 1.96)
        )
    
    # make a new "geometry" column
    if (is_sf) {
        df <- st_as_sf(
            df, coords = c(x_lon, y_lat), crs = crs_orig, remove = FALSE)
        names(df)[names(df) == "geometry"] <- geom_col
        st_geometry(df) <- geom_col
    }
    return(df)
}

# A function to aggregate rasters together (months -> year)
# ARGS:
#   - raster_paths: filepaths to 12 (monthly) rasters of the same data.
#   - buffer: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative) 
#       to avoid discarding useful data.
#   - verbose: a boolean. If TRUE, shows info messages.
#   - fun: a function. Tells how the data should be aggregated together.
# RETURNS:
#   - A single spatRaster with the aggregated data.
monthly_2_yearly_rasters <- function(
        raster_paths, 
        buffer, 
        verbose = TRUE, 
        fun = mean){
    
    # Verify we have 12 paths (for 12 months)
    if (length(raster_paths) != 12){
        stop(paste("List of strings recieved contains", length(raster_paths), "elements, expected 12."))
    }

    # Stack rasters together
    raster_list <- c()
    for (n_path in 1:12){
        # import raster
        raster <- rast(raster_paths[n_path])
        
        if (!same.crs(raster, "EPSG:4326")){
            stop(paste0("Raster (at ", raster_paths[n_path], ") has crs ", 
                crs(raster), ". Expected 'EPSG:4326'."))
        }
        
        cropped_raster <- clip_raster_france_wgs84_crs(
            raster, verbose = verbose, buffer = buffer)     

        if (n_path == 1){
            rasters_together <- cropped_raster
        } else {
            raster_list <- c(rasters_together, cropped_raster)
        }
    }
    
    # Compute and return average/median/other on list of cropped raster
    annual_raster <- app(rasters_together, fun = fun)

    return(annual_raster)
}

# A function to standardize monthly raster aggregation and save path.
# ARGS:
#   - folder_path: path to a folder containing 12 raster files.
#   - save_to: a filepath to save the new raster (must end with .tif).
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
# RETURNS:
#   - A single spatRaster with the aggregated data.
average_monthly_quantitative_rasters <- function(
        folder_path,
        save_to,
        res_km = RES_KM) {
    cli_alert_info("Aggregation of monthly data.")
    
    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0(
            "Skipping pre-processing of monthly '", 
            folder_path, "' data.\n\n"))

    } else {
        files <- list.files(folder_path)
        twelve_paths <- file.path(folder_path, files)
            
        # default format is per month, we need yearly data
        cli_alert_info("Averaging monthly rasters over a year...")
        annual_raster <- monthly_2_yearly_rasters(
            raster_paths = twelve_paths, 
            buffer = res_km*5, 
            verbose = FALSE,
            fun = mean)
        
        cli_alert_info("Saving file...")
        writeRaster(annual_raster,  save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created '", save_to, "' file!\n\n"))
    }
}

# A function that extracts features from all preprocessed rasters (hard coded 
# paths) and matches them to the coordinates in the given (preprocessed as well)
# dataframe.
# ARGS:
#   - file: a string or a data.frame.
#       Must contain columns "LON" and "LAT" (crs = 4326)
#       If a string, reads it as a .csv file. If a data.frame, uses it directly. 
#   - save_to: the path where the output sf will be saved (use ".gpkg")
#   - mode: a string. Three modes available:
#       - "raw": extracts features from each raster best available resolution
#       - "raster": extracts features from each raster standardised resolution
#       - "hexagonal": extracts features from each hexagonal projections
#   - blur: a boolean. If TRUE, blurs the coordinates of the original dataframe
#       before saving them/returning them (but after the extraction of the
#       features to be as precise as possible).
# RETURNS:
#   - a sf object, the dataframe with its added features (= columns)
save_features_from_obs <- function(file, save_to, mode = "raw", blur = TRUE) {
    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0("Skipping the addition of features to ", 
        basename(save_to), ".\n\n"))

    } else {
        # load file (either a datframe or a filepath)
        if (!is.data.frame(file))  {
            df <- read_csv(file, show_col_types = FALSE)
        } else {
            df <- file
        }

        # convert dataframe to geo dataframe
        cli_alert_info(paste0("Conversion to shapefile..."))
        points <- vect(st_as_sf(
            df, coords = c("LON", "LAT"), crs = 4326, remove = FALSE))

        # the user can choose which type of file the features must be extracted
        # from. If raw, uses the initial raw resolution of the files (which are
        # usually finer than 1km² in resolution). If using raster or hexagonal,
        # uses pre-processed data that was made to be pixels/hexagons space out
        # by RES_KM.
        if (mode == "raw") {
            environmental_data_paths <- c(
                CORINE_PATH_PREPROCESSED, CHELSA_ANNUAL_PATHS, NDVI_ANNUAL_PATH,
                ELEVATION_PATH_PREPROCESSED, LIGHT_POLLUTION_FILEPATH)
            
        } else if (mode == "raster") {
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_res10km-WGS84.tif")
            
        } else if (mode == "hexagonal") {
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_hexagons_res10km-WGS84.gpkg")
            
        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }

        # new names are needed when adding columns to 'points' (geo dataframe) 
        environmental_data_names <- c(
            "clc", paste0("chelsa_", CHELSA_DATASETS), "ndvi", "elevation",
            "light_pollution")   
        
        cli_bullets(c(
            "!" = "[Fetching features]: Using hard coded paths:",
            setNames(environmental_data_paths, rep("*", length(environmental_data_paths)))
        ))

        for (i in 1:length(environmental_data_paths)) {
            
            if (mode %in% c("raw", "raster")) {
                # load envrionmental data
                env_data <- rast(environmental_data_paths[i])

                # drop the leading ID column, keep one column per layer
                extracted_values <- terra::extract(
                    env_data, points)[, -1, drop = FALSE]
                
                # Extract layer names
                n_layers <- ncol(extracted_values)
                if (n_layers == 1) {
                    col_names <- environmental_data_names[i]
                } else {
                    # if many layers, take the name of the layers in the raster
                    layer_ids <- names(env_data)
                    col_names <- layer_ids
                }
                colnames(extracted_values) <- col_names

                for (col in col_names) {
                    points[[col]] <- extracted_values[[col]]
                }

            } else if (mode == "hexagonal") {
                env_data <- vect(environmental_data_paths[i])

                # extract values of points projected onto env_data cells
                extracted_values <- terra::extract(
                    env_data, points)[, "dominant_class"]

                if (length(extracted_values) != nrow(points)) {
                    stop("Overlaping shapes in .gpkg file.")
                }
            }
            points[[environmental_data_names[i]]] <- extracted_values
        }

        if (blur) {
            # replace raw coordinates with blurred coordinates
            cli_alert_info("Data anonymization...")
            source(here::here(file.path("data", "config", "seed.R"))) 
            anonym_df <- blur_coordinates(
                points, "LON", "LAT", BLUR_SIZE, BLUR_SEED)
        } else {
            anonym_df <- points
        }

        # Save resulting dataframe
        cli_alert_info("Saving file...")
        writeVector(anonym_df, save_to, overwrite = TRUE)
        cli_alert_success("Saved anonymized dataset with features!")

        return(anonym_df)
    }
}

# A function that stacks together all preprocessed rasters/sf (hard coded paths)
# into one single raster/sf.
# ARGS:
#   - save_to: the path where the output file will be saved 
#       (use ".gpkg" if mode is "hexagonal", ".tif" for "raw" and "raster").
#   - mode: a string. Three modes available:
#       - "raw": builds a raster with the best possible resolution (uses the 
#           highest resolution of all preprocessed rasters).
#       - "raster": builds a raster with from standardised resolutions
#       - "hexagonal": extracts features from each standardised hexagonal sf
#   - lon_min, lon_max, lat_min, lat_max: coordinate extent of df. Used to 
#       define mean latitude and approximate res_km.
# RETURNS:
#   - a sf object, the dataframe with its added features (= columns)
create_master_files_features <- function(
        save_to, mode = "raw",
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX
    ) {
    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0(
            "Skipping the creation of a raster containing containing all the ",
            "features together (same CRS, same pixels).\n\n"))

    } else {
        # the user can choose which type of file the features must be extracted
        # from. If raw, uses the initial raw resolution of the files (which are
        # usually finer than 1km² in resolution). If using raster or hexagonal,
        # uses pre-processed data that was made to be pixels/hexagons space out
        # by RES_KM.

        # Get paths of datasets (outside the main "merging" if/else to be able 
        # to plot the paths used before beginning the computation).
        if (mode == "raw" | mode == "raster") {
            if (mode == "raw") {
                environmental_data_paths <- c(
                    CORINE_PATH_PREPROCESSED, CHELSA_ANNUAL_PATHS, NDVI_ANNUAL_PATH,
                    ELEVATION_PATH_PREPROCESSED, LIGHT_POLLUTION_FILEPATH)
                
            } else { # mode == "raster"
                environmental_data_paths <- c(
                    CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                    NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                    LIGHT_POLLUTION_PROJECTION_BASEPATH)  
                environmental_data_paths <- paste0(
                    environmental_data_paths, 
                    "_projection_france_res", RES_KM, "km-WGS84.tif")
            }
        } else if (mode == "hexagonal") {
            # Shapefiles must be handled differently than rasters
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg")
        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }

        cli_bullets(c(
            "!" = "[Fetching features]: Using hard coded paths:",
            setNames(environmental_data_paths, rep("*", length(environmental_data_paths)))
        ))

        # Merge them together
        if (mode == "raw" | mode == "raster") {
            # Rasters can be handled the same way (not depending on the mode)
            # define target grid 
            target_ext <- ext(lon_min, lon_max, lat_min, lat_max)
            resolutions <- sapply(environmental_data_paths, function(p) {
                r <- rast(p)
                res(r)[1]
            })
            
            # get coarsest resolution (aka the most pixelated)
            target_res_deg <- max(resolutions)  
            template <- rast(
                target_ext, 
                resolution = target_res_deg,
                crs = crs(rast(environmental_data_paths[1]))) 
            cli_alert_info(paste0(
                "Selecting raster with coarsest resolution (res: ", 
                round(target_res_deg, 2), "°)."))
            
            # Crop and resample each raster onto the template
            cli_alert_info("Resampling rasters, can take several minutes...")
            raster_list <- list()
            for (path in environmental_data_paths) {
                env_raster <- rast(path)
                raster_cropped <- crop(env_raster, target_ext)

                # categorical vs continuous resampling method
                if (grepl("CLC", path, ignore.case = TRUE)) {
                    method <- "near"
                } else {
                    method <- "bilinear"
                }

                raster_resampled <- resample(
                    raster_cropped, template, method = method)
                if (grepl("CHELSA", path, ignore.case = TRUE)) {
                    name <- paste(
                        str_split(basename(path), "_")[[1]][c(1,3)], 
                        collapse = "_")
                } else if (grepl("CLC", path, ignore.case = TRUE)) {
                    name <- "CLC"
                } else if (grepl("light", path, ignore.case = TRUE)) {
                    name <- "light_pollution"
                } else if (grepl("Elevation", path, ignore.case = TRUE)) {
                    name <- "Elevation"
                } else if (grepl("ndvi", path, ignore.case = TRUE)) {
                    name <- "ndvi"
                } else {
                    stop(paste0("I'm a dummy, I did not expect new names, ",
                    "I gotta modify this function..."))
                }

                if (length(names(raster_resampled)) > 1) {
                    names(raster_resampled) <- c(
                        name, 
                        names(raster_resampled)[2:length(names(raster_resampled))])
                } else {
                    names(raster_resampled) <- name 
                }
                
                # also rename the active category column so it survives the round trip
                if (length(names(raster_resampled)) > 1) {
                    if (is.factor(raster_resampled[[1]])) {
                        cat_table <- cats(raster_resampled[[1]])[[1]]
                        active_col <- activeCat(raster_resampled[[1]]) + 1  # +1: activeCat is 0-indexed, col 1 is "value"
                        colnames(cat_table)[active_col] <- name
                        levels(raster_resampled[[1]]) <- cat_table
                    }
                } else {
                    if (is.factor(raster_resampled)) {
                        cat_table <- cats(raster_resampled)[[1]]
                        active_col <- activeCat(raster_resampled) + 1  # +1: activeCat is 0-indexed, col 1 is "value"
                        colnames(cat_table)[active_col] <- name
                        levels(raster_resampled) <- cat_table
                    }
                }

                raster_list[[name]] <- raster_resampled
            }

            # Stack into one raster
            names(raster_list) <- NULL # reset names to avoid auto-renaming
            master_raster <- rast(raster_list)
            
            # Save
            cli_alert_success("Saving raster...")
            writeRaster(master_raster, save_to, overwrite = TRUE)
            cli_alert_success("Saved stacked raster!\n\n")
            return(master_raster)

        } else if (mode == "hexagonal") {
            cli_alert_info("Reading hexagon layers and merging attributes...")
            hex_list <- lapply(environmental_data_paths, vect) 

            # keep "seqnum" for hexagon ID and keep_col for one useful column
            join_col <- "seqnum"  
            keep_col <- "dominant_class"

            # Make a list of polygon objects
            hex_list <- Map(
                function(h, path) {
                    df <- as.data.frame(h)
                    df <- df %>% select(contains(join_col), contains(keep_col))

                    if (grepl("CHELSA", path, ignore.case = TRUE)) {
                        prefix <- paste(
                            str_split(basename(path), "_")[[1]][c(1,3)], 
                            collapse = "_")
                    } else if (grepl("CLC", path, ignore.case = TRUE)) {
                        prefix <- "CLC"
                    } else if (grepl("light", path, ignore.case = TRUE)) {
                        prefix <- "light_pollution"
                    } else if (grepl("Elevation", path, ignore.case = TRUE)) {
                        prefix <- "Elevation"
                    } else if (grepl("ndvi", path, ignore.case = TRUE)) {
                        prefix <- "ndvi"
                    } else {
                        stop(paste0("I'm a dummy, I did not excpect new names, ",
                        "you gotta modify the function..."))
                    }

                    # rename columns
                    names(df)[names(df) == keep_col] <- paste0(prefix, "_", keep_col)

                    # put renamed attributes back onto the SpatVector's geometry
                    h_renamed <- h
                    values(h_renamed) <- df
                    h_renamed
                }, 
                hex_list, 
                environmental_data_paths)
            
            # sanity check: same hexagon grid across all layers
            n_geom <- sapply(hex_list, nrow)
            if (length(unique(n_geom)) > 1) {
                stop("Hexagon layers do not share the same grid (different feature counts).")
            }

            master_hex <- hex_list[[1]]
            for (h in hex_list[-1]) {
                # drop geometry duplication, join by hexagon id
                df <- as.data.frame(h)
                master_hex <- merge(master_hex, df, by = "seqnum")
            }

            cli_alert_success("Saving hexagon layer...")
            writeVector(master_hex, save_to, overwrite = TRUE)
            cli_alert_success("Saved stacked hexagons!\n\n")
            return(master_hex)

        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }
    }
}

# A function that automatically loads and prepares a "master" raster (i.e. a
# multi-layered raster loaded with environmental features).
# ARGS:
#   - mode: a string, either 'raw' or 'res'. Loads either the master_raster
#       with the best resolution available (raw) or the downsampled one (res).
#   - rm.lyr: a vector of string, the names of the layers that must be removed.
#   - rm.cor.lyr: a list of vector of >=2 strings. It contains names of columns
#       that are correlated together. The function keep the first name of each
#       vector, and removes the others.
# RETURNS:
#   - a multi-layer SpatRaster.
get_master_raster <- function(mode = "raw", rm.lyr = NULL, rm.cor.lyr = NULL) {
    # Load master raster
    if (mode == "raw") {
        cli_alert_info("Loading multi-layer raster with best resolution.")
        master_raster <- rast(FEATURES_RAW_FULL)
    } else if (mode == "res") {
        cli_alert_info(paste0(
            "Loading multi-layer raster (res = ", RES_KM, "km)."))
        master_raster <- rast(FEATURES_RES_FULL)
    } else {
        stop("The 'mode' parameter should be one of: 'raw', 'res'.")
    }
    
    # reduce size for faster computation
    master_raster <- clip_raster_france_wgs84_crs(  
            master_raster,  
            buffer = 5*RES_KM,   
            verbose = FALSE,
            save_to = NULL)
    
    # Remove useless layers
    for (layer in rm.lyr) {
        master_raster <- subset(
            master_raster, 
            names(master_raster)[names(master_raster) != layer])
    }

    # Remove correlated layers 
    for (correlated_layers in rm.cor.lyr) {
        for (layer in correlated_layers[2:length(correlated_layers)]) {
            master_raster <- subset(
                master_raster, 
                names(master_raster)[names(master_raster) != layer])
        }
    }
 
    return(fix_names(master_raster))
}


##### Rasters/SF functions #####
# A function that creates an empty raster that matches France's extent (WGS 84)
# ARGS:
#   - res_km: an integer/a float. 
#       The size for each pixel of the raster (in km, approximative).
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in EPSG:4326).
# RETURNS:
#   - a SpatRaster with empty values to the extent of France's metropolitan
#       area.
get_france_raster_template <- function(
        res_km,     
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
    # get mean latitude and longitude resolutions
    lat_mean <- (lat_min + lat_max) / 2
    res_lat <- res_km / 111.0
    res_lon <- res_km / (111.0 * cos(lat_mean * pi / 180))
    
    # create empty raster with template grid
    template <- rast(
        extent = ext(lon_min, lon_max, lat_min, lat_max),
        resolution = c(res_lon, res_lat),
        crs = "EPSG:4326"
    )
    return(template)
}

# A function that imports a shapefile of france's borders
# ARGS:
#   - borders: a string. Either "national" (default) or "regional". 
#       If "regional", draws highest level inner borders ofthe country.
# RETURNS:
#   - a sf object of France's metropolitan borders.
# WARNING:
#   - Hard coded extent values exclude all* overseas territories (*but corsica).
get_metropolitan_france_shapefile <- function(borders = "national") {
    # import from different functions depending on borders
    if (borders == "regional"){
        france_sf <- ne_states(country = "France", returnclass = "sf")
        
    } else if (borders == "national") {
        france_sf <- ne_countries(
            country = "France", scale = "large", returnclass = "sf")    
        
    } else {
        stop(paste0(borders, " is not recognised (should be one of 'regional' or 'national')"))
    }
    
    # Use hard coded extent of france to exclude overseas territories from shp
    st_crop(
        france_sf, 
        xmin = -5.5, 
        xmax = 9.7, 
        ymin = 41.2, 
        ymax = 51.2)
}

# A function that clips values of a raster outside of a shapefile
# ARGS:
#   - raster: a SpatRaster.
#   - shapefile: a sf object.
# RETURNS:
#   - a SpatRaster clipped to sf borders.
clip_raster_from_shapefile <- function(raster, shapefile) {
    # verify raster and shapefile crs are marching
    if (!same.crs(raster, shapefile)) {
        shapefile <- st_transform(shapefile, crs(raster))  # fix silently
    }
    # Clip the raster using the shapefile
    mask(raster, shapefile, touches = TRUE)
}

# A function to transform a raster into a (raster) grid layered on france
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
# RETURNS:
#   - NULL, saves the re-projected SpatRaster silently.
project_to_france_custom_grid <- function(raster, save_to, res_km, verbose = TRUE) {   
    # fit raster to new grid
    if (verbose) {
        cli_alert_info("Fitting new grid...\n")
    }
    france_grid_template <- get_france_raster_template(res_km)
    raster_new_grid <- project(
        raster, france_grid_template, method = "near")

    if (verbose) {
        cli_alert_info("Masking...\n")
    }
    # remove values outside of france borders
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = 2*res_km*1000) |> # buffer distance in meters
        st_transform(crs(raster_new_grid)) # transform to match raster's crs
    raster_new_grid_clean <- clip_raster_from_shapefile(
        raster_new_grid, france_sf_buffered)

    if (verbose) {
        cli_alert_info("Saving file...\n")
    }
    writeRaster(raster_new_grid_clean, save_to, overwrite = TRUE)    
    return(NULL)
}

# A function to transform a raster into a (shapefile) grid of hexagons layered 
# on france metropolitan area.
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The filepath where the shapefile will be saved. 
#       Filepath must end in ".gpkg" or ".shp" to save sf.
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
#   - numeric_stat: a string. 
#       If raster contains numerical values, either "mean" or "median".
#   - lat_min, lat_max: latitudinal extent for the output (in EPSG:4326).
# RETURNS:
#   - NULL, saves the shapefile silently.
project_to_hexagons <- function(
        raster, save_to, res_km, numeric_stat, 
        lat_min = LAT_MIN, lat_max = LAT_MAX, 
        verbose = TRUE) {
    if (!numeric_stat %in% c("mean", "median")) {
        stop(paste0("Unknown numeric_stat: '", numeric_stat, 
            "'. Use 'mean' or 'median'."))
    }    
    
    # Determine, per layer, whether it's categorical
    is_cat <- is.factor(raster)
    layer_names <- names(raster)
    n_layers <- terra::nlyr(raster)

    # 1. Define grid of hexagons
    if (verbose) {
        cli_alert_info("Creating grid of hexagons...\n")
    }
    dggs <- dgconstruct(
        spacing = res_km, metric = TRUE, resround = 'nearest')
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = 2*res_km*1000) |> # buffer distance in meters
        st_transform(crs("EPSG:4326"))     # transform to match raster's crs

    # get mean latitude and longitude resolutions
    lat_mean <- (lat_min + lat_max) / 2
    res_lat <- res_km / 111.0
    res_lon <- res_km / (111.0 * cos(lat_mean * pi / 180))

    full_grid <- dgshptogrid(dggs, france_sf_buffered, 
        cellsize = min(res_lat, res_lon)/2) # ensures no hexagon is forgotten

    # compute areas
    full_grid$hex_area <- st_area(full_grid)
    hex_clipped <- st_intersection(full_grid, france_sf_buffered)
    hex_clipped$clipped_area <- st_area(hex_clipped)

    # filter (keep hexagons with at least 50% area within shape)
    hex_filtered <- hex_clipped[
        as.numeric(hex_clipped$clipped_area / hex_clipped$hex_area) >= 0.5, 
    ]
    hex_grid <- full_grid[full_grid$seqnum %in% hex_filtered$seqnum, ]

    # 2. Associate each hexagon to the most represented value within
    # its area in the raster
    if (verbose) {
        cli_alert_info(paste0("Extracting values from ", n_layers, " layer(s), this might take some time...\n"))
    }

    extracted <- exact_extract(
        raster,
        hex_grid,
        function(values, coverage_fracs) {
            # Normalize to data.frame so single- and multi-layer rasters are 
            # handled the same way
            values <- as.data.frame(values)
            keep <- coverage_fracs >= 0.5

            vapply(seq_len(n_layers), function(i) {
                v <- values[keep, i]
                v <- v[!is.na(v)]
                if (length(v) == 0) return(NA_real_)

                if (is_cat[i]) {
                    # Return most represented value
                    as.numeric(names(sort(table(v), decreasing = TRUE))[1])
                } else if (numeric_stat == "mean") {
                    mean(v)
                } else {
                    median(v)
                }
            }, numeric(1))
        }
    )

    # exact_extract stacks per-feature vectors as columns -> dim = c(n_layers, n_features)
    if (is.null(dim(extracted))) {
        # only happens when n_layers == 1 (FUN returns a scalar, not a vector)
        extracted <- matrix(extracted, ncol = n_layers)
    } else {
        extracted <- t(extracted)
    }

    stopifnot(nrow(extracted) == nrow(hex_grid), ncol(extracted) == n_layers)
    colnames(extracted) <- layer_names

    for (i in seq_len(n_layers)) {
        hex_grid[[layer_names[i]]] <- extracted[, i]
    }

    # 3. Polish and save
    if (verbose) {
        cli_alert_info("Saving files...\n")
    }
    st_write(hex_grid, save_to, delete_dsn = TRUE, quiet = TRUE)

    # Export legend for each categorical layer
    if (any(is_cat)) {
        cats_list <- cats(raster)
        for (i in which(is_cat)) {
            category_table <- cats_list[[i]]
            if (is.null(category_table)) next
            legend <- category_table %>%
                mutate(hex = rgb(Red, Green, Blue, maxColorValue = 1))
            out_csv <- paste0(
                file_path_sans_ext(save_to), "_", layer_names[i], ".csv")
            write_csv(legend, out_csv)
        }
    }
    return(NULL)
}

# A function that collapses categories of land cover for the original CORINE 
# SpatRaster
# WARNING:
#   - This function will only work correctly with the ORIGINAL CLC SpatRaster 
#       (downloaded from https://www.data.gouv.fr/datasets/corine-land-cover-edition-2018-france-metropolitaine).
# ARGS:
#   - clc_raster: The original CLC SpatRaster.
#   - level_*: an integer. For each category, it decides if the subcategories
#       should be collapsed together at level 1, 2 or 3 
#       (default is 3: all subcategories are kept as they are).
#   - verbose: a boolean. If TRUE (default), shows info messages.
# RETURNS:
#   - the updated CLC SpatRaster
simplify_CLC <- function(
        clc_raster,
        level_urban = 1,
        level_crops = 1,
        level_forests = 1,
        level_wetlands = 1,
        level_water = 1,
        verbose = TRUE){
    if (verbose) {
        cli_alert_warning("[Harmless warning] This function should only be used with CORINE Land Cover rasters!")
        cli_alert_info("Loading categories...")
    }

    # Get dataframe
    cat_table <- cats(clc_raster)[[1]]   

    # Only subcategories are named in raster. Adding other names.
    all_names <- setNames(cat_table$LABEL3, cat_table$CODE_18)
    all_names <- c(
        all_names,
        "1" = "Artificial surfaces",
        "11" = "Urban fabric",
        "12" = "Industrial, commercial and transport units",
        "13" = "Mine, dump and construction sites",
        "14" = "Artificial, non-agricultural vegetated areas",
        "2" = "Agricultural areas",
        "21" = "Arable land",
        "22" = "Permanent crops",
        "23" = "Pastures",
        "24" = "Heterogeneous agricultural areas",
        "3" = "Forest and semi-natural areas",
        "31" = "Forests",
        "32" = "Shrub and/or herbaceous vegetation associations",
        "33" = "Open spaces with little or no vegetation",
        "4" = "Wetlants",
        "41" = "Inland wetlands",
        "42" = "Coastal wetlans",
        "5" = "Water bodies",
        "51" = "Inland waters",
        "52" = "Marine waters",
        "999" = "NO_DATA"
    )

    # List for levels of collapse
    levels <- list(
        as.integer(substr(as.character(cat_table$CODE_18), 1, 1)),
        as.integer(substr(as.character(cat_table$CODE_18), 1, 2)),
        as.integer(substr(as.character(cat_table$CODE_18), 1, 3))
    )
    groups <- c(level_urban, level_crops, level_forests, 
                level_wetlands, level_water)

    if (verbose) {
        cli_alert_info("Creating new table...")
    }
    # Create new table
    cat_table <- cat_table |>
        mutate(NEW_CODE_18 = case_when(
            grepl("^1", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 1 ]] ]]),
            grepl("^2", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 2 ]] ]]),
            grepl("^3", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 3 ]] ]]),
            grepl("^4", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 4 ]] ]]),
            grepl("^5", cat_table$CODE_18) ~ as.character(levels[[ groups[[ 5 ]] ]]),
            # Default case (if none of the above match)
            TRUE ~ "999")) |>
        group_by(NEW_CODE_18) |>
        mutate(NEW_LABEL3 = all_names[NEW_CODE_18]) |>
        mutate(
            NEW_Red = mean(Red, na.rm = TRUE),
            NEW_Green = mean(Green, na.rm = TRUE),
            NEW_Blue = mean(Blue, na.rm = TRUE)) |>
        ungroup()    
    
    # Build a reclassification matrix: from Value → new integer code
    rcl <- cat_table |>
        select(Value, NEW_CODE_18) |>
        distinct() |>
        as.matrix()

    clc_raster_collapsed <- classify(clc_raster, rcl)

    # An assign new table table keyed on the new codes
    new_levels_simple <- cat_table |>
        select(NEW_CODE_18, NEW_LABEL3, NEW_Red, NEW_Green, NEW_Blue) |>
        distinct() |>
        mutate(
            Value = as.integer(NEW_CODE_18),
            Red = NEW_Red,
            Green = NEW_Green,
            Blue = NEW_Blue) |>
        select(Value, NEW_LABEL3, NEW_CODE_18, Red, Green, Blue)

    levels(clc_raster_collapsed) <- new_levels_simple 

    if (verbose) {
        cli_alert_info("Assigning new colors...")
    }
    # Build back the color table
    new_coltab <- cat_table |>
        select(NEW_CODE_18, NEW_Red, NEW_Green, NEW_Blue) |>
        distinct() |>
        mutate(
            Value = as.integer(NEW_CODE_18),
            R = round(NEW_Red * 255),
            G = round(NEW_Green * 255),
            B = round(NEW_Blue * 255),
            A = 255) |>
        select(Value, R, G, B, A) |>
        as.data.frame()
    coltab(clc_raster_collapsed) <- new_coltab

    # Show layer of names first
    activeCat(clc_raster_collapsed) <- "NEW_LABEL3" 

    if (verbose) {
        cli_alert_success("Re-classified CLC2018 is ready!")
    }
    return(clc_raster_collapsed)
}

# A function that transforms a layer of categorical values from a raster
# into as many layers as there were categories in the original layer. Each
# new layer contains a raster with each cell containing the distance to the 
# closest category in the original raster.
# ARGS:
#   - raster_layer: a SpatRaster containing a single layer.
# RETURNS:
#   - a SpatRaster with one layer of distances per category in raster_layer.
factorial_to_distance2factor <- function(raster_layer) {
    # raster reprojection for faster computation
    original_template <- raster_layer
    r_proj <- project(raster_layer, "EPSG:3035", method = "near")

    # Get the category labels
    cat_df <- levels(r_proj)[[1]]
    cat_labels <- cat_df[[2]]
    cat_codes  <- cat_df[[1]]

    # 2. Create a mask layer for each category 
    # (1 = category present, NA = everything else)
    bin_stack <- segregate(
        r_proj, classes = cat_codes, keep = TRUE, other = NA)

    # rename to the readable labels
    names(bin_stack) <- cat_labels

    # For each binary layer, compute distance to nearest non-NA cell
    n <- nlyr(bin_stack)
    dist_list <- vector("list", n)

    for (i in seq_len(n)) {
        lyr_name <- names(bin_stack)[i]
        d <- distance(bin_stack[[i]])
        names(d) <- paste0("distance to ", lyr_name)
        dist_list[[i]] <- d
    }

    # Combine into one multi-layer SpatRaster
    dist_stack <- rast(dist_list)

    # raster reprojection to original crs
    dist_stack <- project(dist_stack, original_template, method = "bilinear")
    return(dist_stack)
}

# A function to automatically simplify CLC categories and reduce its extent.
# ARGS:
#   - raster_path: a filepaths to a raster (.tif) file.
#   - save_to: a filepath to save the new raster (must end with .tif).
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in EPSG:4326).
# RETURNS:
#   - a SpatRaster with several layers
#       - first layer:  
#       - next layers:
save_simplified_clc <- function(
        raster_path,
        save_to,
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
    
    # load and simplify CORINE raster
    cli_alert_info("Simplification of CLC2018 raster.")
    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning("Skipping simplification of CLC2018.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        raster <- rast(raster_path)

        # crop raster for faster computation
        reduced_extent <- ext(
                lon_min, lon_max, 
                lat_min, lat_max) |>
            as.polygons(crs = "EPSG:4326") |>
            densify(interval = 0.1) |>
            project("EPSG:3035")

        cropped_raster <- crop(raster, reduced_extent)


        # CLC contains (too) many categories. Subcategories can be grouped
        # together to avoid having too many modalities.
        # We keep all categories at level 1 (i.e. lowest information) ...
        # ... except crops, for which we get level 2 subcategories 
        # (level 3 is also available)
        cli_alert_info("Simplifying categories...")
        simplified_clc_raster <- simplify_CLC(
            clc_raster = cropped_raster,
            level_urban = 1,
            level_crops = 2,
            level_forests = 1,
            level_wetlands = 1,
            level_water = 1)

        # Once the categories are set, we need to convert this layer into a 
        # set of numerical layers : each corresponding to the distance
        # to the closest category in the "simplified_clc_raster"
        cli_alert_info("Adding one layer per categories (distance layers)...")
        distance_layers <- factorial_to_distance2factor(simplified_clc_raster)

        clc_stacked_raster <- rast(list(simplified_clc_raster, distance_layers))
        
        cli_alert_info("Saving file...")

        writeRaster(
            project(clc_stacked_raster, "EPSG:4326"),  
            save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created simplified CORINE land cover raster file!\n\n"))
        
        return(clc_stacked_raster)
    }
}

# A function to clip a raster to france metropolitan area and ensure crs WGS84
# ARGS:
#   - raster: a SpatRaster.
#   - save_to: a string or NULL. 
#       If NULL, returns the raster after processing. 
#       If not, the string given is the path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - verbose: a boolean. If TRUE, shows info messages.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
# RETURNS:
#   - clipped SpatRaster
clip_raster_france_wgs84_crs <- function(
        raster,  
        buffer, 
        verbose = TRUE,
        save_to = NULL) {
    
    # import france shapefile with buffer to avoid clipping important data
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(2154) |>              # EPSG:2154 = RGF93, better conservation of distances
        st_buffer(dist = buffer*1000) |> # buffer distance in meters
        st_transform(crs(raster))   # transform to match raster's crs
    
    # reduce size of raster
    raster_cropped <- crop(raster, ext(france_sf_buffered))
    raster_clipped <- clip_raster_from_shapefile(
        raster_cropped, france_sf_buffered)
    
    if (verbose) {
        cli_alert_info("Updating CRS, this might take some time...\n")
    }
    raster_wgs84 <- project(raster_clipped, "EPSG:4326")

    if (verbose) {
        cli_alert_info("Saving file...\n")
    }

    if (is.null(save_to)) {
        return(raster_wgs84)
    } else {
        writeRaster(raster_wgs84, save_to, overwrite = TRUE)  
    }   
}


# A function to standardize raster projection and save path.
# ARGS:
#   - raster_path: filepaths to a raster (.tif) file.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
# RETURNS:
#   - The reprojected SpatRaster
standardized_raster_projection <- function(
        raster_path,
        save_to_basename,
        res_km = RES_KM) {
    cli_alert_info(paste0("Pre-processing ", basename(raster_path), "..."))

    # Creating the paths that will be needed in the function
    raster_reprojected_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_projection_france_res", res_km,"km-WGS84",
        ".tif")
        
    if (!authorise_overwrite(raster_reprojected_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing of raster.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        # import and reproject to a standardized grid
        raster <- rast(raster_path)

        . <- project_to_france_custom_grid(
                raster = raster, 
                save_to = raster_reprojected_path,  
                res_km = res_km) 
        cli_alert_success("Modified raster file saved!\n\n")
    }

    return(raster_reprojected_path)
}

# A function to standardize hexagonal projection and save path.
# ARGS:
#   - raster_path: filepaths to a raster (.tif) file.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - extraction_mode: a string. How to handle numeric data aggregation,
#       see 'project_to_hexagons' function.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km, approximative).
# RETURNS:
#   - The reprojected SpatRaster
standardized_hexagonal_projection <- function(
        raster_path,
        save_to_basename,
        extraction_mode,
        res_km = RES_KM) {
    cli_alert_info(paste0("Pre-processing of ", basename(raster_path), "..."))

    # Creating the paths that will be needed in the function
    shapefile_hexagons_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_projection_france_hexagons_res", res_km,"km-WGS84",
        ".gpkg")


    if (!authorise_overwrite(shapefile_hexagons_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        # import and reproject to a standardized grid of hexagons
        raster <- rast(raster_path)
        . <- project_to_hexagons(
                raster = raster, 
                save_to = shapefile_hexagons_path, 
                res_km = res_km,
                numeric_stat = extraction_mode) 
        cli_alert_success("Saved shapefile and CSV legend!\n\n")
    }

    return(shapefile_hexagons_path)
}