# Set of utility functions used accross scripts.

##### Liraries ##### ----------------------------------------------------------
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


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # all parameters are grouped together


##### Global functions ##### --------------------------------------------------
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


##### Datasets ##### ----------------------------------------------------------
# A function that returns the project extent (in PROJECT_CRS units) as a polygon
# expressed in any CRS. Used to crop rasters that are stored in another CRS.
# ARGS:
#   - crs_out: the CRS of the returned polygon (default: PROJECT_CRS).
#   - xmin, xmax, ymin, ymax: extent, in PROJECT_CRS units (meters).
# RETURNS:
#   - a SpatVector (polygon).
get_extent_polygon <- function(
        crs_out = PROJECT_CRS,
        xmin = LON_MIN, xmax = LON_MAX, 
        ymin = LAT_MIN, ymax = LAT_MAX) {
    as.polygons(ext(xmin, xmax, ymin, ymax), crs = PROJECT_CRS) |>
        densify(interval = 10000) |>   # 10 km (PROJECT_CRS is in meters)
        project(crs_out)
}

# A helper that makes sure a categorical raster is in PROJECT_CRS and written 
# on disk as integers (CORINE codes are small ints). If the raster is already in
# PROJECT_CRS it is only written to disk (no resampling, so no distortion).
# ARGS:
#   - raster: a SpatRaster (single categorical layer).
#   - tmp_dir: a temporary directory to save the intermediary raster.
# RETURNS:
#   - a SpatRaster (file-backed) in PROJECT_CRS.
to_project_crs_on_disk <- function(raster, tmp_dir = tempdir()) {
    filename <- file.path(tmp_dir, "r_proj.tif")
    if (same.crs(raster, PROJECT_CRS)) {
        writeRaster(raster, filename, datatype = "INT2U", overwrite = TRUE)
    } else {
        project(raster, PROJECT_CRS, method = "near",
                filename = filename, datatype = "INT2U", overwrite = TRUE)
    }
}

# A function that "blurs" coordinate within a dataframe by randomly shifting
#   each longitude and latitude coordinates.
# ARGS:
#   - df: a dataframe, sf object or SpatVector with gps coordinates.
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
        x_lon, y_lat, blur_km, 
        seed = NULL, 
        lon_min = LON_MIN, lon_max = LON_MAX, 
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
    # for reproducible results  
    if (!is.null(seed)){
        set.seed(seed) 
    }

    # Check that the dataframe is in EPSG:4326 coordinates (approx.)
    if (!(all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[x_lon])) &
        all(as.data.frame(df[x_lon]) <= FRANCE_BBOX_4326["xmax"])   &
        all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[y_lat]))   &
        all(as.data.frame(df[y_lat]) <= FRANCE_BBOX_4326["ymax"])    
    )) {
        stop(paste0(
            "save_features_from_obs() assumes that your data is in",
            "EPSG:4326 over France's extent. It looks like 'df' is in ",
            "another CRS, or on a different area. Please check your data."
        ))
    }


    # if sf object, remove geometry field (we're gonna make a new one)
    is_sf <- inherits(df, "sf")
    is_sv <- inherits(df, "SpatVector")
    if (is_sf) {
        crs_orig <- sf::st_crs(df)
        geom_col <- attr(df, "sf_column")
        df <- sf::st_drop_geometry(df) 
    } else if (is_sv) {
        crs_orig <- terra::crs(df)
        df <- as.data.frame(df)   # attributes only, geometry dropped
    }

    # Move to a metric CRS: the blur is then exactly in meters, whatever the latitude
    xy_m <- data.frame(x = df[[x_lon]], y = df[[y_lat]]) |>
        st_as_sf(coords = c("x", "y"), crs = 4326) |>
        st_transform(PROJECT_CRS) |>
        st_coordinates()

    # Dividing by 2.45 so ~95% of points stay within a circle of blur_km
    # FYI this is the value for 2D, in 1D it is the usual 1.96
    sd_m <- blur_km * 1000 / 2.45
    xy_m <- xy_m + matrix(rnorm(2 * nrow(xy_m), mean = 0, sd = sd_m), ncol = 2)

    # Back to lon/lat
    xy_deg <- data.frame(x = xy_m[, 1], y = xy_m[, 2]) |>
        st_as_sf(coords = c("x", "y"), crs = PROJECT_CRS) |>
        st_transform(4326) |>
        st_coordinates()
    df[[x_lon]] <- xy_deg[, 1]
    df[[y_lat]] <- xy_deg[, 2]
    
    # make a new "geometry" column
    if (is_sf) {
        df <- st_as_sf(
            df, coords = c(x_lon, y_lat), crs = crs_orig, remove = FALSE)
        names(df)[names(df) == "geometry"] <- geom_col
        st_geometry(df) <- geom_col
    } else if (is_sv) {
        df <- terra::vect(
            df,
            geom = c(x_lon, y_lat),
            crs = crs_orig,
            keepgeom = TRUE   # keeps the lon/lat columns as attributes, like remove = FALSE
        )
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
    for (n_path in 1:12){
        # import raster
        raster <- rast(raster_paths[n_path])
               
        cropped_raster <- clip_raster_france(
            raster, verbose = verbose, buffer = buffer)     

        if (n_path == 1){
            rasters_together <- cropped_raster
        } else {
            rasters_together <- c(rasters_together, cropped_raster)
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
    cli_alert_info("Aggregation of monthly data...")
    
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
        writeRaster(annual_raster, save_to, overwrite = TRUE)
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
    x_lon <- "LON"
    y_lat <- "LAT"
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

        # Check that the dataframe is in EPSG:4326 coordinates (approx.)
        if (!(all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[x_lon])) &
            all(as.data.frame(df[x_lon]) <= FRANCE_BBOX_4326["xmax"])   &
            all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[y_lat]))   &
            all(as.data.frame(df[y_lat]) <= FRANCE_BBOX_4326["ymax"])    
        )) {
            stop(paste0(
                "save_features_from_obs() assumes that your data is in",
                "EPSG:4326 over France's extent. It looks like 'df' is in ",
                "another CRS, or on a different area. Please check your data."
            ))
        }

        # convert dataframe to geo dataframe
        cli_alert_info(paste0("Conversion to shapefile..."))
        points <- vect(st_as_sf(
            df, coords = c(x_lon, y_lat), crs = 4326, remove = FALSE))

        # the user can choose which type of file the features must be extracted
        # from. If raw, uses the initial raw resolution of the files (which are
        # usually finer than 1km² in resolution). If using raster or hexagonal,
        # uses pre-processed data that was made to be pixels/hexagons space out
        # by RES_KM.
        if (mode == "raw") {
            environmental_data_paths <- c(
                CORINE_PATH_PREPROCESSED, CHELSA_ANNUAL_PATHS, NDVI_ANNUAL_PATH,
                ELEVATION_PATH_PREPROCESSED, LIGHT_POLLUTION_FILEPATH, SWF_PATH_PREPROCESSED)
            
        } else if (mode == "raster") {
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH, SWF_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_res10km-", PROJECT_CRS_TAG, ".tif")
            
        } else if (mode == "hexagonal") {
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH, SWF_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_hexagons_res10km-", PROJECT_CRS_TAG, ".gpkg")
            
        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }

        # new names are needed when adding columns to 'points' (geo dataframe) 
        environmental_data_names <- c(
            "clc", paste0("chelsa_", CHELSA_DATASETS), "ndvi", "elevation",
            "light_pollution", "swf")   
        
        cli_bullets(c(
            "!" = "[Fetching features]: Using hard coded paths:",
            setNames(environmental_data_paths, rep("*", length(environmental_data_paths)))
        ))

        for (i in 1:length(environmental_data_paths)) {
            
            # Step 1: extract the values of each layer at each point, as a 
            # dataframe with one row per point and one column per layer
            if (mode %in% c("raw", "raster")) {
                # load envrionmental data
                env_data <- rast(environmental_data_paths[i])
                layer_ids <- names(env_data)
                
                # `points` stay in EPSG:4326; use a copy in the raster's CRS for extraction
                pts_i <- project(points, crs(env_data))

                # drop the leading ID column, keep one column per layer
                extracted_values <- terra::extract(
                    env_data, pts_i)[, -1, drop = FALSE]

            } else if (mode == "hexagonal") {
                # hexagons have one attribute per layer of the raster they 
                # were made from (see project_to_hexagons)
                env_data <- vect(environmental_data_paths[i])

                # extract values of points projected onto env_data cells
                extracted_values <- terra::extract(
                    env_data, project(points, crs(env_data)))

                if (nrow(extracted_values) != nrow(points)) {
                    stop("Overlaping shapes in .gpkg file.")
                }
                # drop the leading ID column, keep one column per layer
                extracted_values <- extracted_values[, -1, drop = FALSE]
                layer_ids <- colnames(extracted_values)
            }

            # Step 2 (same for all modes): add one column to `points` per layer.
            # A single layer is named after environmental_data_names, many 
            # layers keep their own names.
            if (ncol(extracted_values) == 1) {
                col_names <- environmental_data_names[i]
            } else {
                col_names <- layer_ids
            }
            colnames(extracted_values) <- col_names

            for (col in col_names) {
                points[[col]] <- extracted_values[[col]]
            }
        }

        if (blur) {
            # replace raw coordinates with blurred coordinates
            cli_alert_info("Data anonymization...")
            source(here::here(file.path("data", "config", "seed.R"))) 
            anonym_df <- blur_coordinates(
                points, "LON", "LAT", BLUR_SIZE, BLUR_SEED)
        } else {
            # Rebuild geometry explicitly from the LON/LAT columns
            coords_df <- as.data.frame(points)

            anonym_df <- vect(
                coords_df,
                geom = c("LON", "LAT"),
                crs = "EPSG:4326",
                keepgeom = TRUE
            )
        }

        # Save resulting dataframe (reproject in target CRS)
        cli_alert_info("Saving file...")
        writeVector(project(anonym_df, PROJECT_CRS), save_to, overwrite = TRUE)
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
                    ELEVATION_PATH_PREPROCESSED, LIGHT_POLLUTION_FILEPATH, SWF_PATH_PREPROCESSED)
                
            } else { # mode == "raster"
                environmental_data_paths <- c(
                    CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                    NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                    LIGHT_POLLUTION_PROJECTION_BASEPATH, SWF_PROJECTION_BASEPATH)  
                environmental_data_paths <- paste0(
                    environmental_data_paths, 
                    "_projection_france_res", RES_KM, "km-", PROJECT_CRS_TAG, ".tif")
            }
        } else if (mode == "hexagonal") {
            # Shapefiles must be handled differently than rasters
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH, SWF_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_projection_france_hexagons_res", RES_KM, "km-", PROJECT_CRS_TAG, ".gpkg")
        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }

        cli_bullets(c(
            "!" = "[Fetching features]: Using hard coded paths:",
            setNames(environmental_data_paths, rep("*", length(environmental_data_paths)))
        ))

        # Merge them together
        if (mode == "raw" | mode == "raster") {            
            target_ext <- ext(lon_min, lon_max, lat_min, lat_max)

            # Resolutions in meters (degrees converted via the N-S length of a degree)
            resolutions <- sapply(environmental_data_paths, function(p) {
                r <- rast(p)
                if (is.lonlat(r)) res(r)[2] * 111320 else res(r)[1]
            })
            target_res <- ceiling(max(resolutions) / 100) * 100   # round up to next 100 m
            template <- rast(target_ext, resolution = target_res, crs = PROJECT_CRS)
            cli_alert_info(paste0(
                "Selecting raster with coarsest resolution (res: ", 
                round(target_res, 2), " m)."))
           
    
            # Crop and resample each raster onto the template
            cli_alert_info("Resampling rasters, can take several minutes...")
            raster_list <- list()
            for (path in environmental_data_paths) {
                env_raster <- rast(path)

                # cheap pre-crop
                if (!same.crs(env_raster, PROJECT_CRS)) {
                    ext_native <- project(target_ext, from = PROJECT_CRS, to = crs(env_raster))
                    env_raster <- crop(env_raster, ext_native * 1.15)  # small buffer
                } else {
                    env_raster <- crop(env_raster, target_ext)
                }

                # categorical vs continuous resampling method
                if (grepl("CLC", path, ignore.case = TRUE)) {
                    method <- "near"
                } else {
                    method <- "bilinear"
                }

                # takes care of the projection as well
                raster_resampled <- terra::project(
                    env_raster, template, method = method)

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
                } else if (grepl("SWF", path, ignore.case = TRUE)) {
                    name <- "swf"
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
            master_raster <- fix_names(rast(raster_list))

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

                    if (grepl("CLC", path, ignore.case = TRUE)) {
                        # Categorical layer: keep every attribute, no prefix
                        df <- df[, c(join_col, setdiff(names(df), join_col)), drop = FALSE]

                    } else {
                        # Numeric layers: name the output column after the variable
                        name <- if (grepl("CHELSA", path, ignore.case = TRUE)) {
                            paste(str_split(basename(path), "_")[[1]][c(1, 3)], collapse = "_")
                        } else if (grepl("light", path, ignore.case = TRUE)) "light_pollution"
                        else if (grepl("Elevation", path, ignore.case = TRUE)) "Elevation"
                        else if (grepl("ndvi", path, ignore.case = TRUE)) "ndvi"
                        else stop("Unexpected layer name, modify the function.")

                        # drop the id and hex_area, what remains is the value column
                        value_cols <- setdiff(names(df), c(join_col, "hex_area"))
                        value_cols <- value_cols[!grepl("hex_area$", value_cols)]

                        if (length(value_cols) != 1) {
                            stop(paste0(
                                "Expected exactly one value column in ", basename(path),
                                " but found: ", paste(value_cols, collapse = ", ")))
                        }

                        df <- df[, c(join_col, value_cols), drop = FALSE]
                        names(df)[2] <- name
                    }

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

            master_df <- Reduce(
                function(a, b) merge(a, b, by = join_col),
                lapply(hex_list, as.data.frame))
            master_hex <- hex_list[[1]]
            values(master_hex) <- master_df[match(master_hex$seqnum, master_df$seqnum), ]

            cli_alert_success("Saving hexagon layer...")
            writeVector(master_hex, save_to, overwrite = TRUE)
            cli_alert_success("Saved stacked hexagons!\n\n")
            return(fix_names(master_hex))

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
#   - keep.lyr: a vector of strings, the names of the layers to keep. 
#               Default is NULL (returns all layers).
# RETURNS:
#   - a multi-layer SpatRaster.
get_master_t_rasters <- function(mode = "raw", keep.lyr = NULL) {
    # Load master raster
    if (mode == "raw") {
        cli_alert_info("Loading multi-layer raster with best resolution.")
        master_raster <- rast(FEATURES_RAW_TRF)
    } else if (mode == "res") {
        cli_alert_info(paste0(
            "Loading multi-layer raster (res = ", RES_KM, "km)."))
        master_raster <- rast(FEATURES_RES_TRF)
    } else {
        stop("The 'mode' parameter should be one of: 'raw', 'res'.")
    }
    
    # format names to common format
    master_raster <- fix_names(master_raster)
    
    # Keep only the necessary layers
    if (!is.null(keep.lyr)) {
        master_raster <- master_raster[[keep.lyr]]
    }
    
    # reduce size for faster computation
    master_raster <- clip_raster_france(  
            master_raster,  
            buffer = 5*RES_KM,   
            verbose = FALSE,
            save_to = NULL)
 
    return(master_raster)
}


##### Rasters/SF functions ##### ----------------------------------------------
# A function that checks whether a SpatRaster or SpatVector has a metric CRS.
# Inspired from terra::is.lonlat
# ARGS:
#   - r: a SpatRaster or SpatVector
# RETURNS:
#   - a boolean
is_metric <- function(r) {
    crs_text <- crs(r)
    !is.lonlat(r) &&
        grepl('UNIT\\["metre"|LENGTHUNIT\\["metre"',
            crs_text, ignore.case = TRUE)
}

# A function that replaces NA values within a buffer around actual values in a
# raster file.
# ARGS:
#   - raster_path: a filepath to a raster (.tif) file.
#   - save_to: a filepath to save the new raster (must end with .tif).
#   - replace_with: the replacement value for selected NAs
#   - buffer_km: a numeric, the size of the buffer (in km) which selects NAs.
# RETURNS:
#   - The buffered raster
add_buffer_to_raster <- function(
    raster_path, save_to, replace_with = 0, buffer_km = RES_KM*5) {
    # load raster
    raster <- rast(raster_path)
    if (!is_metric(raster)) {
        stop("This function can only be used with rasters in metric CRS.")
    }

    # compute distance to closest not NA cell
    d <- distance(raster)

    # Create a buffer where NAs are replaced by "replace_with"
    raster_with_buffer <- ifel(
        (is.na(raster) & (d <= (buffer_km*1000))), 
        replace_with, 
        raster)
    writeRaster(raster_with_buffer, save_to, overwrite = TRUE)    
    
    return(raster_with_buffer)
}

# A function that creates an empty raster that matches France's extent (WGS 84)
# ARGS:
#   - res_km: an integer/a float. 
#       The size for each pixel of the raster (in km, approximative).
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
# RETURNS:
#   - a SpatRaster with empty values to the extent of France's metropolitan
#       area.
get_france_raster_template <- function(
        res_km, lon_min = LON_MIN, lon_max = LON_MAX,
        lat_min = LAT_MIN, lat_max = LAT_MAX) {
    rast(extent = ext(lon_min, lon_max, lat_min, lat_max),
         resolution = res_km * 1000, crs = PROJECT_CRS)
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
    france_sf <- st_set_agr(france_sf, "constant")
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
    # NOTE: method = "near" keeps a single source pixel per cell. When the 
    # source is much finer than res_km, "average" (continuous layers) or 
    # "mode" (categorical layers) would be more representative.
    raster_new_grid <- project(
        raster, france_grid_template, method = "near")

    if (verbose) {
        cli_alert_info("Masking...\n")
    }
    # remove values outside of france borders
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
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
#   - lat_min, lat_max: latitudinal extent for the output (in PROJECT_CRS).
# RETURNS:
#   - NULL, saves the shapefile silently.
project_to_hexagons <- function(
        raster, save_to, res_km, numeric_stat, 
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
        st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
        st_buffer(dist = 2*res_km*1000) |> # buffer distance in meters
        st_transform(4326)                 # dggridR works only with lon/lat

    # get mean latitude and longitude resolutions
    lat_mean <- mean(FRANCE_BBOX_4326[c("ymin", "ymax")])  # degrees
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

    # Polygons must be in the raster's CRS for the extraction
    hex_for_extract <- st_transform(hex_grid, crs(raster))

    extracted <- exact_extract(
        raster,
        hex_for_extract,
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
    hex_grid <- st_transform(hex_grid, PROJECT_CRS)  # saved in PROJECT_CRS
    st_write(hex_grid, save_to, delete_dsn = TRUE, quiet = TRUE)

    # Export legend for each categorical layer
    if (any(is_cat)) {
        cats_list <- cats(raster)
        for (i in which(is_cat)) {
            category_table <- cats_list[[i]]
            if (is.null(category_table)) next
            legend <- category_table %>%
                mutate(hex = rgb(Red, Green, Blue, maxColorValue = 1))
            print(paste("LAYER NAME:", layer_names[i]))
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
        "4" = "Wetlands",
        "41" = "Inland wetlands",
        "42" = "Coastal wetlands",
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

    # Set NO_DATA (code 999) to NA
    clc_raster_collapsed[clc_raster_collapsed == 999] <- NA

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

# A function that uses a raster layer of categorical values to generate
# sub-layers : each sub-layer shows the distance of each cell to the closest 
# category in the input raster layer. 
# ARGS:
#   - raster: a SpatRaster (single layer) of categorical values.
#   - temp_dir: a temporary directory to save intermediary rasters
# RETURNS:
#   - a SpatRaster with one layer of distances per category in raster
factorial_to_distances <- function(
        raster,
        tmp_dir = tempdir()) {

    # Tell terra to be conservative with RAM and use of disk
    terraOptions(memfrac = 0.5, todisk = TRUE, tempdir = tmp_dir, progress = 1)

    raster_clipped <- clip_raster_france(
        raster, buffer = 5*RES_KM, verbose = FALSE, save_to = NULL)

    # Write projected raster to disk, as integers (CORINE codes are small ints)
    r_proj <- to_project_crs_on_disk(raster_clipped, tmp_dir)
    rm(raster_clipped); gc()

    present <- freq(r_proj)$value
    present <- present[!is.na(present)]
    n <- length(present)

    # Process each class, writing each distance layer to its own file
    dist_files <- character(n)
    for (i in seq_len(n)) {
        lab <- present[i]
        cli_alert_info("Processing distances for category {i} of {n}")

        m <- ifel(r_proj == lab, 1, NA,
                  filename = file.path(tmp_dir, sprintf("m_%s.tif", lab)),
                  datatype = "INT1U", overwrite = TRUE)

        dist_files[i] <- file.path(tmp_dir, sprintf("dist_%s.tif", lab))
        d <- distance(m, filename = dist_files[i],
                      datatype = "FLT4S", overwrite = TRUE)
        names(d) <- paste0("distance to ", lab)

        rm(m, d); gc()
        unlink(file.path(tmp_dir, sprintf("m_%s.tif", lab)))
    }

    # Virtual stack: reads lazily from disk, no RAM cost
    dist_stack <- rast(dist_files)
    names(dist_stack) <- paste0("distance to ", present)

    # Project back to the original grid, writing to disk
    out <- project(dist_stack, raster, method = "bilinear",
                   filename = file.path(tmp_dir, "dist_final.tif"),
                   datatype = "FLT4S", overwrite = TRUE)
    out
}

# A function that uses a raster layer of categorical values to generate
# sub-layers : each sub-layer shows the proportion of cells surrounding each 
# cell that belong to a category of the input raster layer.
#
# The proportion of a category within a circle is a *convolution* of the 
# category's 0/1 indicator layer with a circular kernel (way faster than 
# terra::focal). 
# ARGS:
#   - raster: a SpatRaster (single layer) of categorical values.
#   - radius_m: a numeric. The radius in meters count categories in.
#   - tmp_dir: a temporary directory to save intermediary rasters
#   - max_fft_cells: a numeric. Maximum number of cells (approximative) of the
#       complex arrays handled at once by the FFT. Lower it if you run out of 
#       RAM (one complex cell = 16 bytes, and ~5 such arrays are alive at once).
# RETURNS:
#   - a SpatRaster with one layer of proportion per category in raster_layer.
factorial_to_proportion <- function(
        raster,
        radius_m,
        tmp_dir = tempdir(),
        max_fft_cells = 8e6) {

    # Tell terra to be conservative with RAM and use of disk
    terraOptions(memfrac = 0.5, todisk = TRUE, tempdir = tmp_dir, progress = 1)

    # Clip with a buffer at least as large as the radius, to avoid edge effects
    # (assumes `buffer` is in km, as in your original call)
    raster_clipped <- clip_raster_france(
        raster,
        buffer  = 5 * RES_KM,
        verbose = FALSE, save_to = NULL)

    # Project to a metric CRS so the radius is in meters
    r_proj <- to_project_crs_on_disk(raster_clipped, tmp_dir)
    rm(raster_clipped); gc()

    # Pixels hold numeric CODES (for categorical rasters, freq() would return 
    # LABELS instead, which cannot be compared to pixel values). Codes are 
    # therefore read directly from the pixels, and labels (when the raster has 
    # a category table) are only used to name the layers.
    n_row <- nrow(r_proj)
    n_col <- ncol(r_proj)
    read_codes <- function(first_row, n_read) {
        v <- readValues(r_proj, row = first_row, nrows = n_read)
        if (!is.numeric(v)) {   # some terra versions can return labels
            cat_table <- cats(r_proj)[[1]]
            v <- cat_table[[1]][match(as.character(v), 
                as.character(cat_table[[activeCat(r_proj) + 1]]))]
        }
        matrix(v, nrow = n_read, ncol = n_col, byrow = TRUE)
    }
    readStart(r_proj)   # r_proj is file-backed: it must be opened to read blocks
    on.exit(try(readStop(r_proj), silent = TRUE), add = TRUE)

    present_codes <- numeric(0)
    for (first_row in seq(1, n_row, by = 500)) {
        n_read <- min(500, n_row - first_row + 1)
        present_codes <- union(
            present_codes, unique(as.vector(read_codes(first_row, n_read))))
    }
    present_codes <- sort(present_codes[!is.na(present_codes)])
    n <- length(present_codes)
    if (n == 0) {
        stop("The raster contains no data (only NAs) after clipping.")
    }

    # Names of the layers: category label if available, else the code
    present <- present_codes
    if (is.factor(r_proj)) {
        cat_table <- cats(r_proj)[[1]]
        labels <- cat_table[[activeCat(r_proj) + 1]][
            match(present_codes, cat_table[[1]])]
        if (anyNA(labels)) {
            cli_alert_warning(paste0(
                "Pixel codes not found in the category table (",
                paste(present_codes[is.na(labels)], collapse = ", "), 
                "): using the codes as names."))
            labels[is.na(labels)] <- present_codes[is.na(labels)]
        }
        present <- labels
    }

    # Circular window: 1 inside the circle, 0 outside (focalMat is used so the 
    # window is exactly the one terra would use).
    kernel <- (focalMat(r_proj, d = radius_m, type = "circle") > 0) * 1
    half_r <- (nrow(kernel) - 1) / 2
    half_c <- (ncol(kernel) - 1) / 2

    # Blocks of rows: each block is read with `half_r` extra rows above and 
    # below (needed to compute the proportions of the block's edge rows).
    pad_c <- nextn(n_col + half_c, factors = c(2, 3, 5))  # fast FFT sizes
    block_rows <- max(1, floor(max_fft_cells / pad_c) - 2 * half_r)
    block_starts <- seq(1, n_row, by = block_rows)

    # FFT of the zero-padded kernel, with its centre moved to the corner [1, 1]
    # (circular convention). Cached because most blocks share the same size.
    kernel_fft_cache <- list()
    get_kernel_fft <- function(pad_r) {
        key <- as.character(pad_r)
        if (is.null(kernel_fft_cache[[key]])) {
            k_pad <- matrix(0, pad_r, pad_c)
            row_idx <- ((seq_len(nrow(kernel)) - 1 - half_r) %% pad_r) + 1
            col_idx <- ((seq_len(ncol(kernel)) - 1 - half_c) %% pad_c) + 1
            k_pad[row_idx, col_idx] <- kernel
            kernel_fft_cache[[key]] <<- fft(k_pad)
        }
        kernel_fft_cache[[key]]
    }

    # Convolves a list of real matrices with the kernel (two matrices per FFT)
    convolve_layers <- function(layers, pad_r) {
        k_fft <- get_kernel_fft(pad_r)
        n_r <- nrow(layers[[1]])
        n_c <- ncol(layers[[1]])
        convolved <- vector("list", length(layers))

        for (s in seq(1, length(layers), by = 2)) {
            has_pair <- (s + 1) <= length(layers)
            z <- matrix(0 + 0i, pad_r, pad_c)
            z[seq_len(n_r), seq_len(n_c)] <- complex(
                real = layers[[s]],
                imaginary = if (has_pair) layers[[s + 1]] else 0)
            res <- fft(fft(z) * k_fft, inverse = TRUE) / (pad_r * pad_c)
            
            convolved[[s]] <- Re(res)[seq_len(n_r), seq_len(n_c), drop = FALSE]
            if (has_pair) {
                convolved[[s + 1]] <- Im(res)[seq_len(n_r), seq_len(n_c), drop = FALSE]
            }
        }
        convolved
    }

    # Output: one layer of proportion per category, written block by block
    prop_file <- file.path(tmp_dir, "prop_stack.tif")
    prop_stack <- rast(r_proj, nlyrs = n)
    writeStart(prop_stack, prop_file, overwrite = TRUE, datatype = "FLT4S")
    
    for (b in seq_along(block_starts)) {
        cli_alert_info("Processing proportions, block {b} of {length(block_starts)}")
        out_first <- block_starts[b]
        out_last <- min(out_first + block_rows - 1, n_row)
        in_first <- max(1, out_first - half_r)
        in_last <- min(n_row, out_last + half_r)

        values_in <- read_codes(in_first, in_last - in_first + 1)
        is_valid <- !is.na(values_in)

        # layer 1: valid (non-NA) pixels, layers 2..n+1: indicator of each class
        layers <- c(
            list(is_valid * 1),
            lapply(present_codes, function(code) (is_valid & (values_in == code)) * 1))
        pad_r <- nextn(nrow(values_in) + half_r, factors = c(2, 3, 5))
        convolved <- convolve_layers(layers, pad_r)

        # keep only the rows of the block (without the margins)
        keep <- (out_first - in_first + 1):(out_last - in_first + 1)
        n_valid <- convolved[[1]][keep, , drop = FALSE]
        no_data <- n_valid < 0.5   # no valid pixel in the circle -> NA

        block_values <- matrix(NA_real_, nrow = length(keep) * n_col, ncol = n)
        for (i in seq_len(n)) {
            prop <- convolved[[i + 1]][keep, , drop = FALSE] / n_valid
            prop[no_data] <- NA
            prop <- pmin(pmax(prop, 0), 1)   # removes floating point noise
            block_values[, i] <- as.vector(t(prop))   # terra wants cells by row
        }
        writeValues(prop_stack, block_values, out_first, length(keep))
        rm(values_in, is_valid, layers, convolved, n_valid, block_values); gc()
    }
    prop_stack <- writeStop(prop_stack)
    readStop(r_proj)

    # Sanity check: at pixels with data, proportions of all categories sum to 1
    sum_check <- global(app(prop_stack, sum), "range", na.rm = TRUE)
    if (any(abs(unlist(sum_check) - 1) > 1e-3)) {
        stop("Proportions do not sum to 1 (range: ", 
             paste(round(unlist(sum_check), 4), collapse = " - "), 
             "). Codes found in pixels: ", paste(present_codes, collapse = ", "))
    }
    names(prop_stack) <- paste0(
        "prop_", present, "_within_", as.integer(radius_m), "m")

    # Back to the original grid
    out <- project(prop_stack, raster, method = "bilinear",
                    filename = file.path(tmp_dir, "prop_final.tif"),
                    datatype = "FLT4S", overwrite = TRUE)
    names(out) <- names(prop_stack)
    out
}

# A function to automatically simplify CLC categories and reduce its extent.
# ARGS:
#   - raster_path: a filepaths to a raster (.tif) file.
#   - save_to: a filepath to save the new raster (must end with .tif).
#   - lon_min, lon_max, lat_min, lat_max: extent for the output (in PROJECT_CRS).
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
        cropped_raster <- crop(
            raster, 
            ext(get_extent_polygon(
                crs(raster), 
                lon_min, lon_max, lat_min, lat_max
            ))
        )

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
        
        # replace placeholder label
        names(simplified_clc_raster) <- "land.cover"
        
        cli_alert_info(paste0(
            "New categories: ",
            paste(
                unique(simplified_clc_raster$land.cover)$land.cover, 
                collapse = ", "
            ), 
            "."
        ))

        # Once the categories are set, we can compute distances to each
        # and proportion of each around pixels, this gives contextual 
        # information that the model cannot have alternatively.
        cli_alert_info("Adding one layer per categories (distance layers)...")
        distance_layers <- factorial_to_distances(simplified_clc_raster)

        cli_alert_info("Adding one layer per categories (proportion layers)...")
        # Using square diagonal, not its length for radius
        rad <- sqrt(2*((1000*RES_KM)**2))
        proportion_layers <- factorial_to_proportion(
            simplified_clc_raster, radius_m = rad)
        # Exclude proportion of "0" -> artifact coming from NA values.
        proportion_layers[[paste0("prop_0_within_", as.integer(rad),"m")]] <- NULL

        cli_alert_info("Stacking new layers together...")
        clc_stacked_raster <- rast(
            list(simplified_clc_raster, distance_layers, proportion_layers))
        
        cli_alert_info("Saving file...")
        writeRaster(clc_stacked_raster, save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created simplified CORINE land cover raster file!\n\n"))
        
        return(fix_names(clc_stacked_raster))
    }
}

# A function to clip a raster to france metropolitan area
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
clip_raster_france <- function(
        raster,  
        buffer, 
        verbose = TRUE,
        save_to = NULL) {
    
    # import france shapefile with buffer to avoid clipping important data
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>    
        st_buffer(dist = buffer*1000) |> # buffer distance in meters
        st_transform(crs(raster))        # transform to match raster's crs

    # reduce size of raster
    raster_cropped <- crop(raster, ext(france_sf_buffered))
    raster_clipped <- clip_raster_from_shapefile(
        raster_cropped, france_sf_buffered)

    if (verbose) {
        cli_alert_info("Saving file...\n")
    }

    if (is.null(save_to)) {
        return(raster_clipped)
    } else {
        writeRaster(raster_clipped, save_to, overwrite = TRUE)  
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
        "_projection_france_res", res_km, "km-", PROJECT_CRS_TAG, ".tif")
        
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
    
    cli_alert_info(paste0("Pre-processing of raster..."))

    # Creating the paths that will be needed in the function
    shapefile_hexagons_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_projection_france_hexagons_res", res_km, "km-", PROJECT_CRS_TAG, ".gpkg")


    if (!authorise_overwrite(shapefile_hexagons_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        # import and reproject to a standardized grid of hexagons
        if (!inherits(raster_path, "SpatRaster")) {
            raster <- rast(raster_path)
        } else {
            raster <- raster_path
        }
        
        . <- project_to_hexagons(
                raster = raster, 
                save_to = shapefile_hexagons_path, 
                res_km = res_km,
                numeric_stat = extraction_mode) 
        cli_alert_success("Saved shapefile and CSV legend!\n\n")
    }

    return(shapefile_hexagons_path)
}
