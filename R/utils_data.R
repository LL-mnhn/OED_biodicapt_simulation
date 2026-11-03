# Set of utility functions used accross scripts.

##### Liraries ##### ----------------------------------------------------------
library(sf)
library(cli)
library(FNN)
library(readr)
library(dplyr)
library(terra)
library(tools)
library(stringr)
library(dggridR)
library(rnaturalearth)
library(exactextractr)


##### Parameters ##### --------------------------------------------------------
source(here::here("data/config/config.R")) # all parameters are grouped together
BUFFER_SIZE = 5*RES_KM*1000  # Common buffer size for all functions (in metres)


##### Global functions ##### --------------------------------------------------
# A function to ask for a inputs by a user, which works in both 
# interactive (with user inputs) and batch mode (on remote cluster)
# ARGS:
#   - prompt: a string. The message to display before asking for input.
#   - remote_username : a string. The name of the user for which there should 
#       not be an interactive session (check with Sys.info()["user"] in remote 
#       terminal).
#   - default_remote_behaviour: a string. The default answer when remote.
# RETURNS:
#   - a string with the user answer or auto answer when on a remote machine.
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
        stop("Unknown environment.")
    }
    return(txt)
}

# A function to ask the user if he wants to overwrite a file.
# ARGS:
#   - path: a path (folder or file) to check.
#   - remote_username : the name of the user for which there should not be an 
#       interactive session (check with Sys.info()["user"] in remote terminal).
#   - default_remote_behaviour: a string. The default answer when remote.
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

# A function that checks if two objects have the same CRS
# ARGS:
#   - x1: a terra object or a string (ex: "EPSG:4326")
#   - x2: a terra object or a string (ex: "EPSG:4326")
# RETURNS:
#   A boolean: TRUE if the two CRS are equivalent.
check_crs <- function(x1, x2) {
    same.crs(crs(x1, proj = TRUE), crs(x2, proj = TRUE))
}



##### Datasets ##### ----------------------------------------------------------
# A helper that makes sure a categorical raster is in PROJECT_CRS and written 
# on disk as integers (CORINE codes are small ints). If the raster is already in
# PROJECT_CRS it is only written to disk (no resampling, so no distortion).
# ARGS:
#   - raster: a SpatRaster (single categorical layer).
# RETURNS:
#   - a SpatRaster (file-backed) in PROJECT_CRS.
to_project_crs_on_disk <- function(raster) {
    filename <- file.path(tempdir(), "r_proj.tif")
    if (check_crs(raster, PROJECT_CRS)) {
        writeRaster(
            fix_names(raster), 
            filename, 
            datatype = "INT2U", overwrite = TRUE)
    } else {
        project(
            fix_names(raster), 
            PROJECT_CRS, 
            method = "near",
            filename = filename, 
            datatype = "INT2U", 
            overwrite = TRUE)
    }
}

# A function that "blurs" coordinate within a dataframe by randomly shifting
# each longitude and latitude coordinates.
# ARGS:
#   - df: a dataframe, sf object or SpatVector with gps coordinates in EPSG:4326
#   - x_lon: a string. The name of the column with longitude values.
#   - y_lat: a string. The name of the colum with latitude values.
#   - blur_km: an integer/a float. 
#       The range within which the blurring occurs (in km).
#   - seed: an integer. When given, sets a seed at function level.
# RETURNS:
#   - orginal df with shifted coordinates.
blur_coordinates <- function(df, x_lon, y_lat, blur_km, seed = NULL) {
    fun_crs <- "EPSG:4326"
    
    # for reproducible results  
    if (!is.null(seed)){
        set.seed(seed) 
    }

    # if sf object, remove geometry field (we're gonna make a new one)
    is_sf <- inherits(df, "sf")
    is_sv <- inherits(df, "SpatVector")
    if (is_sf) {
        if (!check_crs(df, fun_crs)) { 
            stop(paste("CRS is not ", fun_crs))
        }
        crs_orig <- sf::st_crs(df)
        geom_col <- attr(df, "sf_column")
        df <- sf::st_drop_geometry(df)  # attributes only, geometry dropped
    } else if (is_sv) {
        if (!check_crs(df, fun_crs)) { 
            stop(paste("CRS is not ", fun_crs))
        }
        crs_orig <- terra::crs(df)
        df <- as.data.frame(df)   # attributes only, geometry dropped
    } else if (!(
        all(FRANCE_BBOX_4326["xmin"] <= as.data.frame(df[x_lon])) &
        all(as.data.frame(df[x_lon]) <= FRANCE_BBOX_4326["xmax"]) &
        all(FRANCE_BBOX_4326["ymin"] <= as.data.frame(df[y_lat])) &
        all(as.data.frame(df[y_lat]) <= FRANCE_BBOX_4326["ymax"])    
    )) {
        # Check that the dataframe is in EPSG:4326 coordinates (quick security)
        stop(paste0(
            "blur_coordinates() assumes that your data is in",
            "EPSG:4326 over France's extent. It looks like 'df' is in ",
            "another CRS, or on a different area. Please check your data."
        ))
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
    return(fix_names(df))
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

# A function that clips the values of a raster outside of an overlapping 
# shapefile
# ARGS:
#   - raster: a SpatRaster.
#   - shapefile: a sf object.
# RETURNS:
#   - a SpatRaster clipped to sf borders.
clip_raster_from_shapefile <- function(raster, shapefile) {
    # verify raster and shapefile crs are marching
    if (!check_crs(raster, shapefile)) {
        shapefile <- st_transform(shapefile, crs(raster))  # fix silently
    }
    # Clip the raster using the shapefile
    fix_names(mask(raster, shapefile, touches = TRUE))
}

# A function to propagate the nearest valid values of a raster into NA cells.
# ARGS:
#   - raster: a SpatRaster.
#   - save_to: a string or NULL. 
#       If NULL, returns the raster after processing. 
#       If not, the string given is the path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - max_distance: an integer/a float. 
#       Is used to define a buffer around the raster (in meters).
# RETURNS:
#   a SpatRaster.
propagate_raster_values <- function(
        raster,
        max_distance = BUFFER_SIZE * 1.05,   # in meters
        save_to = NULL,
        progress = TRUE,
        k = 49,                # number of frontier cells used for IDW (numeric layers)
        fade_power = 0.33,     # 1 = linear fade, >1 = stays near border values longer
        regional_value = NULL  # NULL = mean of valid cells per layer; or a single number /
                               # a vector with one value per layer
    ) {

    # Security checks
    if (!is_metric(raster) && !is.lonlat(raster)) {
        stop("propagate_raster_values() needs a metric or lon/lat CRS.")
    }

    is_cat <- is.factor(raster)   # one logical per layer
    if (sum(global(is.na(raster), "sum")[[1]]) == 0) {
        if (is.null(save_to)) {
            return(raster)
        } else {
            dtype <- if (all(is_cat)) "INT1U" else "FLT4S"
            writeRaster(fix_names(raster), save_to, overwrite = TRUE, datatype = dtype)
            invisible(raster)
        }
    }

    # Silence terra's own per-operation bars; ours replaces them
    old_opt <- terraOptions(progress = 0, print = FALSE)
    on.exit(terraOptions(progress = old_opt$progress), add = TRUE)

    # Class codes per categorical layer
    codes_list <- lapply(seq_len(nlyr(raster)), function(i) {
        if (!is_cat[i]) return(NULL)
        lyr <- raster[[i]]
        levels(lyr) <- NULL
        coltab(lyr) <- NULL
        codes <- unique(lyr)[[1]]
        sort(codes[!is.na(codes)])
    })

    total_steps <- sum(vapply(seq_len(nlyr(raster)), function(i) {
        if (is_cat[i]) length(codes_list[[i]]) else 1
    }, numeric(1)))

    # Progress bar
    step <- 0
    pb <- if (progress && total_steps > 0) {
        txtProgressBar(min = 0, max = total_steps, style = 3)
    }
    tick <- function() {
        step <<- step + 1
        if (!is.null(pb)) setTxtProgressBar(pb, step)
    }
    on.exit(if (!is.null(pb)) close(pb), add = TRUE)

    # Categorical: nearest-class fill, one distance() call per class
    fill_nearest_cat <- function(lyr, codes, max_distance) {
        best_d <- init(lyr, Inf)
        best_c <- init(lyr, NA)

        for (code in codes) {
            d <- distance(ifel(lyr == code, 1, NA))
            closer <- d < best_d
            best_d <- ifel(closer, d, best_d)
            best_c <- ifel(closer, code, best_c)
            tick()
        }

        ifel(is.na(lyr) & best_d <= max_distance, best_c, lyr)
    }

    # Numeric: IDW from the k nearest frontier cells, faded toward a regional
    # value as distance from valid data approaches max_distance
    fill_nearest_numeric <- function(lyr, max_distance, mu) {
        is_na <- is.na(lyr)

        # Distance (m) from each cell to the nearest valid cell
        d      <- distance(ifel(is_na, NA, 1))
        d_vals <- values(d, mat = FALSE)
        target <- which(values(is_na, mat = FALSE) & d_vals <= max_distance)

        # Frontier = valid cells touching an NA cell
        touches_na <- focal(is_na, w = 3, fun = "max", fillvalue = 0)
        src        <- which(values(!is_na & touches_na == 1, mat = FALSE))

        if (length(target) > 0 && length(src) > 0) {
            kk <- min(k, length(src))
            nn <- FNN::get.knnx(
                xyFromCell(lyr, src), xyFromCell(lyr, target), k = kk)

            # IDW of the nearest frontier values
            eps  <- res(lyr)[1] * 1e-3
            w    <- 1 / (nn$nn.dist + eps)^2
            vals <- matrix(values(lyr, mat = FALSE)[src[nn$nn.index]], ncol = kk)
            idw  <- rowSums(vals * w) / rowSums(w)

            # Fade toward regional value: 0 at the border, 1 at max_distance
            f <- pmin(d_vals[target] / max_distance, 1)^fade_power
            lyr[target] <- (1 - f) * idw + f * mu
        }

        tick()
        lyr
    }

    # Regional value per numeric layer
    if (is.null(regional_value)) {
        mu_all <- global(raster, "mean", na.rm = TRUE)[[1]]
    } else {
        mu_all <- rep_len(regional_value, nlyr(raster))
    }

    layers <- lapply(seq_len(nlyr(raster)), function(i) {
        lyr <- raster[[i]]
        nm  <- names(lyr)

        if (is_cat[i]) {
            cat_table <- cats(lyr)[[1]]
            active    <- activeCat(lyr)
            ct        <- coltab(lyr)[[1]]

            levels(lyr) <- NULL
            coltab(lyr) <- NULL

            lyr <- fill_nearest_cat(lyr, codes_list[[i]], max_distance)

            lyr <- as.int(lyr)
            levels(lyr) <- cat_table
            activeCat(lyr) <- active
            if (!is.null(ct)) coltab(lyr) <- ct
        } else {
            lyr <- fill_nearest_numeric(lyr, max_distance, mu_all[i])
        }
        names(lyr) <- nm
        lyr
    })
    new_raster <- rast(layers)

    if (is.null(save_to)) {
        return(new_raster)
    } else {
        dtype <- if (all(is_cat)) "INT1U" else "FLT4S"
        writeRaster(fix_names(new_raster), save_to, overwrite = TRUE, datatype = dtype)
        invisible(new_raster)
    }
}

# A function to clip a raster to france metropolitan area
# ARGS:
#   - raster: a SpatRaster.
#   - save_to: a string or NULL. 
#       If NULL, returns the raster after processing. 
#       If not, the string given is the path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - buffer: an integer/a float. 
#       Is used to define a buffer around the raster (in meters).
#   - fill_missing: a boolean. If TRUE, fills NAs inside buffer by means.
# RETURNS:
#   - clipped SpatRaster
clip_raster_to_france_shape <- function(
        raster,  
        buffer = BUFFER_SIZE, 
        fill_missing = TRUE,
        save_to = NULL) {
    
    # import france shapefile with buffer to avoid clipping important data
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>    
        st_buffer(dist = buffer) |> # buffer distance in meters
        st_transform(crs(raster))   # transform to match raster's crs

    # reduce size of raster (MAKE IT FASTER)
    raster_cropped <- crop(raster, ext(france_sf_buffered), extend = TRUE)

    # Fill missing values
    if (fill_missing) {
        raster_cropped <- propagate_raster_values(
            raster_cropped, max_distance = buffer*1.05)
    }

    # Clip to buffer size
    raster_clipped <- clip_raster_from_shapefile(
        raster_cropped, france_sf_buffered)
    raster_clipped <- fix_names(raster_clipped)

    if (is.null(save_to)) {
        return(raster_clipped)
    } else {
        writeRaster(raster_clipped, save_to, overwrite = TRUE)  
        return(invisible(raster_clipped))
    }   
}

# A function to aggregate rasters together (months -> year)
# ARGS:
#   - raster_paths: filepaths to 12 (monthly) rasters of the same data.
#   - fun: a function. Tells how the data should be aggregated together.
# RETURNS:
#   - A single spatRaster with the aggregated data.
monthly_2_yearly_rasters <- function(
        raster_paths, 
        fun = mean){
    
    # Verify we have 12 paths (for 12 months)
    if (length(raster_paths) != 12){
        stop(paste("List of strings recieved contains", length(raster_paths), "elements, expected 12."))
    }

    # Stack rasters together
    for (n_path in 1:12){
        # import raster
        raster <- rast(raster_paths[n_path])
        
        # Standardise size 
        # (larger buffer because raster_paths are in EPSG:4326)
        cropped_raster <- clip_raster_to_france_shape(
            raster, buffer = BUFFER_SIZE*1.25, fill_missing = FALSE)   

        if (n_path == 1){
            rasters_together <- cropped_raster
        } else {
            rasters_together <- c(rasters_together, cropped_raster)
        }
    }
    
    # Compute and return average/median/other on list of cropped raster
    annual_raster <- app(rasters_together, fun = fun)

    # project to PROJECT_CRS
    annual_raster <- project(annual_raster, PROJECT_CRS)
    annual_raster <- clip_raster_to_france_shape(
            annual_raster, fill_missing = TRUE) 

    return(fix_names(annual_raster))
}

# A function to standardize monthly raster aggregation and save path.
# ARGS:
#   - folder_path: path to a folder containing 12 raster files.
#   - save_to: a filepath to save the new raster (must end with .tif).
# RETURNS:
#   - A single spatRaster with the aggregated data.
average_monthly_quantitative_rasters <- function(folder_path, save_to) {
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
            raster_paths = twelve_paths, fun = mean)
        # monthly_2_yearly_rasters already returns fix_names(raster)
        
        cli_alert_info("Saving file...")
        writeRaster(annual_raster, save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created '", save_to, "' file!\n\n"))
        return(invisible(annual_raster))
    }
}


# A function that stacks together all preprocessed rasters/sf (hard coded paths)
# into one single raster/sf.
# ARGS:
#   - save_to: the path where the output file will be saved 
#       (use ".gpkg" if mode is "hexagonal", ".tif" for "raw" and "raster").
#   - mode: a string. Three modes available:
#       - "raw": builds a raster with the best possible resolution (uses the 
#           coarsest resolution of all preprocessed rasters).
#       - "raster": builds a raster with from standardised resolutions
# RETURNS:
#   - a sf object, the dataframe with its added features (= columns)
create_master_files_features <- function(
        save_to, mode = "raw") {

    # TRUE if `r` can be reduced to `template` by an integer-factor aggregation:
    # same CRS, template resolution is an integer multiple (>= 2) of the source
    # resolution, and the source grid lines up with the template cell boundaries.
    can_aggregate <- function(r, template, tol = 1e-6) {
        if (!check_crs(r, template)) return(FALSE)
 
        fact <- res(template) / res(r)
        if (any(fact < 2 - tol) || any(abs(fact - round(fact)) > tol)) return(FALSE)
 
        rt <- res(template)
        dx <- (xmin(r) - xmin(template)) / rt[1]
        dy <- (ymax(r) - ymax(template)) / rt[2]
        abs(dx - round(dx)) <= tol && abs(dy - round(dy)) <= tol
    }
 
    # Resample `r` onto `template`. `methods` has one entry per layer of `r`:
    # "near" (categorical) or "bilinear" (continuous).
    #   - same grid already          -> returned as is
    #   - aligned, integer factor    -> aggregate() (modal / mean)
    #   - anything else              -> project()
    # Layers are processed per method group, not one by one.
    resample_to_template <- function(r, template, methods) {
        is_near <- methods == "near"
        parts <- list()
 
        if (can_aggregate(r, template)) {
            fact <- as.integer(round(res(template) / res(r)))
            cli_alert_info(paste0(
                "Aligned grid: using aggregate() with factor ",
                paste(fact, collapse = "x"), "."))
 
            if (any(is_near)) {
                parts$near <- terra::aggregate(
                    r[[which(is_near)]], fact = fact, fun = "modal", na.rm = TRUE)
            }
            if (any(!is_near)) {
                parts$cont <- terra::aggregate(
                    r[[which(!is_near)]], fact = fact, fun = "mean", na.rm = TRUE)
            }
        } else {
            cli_alert_info("Grid not aligned or different CRS: using project().")
 
            if (any(is_near)) {
                parts$near <- terra::project(
                    r[[which(is_near)]], template, method = "near")
            }
            if (any(!is_near)) {
                parts$cont <- terra::project(
                    r[[which(!is_near)]], template, method = "bilinear")
            }
        }
 
        out <- do.call(c, unname(parts))
        out <- out[[names(r)]]  # restore original layer order
 
        # snap to the template if aggregate() left a slightly different extent
        same_grid <- isTRUE(all.equal(as.vector(ext(out)), as.vector(ext(template)))) &&
            all(dim(out)[1:2] == dim(template)[1:2])
        if (!same_grid) {
            out <- terra::resample(out, template, method = "near")
        }
 
        # restore factor levels if they were lost along the way
        for (i in which(is.factor(r))) {
            if (!is.factor(out[[i]])) {
                levels(out[[i]]) <- cats(r)[[i]]
            }
        }
        out
    }

    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning(paste0(
            "Skipping the creation of a raster containing containing all the ",
            "features together (same CRS, same pixels).\n\n"))

    } else {
        terraOptions(memfrac = 0.5, progress = 1)
        # the user can choose which type of file the features must be extracted
        # from. If raw, uses the initial raw resolution of the files (which are
        # usually finer than 1km² in resolution). If using raster or hexagonal,
        # uses pre-processed data that was made to be pixels/hexagons space out
        # by RES_KM.

        # Get paths of datasets (outside the main "merging" if/else to be able 
        # to plot the paths used before beginning the computation).
        if (mode == "raw") {
            environmental_data_paths <- c(
                CORINE_PATH_PREPROCESSED, CHELSA_PATHS_PREPROCESSED, 
                NDVI_PATH_PREPROCESSED, ELEVATION_PATH_PREPROCESSED, 
                LIGHT_POLLUTION_PATH_PREPROCESSED, SWF_PATH_PREPROCESSED,
                PESTICIDE_PATH_PREPROCESSED)
            
        } else if (mode == "raster"){
            environmental_data_paths <- c(
                CORINE_PROJECTION_BASEPATH, CHELSA_PROJECTION_BASEPATHS, 
                NDVI_PROJECTION_BASEPATH, ELEVATION_PROJECTION_BASEPATH, 
                LIGHT_POLLUTION_PROJECTION_BASEPATH, SWF_PROJECTION_BASEPATH,
                PESTICIDE_PROJECTION_BASEPATH)  
            environmental_data_paths <- paste0(
                environmental_data_paths, 
                "_res", RES_KM, "km_", PROJECT_CRS_TAG, ".tif")
        } else {
            stop("mode should be one of 'raw', 'raster' or 'hexagonal'.")
        }

        cli_bullets(c(
            "!" = "[Fetching features]: Using hard coded paths:",
            setNames(environmental_data_paths, rep("*", length(environmental_data_paths)))
        ))

        # Create target extent for multi-layer raster
        france_sf <- get_metropolitan_france_shapefile()
        france_sf_buffered <- france_sf |>
            st_transform(PROJECT_CRS) |>    
            st_buffer(dist = BUFFER_SIZE) # buffer distance in meters
        target_ext <- ext(france_sf_buffered)

        # Security check (Right CRS + minimal extent)
        for (path in environmental_data_paths) {
            r <- rast(path)
            if (!check_crs(r, PROJECT_CRS)) {
                # a loser check, just to be sure
                stop(paste("CRS should be", PROJECT_CRS, "for file:", path))
            }
            if (!((abs(xmin(target_ext) - xmin(r)) <= RES_KM*1000) &&
                  (abs(xmax(r) - xmax(target_ext)) <= RES_KM*1000) &&
                  (abs(ymin(target_ext) - ymin(r)) <= RES_KM*1000) &&
                  (abs(ymax(r) - ymax(target_ext)) <= RES_KM*1000))) {
                stop(paste("Raster extent should be larger for file: ", path))
            }
        }

        # Get all resolutions in meters 
        # (LEGACY: degrees converted via the N-S length of a degree)
        resolutions <- vapply(environmental_data_paths, function(p) {
            if (tolower(tools::file_ext(p)) %in% c("tif", "tiff")) {
                r <- rast(p)
                if (is.lonlat(r)) res(r)[2] * 111320 else res(r)[1]
            } else {
                NA_real_
            }
        }, numeric(1))

        # Select coarsest resolution and create template
        target_res <- ceiling(max(resolutions, na.rm = TRUE) / 100) * 100  # round up to next 100 m
        template <- rast(target_ext, resolution = target_res, crs = PROJECT_CRS)
        cli_alert_info(paste0(
            "Selecting raster with coarsest resolution (res: ", 
            round(target_res, 2), " m)."))
        
        # Crop and resample each raster onto the template
        cli_alert_info("Resampling rasters, can take several minutes...")
        raster_list <- list()
        for (path in environmental_data_paths) {
            cli_alert_info(paste0("Treating file: ", path, "..."))
            if (tolower(tools::file_ext(path)) %in% c("tif", "tiff")) {
                env_raster <- rast(path)
            } else {
                stop(paste("File should be a raster:", path))
                # LEGACY for SpatVectors
                # vector <- vect(path)
                # if ("mean_tii" %in% names(vector)) {
                #     vector <- project(vector, PROJECT_CRS)
                #     env_raster <- rasterize(
                #         vector, template, 
                #         field = "mean_tii", fun = "mean", touches = TRUE)
                # } else {
                #     stop(
                #         paste0("Cannot find field name for vector: ", 
                #                 path, "."))
                # }
            }
            
            # cheap pre-crop
            if (!check_crs(env_raster, PROJECT_CRS)) {
                ext_native <- project(target_ext, from = PROJECT_CRS, to = crs(env_raster))
                env_raster <- crop(env_raster, ext_native * 1.05, extend = TRUE)
            } else {
                env_raster <- crop(env_raster, target_ext, extend = TRUE)
            }

            # Project on template grid
            # CLC: land cover layer = "near"; distance/proportion = "bilinear".
            # Everything else is continuous -> "bilinear".
            # if (grepl("CLC", path, ignore.case = TRUE)) {
            #     methods <- ifelse(
            #         grepl("^(distance|prop_)", names(env_raster), ignore.case = TRUE),
            #         "bilinear", "modal")
            # } else {
            #     methods <- rep("bilinear", nlyr(env_raster))
            # }
            # use "near": less precise but way faster to compute
            methods <- rep("near", nlyr(env_raster))
            raster_resampled <- resample_to_template(env_raster, template, methods)

            # add buffer + fill missing values
            raster_resampled <- clip_raster_to_france_shape(raster_resampled, fill_missing = FALSE)
            name <- strsplit(basename(path), "_")[[1]][2]

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
        cli_alert_info("Stacking layers together...")
        names(raster_list) <- NULL # reset names to avoid auto-renaming
        master_raster <- fix_names(rast(raster_list))

        # Save
        cli_alert_success("Saving raster...")
        writeRaster(master_raster, save_to, overwrite = TRUE)
        cli_alert_success("Saved stacked raster!\n\n")
        return(invisible(master_raster))
    }
}

# A function that extracts features from all preprocessed rasters (hard coded 
# paths) and matches them to the coordinates in the given (preprocessed as well)
# dataframe.
# ARGS:
#   - file: a string or a data.frame.
#       Must contain columns "LON" and "LAT" (crs = 4326)
#       If a string, reads it as a .csv file. If a data.frame, uses it directly. 
#   - raster_path: a string. The path to a multi-layer SpatRaster.
#       WARNING: it must have as many columns are in the hard coded
#       'environmental_data_names' variable (and order matters btw, glhf).
#   - save_to: the path where the output sf will be saved (use ".gpkg")
#   - blur: a boolean. If TRUE, blurs the coordinates of the original dataframe
#       before saving them/returning them (but after the extraction of the
#       features to be as precise as possible).
# RETURNS:
#   - a SpatVector object, the dataframe with its added features (= columns)
save_features_from_obs <- function(file, raster_path, save_to, blur = TRUE) {
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
                "save_features_from_obs() assumes that the given `file` is in",
                "EPSG:4326 over France's extent. It looks like 'df' is in ",
                "another CRS, or on a different area. Please check your data."
            ))
        }

        # convert dataframe to geo dataframe
        cli_alert_info(paste0("Conversion to shapefile..."))
        points <- vect(st_as_sf(
            df, coords = c(x_lon, y_lat), crs = 4326, remove = FALSE))

        # Step 1: extract the values of each layer at each point
        env_data <- rast(raster_path)
        pts_i <- project(points, crs(env_data))
        # drop the leading ID column, keep the other columns
        extracted_values <- terra::extract(
            env_data, pts_i)[, -1, drop = FALSE]

        # Security check
        if (anyNA(as.data.frame(extracted_values))) {
            stop("Found NAs in extracted values table, aborting operation.")
        }

        # Step 2: rename columns
        layer_ids <- names(env_data)
        colnames(extracted_values) <- layer_ids
        for (col in layer_ids) {
            points[[col]] <- extracted_values[[col]]
        }

        # Step 3: add blur of locations coordinates
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

        # Step 4: save resulting dataframe (reproject in target CRS)
        cli_alert_info("Saving file...")
        writeVector(
            fix_names(project(anonym_df, PROJECT_CRS)), 
            save_to, overwrite = TRUE)
        cli_alert_success("Saved anonymized dataset with features!")

        return(invisible(fix_names(project(anonym_df, PROJECT_CRS))))
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
    
    # Make sure names are correctly formatted
    master_raster <- fix_names(master_raster)
    
    # Keep only the necessary layers
    if (!is.null(keep.lyr)) {
        master_raster <- master_raster[[keep.lyr]]
    }
    
    # Make sure that values are limited to france's extent
    master_raster <- clip_raster_to_france_shape(master_raster)
 
    return(master_raster)
}


##### Rasters/SF functions ##### ----------------------------------------------
# A function that creates an empty raster that matches France's extent
# ARGS:
#   - resolution: an integer/a float. 
#       The size for each pixel of the raster (in km).
# RETURNS:
#   - a SpatRaster with empty values to the extent of France's metropolitan
#       area in PROJECT_CRS
get_france_raster_template <- function(resolution) {
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
        st_buffer(dist = BUFFER_SIZE)      # buffer distance in meters
    rast(
        extent = ext(france_sf_buffered),
        resolution = resolution * 1000, crs = PROJECT_CRS
    )
}

# A function to transform a raster into a (raster) grid layered on france
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - res_km: an integer/a float. 
#       Is used to define a buffer around the raster (in km).
# RETURNS:
#   - NULL, saves the re-projected SpatRaster silently.
project_raster_to_france_grid <- function(raster, save_to, res_km) {   
    # Validate resampling method for each layer
    method <- c()
    for (i in seq_len(nlyr(raster))) {
        is_cat <- is.factor(raster[[i]])

        if (is_cat) {
            # allowed <- c("near", "mode")
            method <- c(method, "near")
        } else {
            # allowed <- c("near", "bilinear", "cubic", "cubicspline",
            #   "lanczos", "average", "sum", "min", "max",
            #   "median", "mode", "rms")
            method <- c(method, "bilinear")
        }
    }

    # fit raster to new grid
    cli_alert_info("Fitting new grid...\n")
    france_grid_template <- get_france_raster_template(res_km)
    projected_layers <- lapply(
        seq_len(nlyr(raster)), function(i) {
            cli_alert_info("Projecting layer {i}/{nlyr(raster)}: {names(raster)[i]}")
            project(
                raster[[i]], 
                france_grid_template, 
                method = method[i],
                threads = TRUE)
        }
    )
    raster_new_grid <- do.call(c, projected_layers)

    cli_alert_info("Masking...\n")
    # remove values outside of france borders
    raster_new_grid_clean <- clip_raster_to_france_shape(raster_new_grid)

    cli_alert_info("Saving file...\n")
    writeRaster(fix_names(raster_new_grid_clean), save_to, overwrite = TRUE)  
    
    return(invisible(fix_names(raster_new_grid_clean)))
}

# A function to rasterise a shapefile into a (raster) grid layered on france
# ARGS:
#   - vector: a shapefile (spatVector) object.
#   - layer_name: a string.
#   - save_to: a string. The path where the raster will be saved. 
#       Filepath must end in ".tif" to save raster.
#   - touches: a boolean, passed directly to terra::rasterize()
#   - res_km: an integer/a float. Grid resolution in km.
# RETURNS:
#   - NULL, saves the re-projected SpatRaster silently.
rasterise_to_france_grid <- function(
        vector, layer_name, save_to, res_km, touches = TRUE) {   
    
    # Security checks
    if (ncol(vector) != 1) {
        stop(paste0(
            "Expected a shapefile with 1 column, got ", ncol(vector), "."))
    }
    
    # Rasterise
    # fit raster to new grid
    cli_alert_info("Rasterising to a template grid...\n")
    france_grid_template <- get_france_raster_template(res_km)
    raster_new_grid <- rasterize(
        project(vector, crs(france_grid_template)), 
        france_grid_template, 
        field = layer_name, fun = "mean", touches = touches)
    
    cli_alert_info("Masking...\n")
    raster_new_grid_clean <- clip_raster_to_france_shape(raster_new_grid)
    
    cli_alert_info("Saving file...\n")
    if (is.null(save_to)) {
        return(raster_new_grid_clean)
    } else {
        writeRaster(raster_new_grid_clean, save_to, overwrite = TRUE)   
        return(invisible(raster_new_grid_clean))
    }

}

# A function to transform a raster into a (shapefile) grid of hexagons layered 
# on france metropolitan area.
# ARGS:
#   - raster: a raster.
#   - save_to: a string. The filepath where the shapefile will be saved. 
#       Filepath must end in ".gpkg" or ".shp" to save sf.
#   - res_km: an integer/a float. Grid resolution in km.
#   - numeric_stat: a string. 
#       If raster contains numerical values, either "mean" or "median".
#   - lat_min, lat_max: latitudinal extent for the output (in PROJECT_CRS).
# RETURNS:
#   - NULL, saves the shapefile silently.
project_raster_to_hexagons <- function(
        raster, save_to, res_km, numeric_stat = "mean") {
    # Security checks
    if (!numeric_stat %in% c("mean", "median")) {
        stop(paste0("Unknown numeric_stat: '", numeric_stat, 
            "'. Use 'mean' or 'median'."))
    }
    if (!inherits(raster, "SpatRaster")) {
        stop(paste0("raster must be of class SpatRaster, got ", class(raster)))
    }
    
    # Determine, per layer, whether it's categorical
    is_cat <- is.factor(raster)
    layer_names <- names(raster)
    n_layers <- terra::nlyr(raster)

    # 1. Define grid of hexagons
    cli_alert_info("Creating grid of hexagons...\n")
    dggs <- dgconstruct(
        spacing = res_km, metric = TRUE, resround = 'nearest')
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
        st_buffer(dist = BUFFER_SIZE) |>   # buffer distance in meters
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
    cli_alert_info(paste0(
        "Extracting values from ", n_layers, 
        " layer(s), this might take some time...\n"))

    # Polygons must be in the raster's CRS for the extraction
    hex_for_extract <- st_transform(hex_grid, crs(raster))

    extracted <- exact_extract(
        raster, hex_for_extract,
        fun = function(df) {
            w <- df$coverage_fraction
            vapply(seq_len(n_layers), function(i) {
                v <- df[[i]]
                ok <- !is.na(v) & w > 0
                if (!any(ok)) return(NA_real_)
                v <- v[ok]; wi <- w[ok]
                if (is_cat[i]) {
                    tot <- tapply(wi, v, sum)
                    as.numeric(names(tot)[which.max(tot)])
                } else if (numeric_stat == "mean") {
                    weighted.mean(v, wi)
                } else {
                    # weighted median, e.g. via Hmisc/matrixStats::weightedMedian
                    matrixStats::weightedMedian(v, wi)
                }
            }, numeric(1))
        },
        summarize_df = TRUE,
        progress = FALSE
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
    cli_alert_info("Saving files...\n")
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
    return(invisible(list(hex_grid = hex_grid, legend = legend)))
}

# A function to transform a shapefile into a (shapefile) grid of hexagons 
# layered on france's metropolitan area.
# ARGS:
#   - vector: a SpatVector (shapefile) object.
#   - fields: a string or NULL (all columns).
#   - save_to: a string. The filepath where the shapefile will be saved. 
#       Filepath must end in ".gpkg" or ".shp" to save sf.
#   - res_km: an integer/a float. Grid resoution in km.
#   - numeric_stat: a string. 
#       If raster contains numerical values, either "mean" or "median".
#   - lat_min, lat_max: latitudinal extent for the output (in PROJECT_CRS).
# RETURNS:
#   - NULL, saves the shapefile silently.
project_shp_to_hexagons <- function(
        vector, save_to, res_km, numeric_stat = "mean", fields = NULL) {

    # Security checks
    if (!numeric_stat %in% c("mean", "median")) {
        stop(paste0("Unknown numeric_stat: '", numeric_stat,
            "'. Use 'mean' or 'median'."))
    }
    if (!inherits(vector, "SpatVector")) {
        stop("`vector` must be a terra SpatVector.")
    }

    # Convert to sf, fix invalid geometries, and move to the equal-area CRS
    # (so that areas and lengths used as weights are meaningful)
    vec_sf <- sf::st_as_sf(vector) |>
        sf::st_make_valid() |>
        sf::st_transform(PROJECT_CRS)

    # Fields to aggregate (default: all attribute columns)
    geom_col <- attr(vec_sf, "sf_column")
    if (is.null(fields)) {
        fields <- setdiff(names(vec_sf), geom_col)
    }
    missing_fields <- setdiff(fields, names(vec_sf))
    if (length(missing_fields) > 0) {
        stop("Unknown field(s): ", paste(missing_fields, collapse = ", "))
    }
    if (length(fields) == 0) stop("No attribute fields to aggregate.")

    # Determine, per field, whether it's categorical
    is_cat <- vapply(fields, function(f) {
        x <- vec_sf[[f]]
        is.character(x) || is.factor(x) || is.logical(x)
    }, logical(1))

    vec_sf <- vec_sf[, fields]
    for (f in fields[is_cat]) vec_sf[[f]] <- as.character(vec_sf[[f]])
    sf::st_agr(vec_sf) <- "constant"

    # 1. Define grid of hexagons
    cli_alert_info("Creating grid of hexagons...\n")
    dggs <- dgconstruct(
        spacing = res_km, metric = TRUE, resround = 'nearest')
    france_sf <- get_metropolitan_france_shapefile()
    france_sf_buffered <- france_sf |>
        st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
        st_buffer(dist = BUFFER_SIZE) |>   # buffer distance in meters
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

    # 2. Associate each hexagon to the aggregated values of the vector
    # features that fall within it
    cli_alert_info(paste0(
        "Intersecting ", length(fields),
        " field(s) with hexagons, this might take some time...\n"))

    hex_proj <- st_transform(hex_grid, PROJECT_CRS)[, "seqnum"]
    sf::st_agr(hex_proj) <- "constant"
    inter <- suppressWarnings(st_intersection(hex_proj, vec_sf))

    # Weight of each piece: area for polygons, length for lines, 1 for points
    geom_type <- unique(as.character(sf::st_geometry_type(inter)))
    if (all(grepl("POLYGON", geom_type))) {
        inter$.w <- as.numeric(st_area(inter))
    } else if (all(grepl("LINESTRING", geom_type))) {
        inter$.w <- as.numeric(sf::st_length(inter))
    } else {
        inter$.w <- 1
    }
    inter_df <- sf::st_drop_geometry(inter)
    groups <- split(seq_len(nrow(inter_df)), inter_df$seqnum)

    weighted_median <- function(x, w) {
        o <- order(x); x <- x[o]; w <- w[o]
        x[which(cumsum(w) >= sum(w) / 2)[1]]
    }

    aggregate_field <- function(f, categorical) {
        vapply(as.character(hex_grid$seqnum), function(id, ...) {
            idx <- groups[[id]]
            if (is.null(idx)) return(if (categorical) NA_character_ else NA_real_)
            v <- inter_df[[f]][idx]
            w <- inter_df$.w[idx]
            ok <- !is.na(v) & w > 0
            v <- v[ok]; w <- w[ok]
            if (length(v) == 0) return(if (categorical) NA_character_ else NA_real_)

            if (categorical) {
                # Value covering the largest total area/length/count
                totals <- tapply(w, v, sum)
                names(totals)[which.max(totals)]
            } else if (numeric_stat == "mean") {
                weighted.mean(v, w)
            } else {
                weighted_median(v, w)
            }
        }, FUN.VALUE = if (categorical) character(1) else numeric(1),
           USE.NAMES = FALSE)
    }

    for (i in seq_along(fields)) {
        hex_grid[[fields[i]]] <- aggregate_field(fields[i], is_cat[i])
    }

    # 3. Polish and save
    cli_alert_info("Saving files...\n")
    hex_grid <- st_transform(hex_grid, PROJECT_CRS)  # saved in PROJECT_CRS
    st_write(hex_grid, save_to, delete_dsn = TRUE, quiet = TRUE)

    return(invisible(hex_grid))
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
#       (default is 1: categories are collapsed to their primary category).
# RETURNS:
#   - the updated CLC SpatRaster
simplify_CLC <- function(
        clc_raster,
        level_urban = 1,
        level_crops = 1,
        level_forests = 1,
        level_wetlands = 1,
        level_water = 1){
    cli_alert_warning(paste0(
        "[Harmless warning] This function should only be ",
        "used with CORINE Land Cover rasters!"))
    cli_alert_info("Loading categories...")

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

    # Create new table
    cli_alert_info("Creating new table...")
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

    # Build back the color table
    cli_alert_info("Assigning new colors...")
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

    cli_alert_success("Re-classified CLC2018 is ready!")
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
        raster) {

    # Tell terra to be conservative with RAM and use of disk
    terraOptions(memfrac = 0.5, progress = 1)

    raster_clipped <- clip_raster_to_france_shape(
        raster, fill_missing = FALSE)

    # Write projected raster to disk, as integers (CORINE codes are small ints)
    r_proj <- to_project_crs_on_disk(raster_clipped)
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
                  filename = file.path(tempdir(), sprintf("m_%s.tif", lab)),
                  datatype = "INT1U", overwrite = TRUE)

        dist_files[i] <- file.path(tempdir(), sprintf("dist_%s.tif", lab))
        d <- distance(m, filename = dist_files[i],
                      datatype = "FLT4S", overwrite = TRUE)
        names(d) <- paste0("distance to ", lab)

        rm(m, d); gc()
        unlink(file.path(tempdir(), sprintf("m_%s.tif", lab)))
    }

    # Virtual stack: reads lazily from disk, no RAM cost
    dist_stack <- rast(dist_files)
    names(dist_stack) <- paste0("distance to ", present)

    # Project back to the original grid, writing to disk
    out <- project(dist_stack, raster, method = "bilinear",
                   filename = file.path(tempdir(), "dist_final.tif"),
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
#   - max_fft_cells: a numeric. Maximum number of cells (approximative) of the
#       complex arrays handled at once by the FFT. Lower it if you run out of 
#       RAM (one complex cell = 16 bytes, and ~5 such arrays are alive at once).
# RETURNS:
#   - a SpatRaster with one layer of proportion per category in raster_layer.
factorial_to_proportion <- function(
        raster,
        radius_m,
        max_fft_cells = 8e6) {

    # Tell terra to be conservative with RAM and use of disk
    terraOptions(memfrac = 0.5, progress = 1)

    # Clip with a buffer
    raster_clipped <- clip_raster_to_france_shape(
        raster, fill_missing = FALSE)

    # Project to a metric CRS so the radius is in meters
    r_proj <- to_project_crs_on_disk(raster_clipped)
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
    prop_file <- file.path(tempdir(), "prop_stack.tif")
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
        "prop.", present, ".within.", as.integer(radius_m), "m")

    # Back to the original grid
    out <- project(prop_stack, raster, method = "bilinear",
                    filename = file.path(tempdir(), "prop_final.tif"),
                    datatype = "FLT4S", overwrite = TRUE)
    names(out) <- names(prop_stack)
    out
}

# A function to automatically simplify CLC categories and reduce its extent.
# ARGS:
#   - raster_path: a filepaths to a raster (.tif) file.
#   - save_to: a filepath to save the new raster (must end with .tif).
# RETURNS:
#   - a SpatRaster with several layers
#       - first layer:  
#       - next layers:
save_simplified_clc <- function(raster_path, save_to) {
    
    # load and simplify CORINE raster
    cli_alert_info("Simplification of CLC2018 raster.")
    if (!authorise_overwrite(save_to)) {
        # Check if file exists
        cli_alert_warning("Skipping simplification of CLC2018.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        raster <- rast(raster_path)
        if (check_crs(rast(CORINE_FILEPATH), PROJECT_CRS)) {
            stop(paste0("Loaded raster is not in crs ", PROJECT_CRS))
        }

        # crop raster for faster computation
        france_sf <- get_metropolitan_france_shapefile()
        france_sf_buffered <- france_sf |>
            st_transform(PROJECT_CRS) |>       # equal-area CRS (distances in meters)
            st_buffer(dist = BUFFER_SIZE)      # buffer distance in meters
        cropped_raster <- crop(raster, ext(france_sf_buffered), extend = TRUE)

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
        proportion_layers[[paste0("prop.0.within.", as.integer(rad),"m")]] <- NULL

        cli_alert_info("Stacking new layers together...")
        clc_stacked_raster <- rast(
            list(simplified_clc_raster, distance_layers, proportion_layers))
        
        cli_alert_info("Saving file...")
        writeRaster(clc_stacked_raster, save_to, overwrite = TRUE)
        cli_alert_success(paste0(
            "Created simplified CORINE land cover raster file!\n\n"))
        
        return(invisible(clc_stacked_raster))
    }
}

# A function to standardise raster projection and save path.
# ARGS:
#   - raster_path: filepaths to a raster (.tif) file.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - res_km: an integer/a float. Grid resolution in km.
# RETURNS:
#   - The path of the saved raster
raster_standardised_raster_projection <- function(
        raster_path,
        save_to_basename,
        res_km = RES_KM) {
    cli_alert_info(paste0("Pre-processing ", basename(raster_path), "..."))

    # Creating the paths that will be needed in the function
    raster_reprojected_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_res", res_km, "km_", PROJECT_CRS_TAG, ".tif")
        
    if (!authorise_overwrite(raster_reprojected_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing of raster.\n\n")

    } else {
        cli_alert_info("Loading raster...")
        # import and reproject to a standardized grid
        raster <- rast(raster_path)

        . <- project_raster_to_france_grid(
                raster = raster, 
                save_to = raster_reprojected_path,  
                res_km = res_km) 
        cli_alert_success("Modified raster file saved!\n\n")
    }

    return(raster_reprojected_path)
}

# A function to standardise shp projection, rasteresing and save.
# ARGS:
#   - shp_path: filepaths to a shapefile (.gpkg) file.
#   - layer_name: a string.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - res_km: an integer/a float. Grid resolution in km.
# RETURNS:
#   - The path of the saved raster
shp_standardised_raster_projection <- function(
        shp_path,
        layer_name,
        save_to_basename,
        res_km = RES_KM) {
    cli_alert_info(paste0("Pre-processing ", basename(shp_path), "..."))

    # Creating the paths that will be needed in the function
    raster_reprojected_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_res", res_km, "km_", PROJECT_CRS_TAG, ".tif")
        
    if (!authorise_overwrite(raster_reprojected_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing of raster.\n\n")

    } else {
        cli_alert_info("Loading shapefile...")
        # import and reproject to a standardized grid
        shp <- vect(shp_path)[, layer_name]

        . <- rasterise_to_france_grid(
                vector = shp, 
                layer_name = layer_name,
                save_to = raster_reprojected_path,  
                res_km = res_km) 
        cli_alert_success("Modified raster file saved!\n\n")
    }

    return(raster_reprojected_path)
}

# A function to standardise hexagonal projection and save path.
# ARGS:
#   - raster_path: filepaths to a raster (.tif) file.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - res_km: an integer/a float. Grid resolution in km.
# RETURNS:
#   - The reprojected SpatVector
raster_standardised_hexagonal_projection <- function(
        raster_path,
        save_to_basename,
        res_km = RES_KM) {
    
    cli_alert_info(paste0("Pre-processing of raster..."))

    # Creating the paths that will be needed in the function
    shapefile_hexagons_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_hexagons-res", res_km, "km_", PROJECT_CRS_TAG, ".gpkg")


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
        
        . <- project_raster_to_hexagons(
                raster = raster, 
                save_to = shapefile_hexagons_path, 
                res_km = res_km) 
        cli_alert_success("Saved shapefile and CSV legend!\n\n")
    }

    return(shapefile_hexagons_path)
}

# A function to standardise hexagonal projection and save path.
# ARGS:
#   - shp_path: filepaths to a shapefile (.gpkg) file.
#   - layer_name: a string.
#   - save_to_basename: a basename (path without extension) to save the file.
#   - res_km: an integer/a float. Grid resolution in km.
# RETURNS:
#   - The reprojected SpatVector
shp_standardised_hexagonal_projection <- function(
        shp_path,
        layer_name,
        save_to_basename,
        res_km = RES_KM) {
    
    cli_alert_info(paste0("Pre-processing of vector..."))

    # Creating the paths that will be needed in the function
    shapefile_hexagons_path <- paste0(
        file_path_sans_ext(save_to_basename), 
        "_hexagons-res", res_km, "km_", PROJECT_CRS_TAG, ".gpkg")


    if (!authorise_overwrite(shapefile_hexagons_path)) {
        # Check if file exists
        cli_alert_warning("Skipping pre-processing.\n\n")

    } else {
        cli_alert_info("Loading vector...")
        # import and reproject to a standardized grid of hexagons
        if (!inherits(shp_path, "SpatVector")) {
            shp <- vect(shp_path)
        } else {
            shp <- shp_path
        }
        
        . <- project_shp_to_hexagons(
                vector = shp[, layer_name], 
                save_to = shapefile_hexagons_path, 
                res_km = res_km) 
        cli_alert_success("Saved shapefile!\n\n")
    }

    return(shapefile_hexagons_path)
}
