# This scripts checks the presence of files necessary for scripts 
# in ./data/raw_data

##### Libraries ##### ---------------------------------------------------------
library(cli)
library(terra)
library(elevatr)


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data", "config", "config.R"))) # Global parameters


##### Local functions ##### ---------------------------------------------------
check_biodicapt_files <- function() {
    cli_alert_info("Checking BIODICAPT files...")

    # The folder...
    if (!file.exists(BIODICAPT_FOLDER)) {
        cli_alert_danger(
            paste0("Folder '", BIODICAPT_FOLDER, "' was not found!\n"))
        stop("Tip: did you create the necessary folders for the files?")
    }
    
    # ...must contain 1 file with the following name :
    if (!file.exists(file.path(BIODICAPT_FOLDER, BIODICAPT_FILENAME))) {
        cli_alert_danger(paste0(
            "File '", BIODICAPT_FILENAME, 
            "' was not found in '", BIODICAPT_FOLDER,"'!\n"))
        cli_alert_warning(paste(
            "This dataset can only be obtained from the coordinators",
            "of the BIODICAPT project. Contact them directly."))
        stop("Tip: these files must be downloaded manually.")
    }

    cli_alert_success("Raw BIODICAPT file is available.\n\n")
}

check_eni500_files <- function() {
    cli_alert_info("Checking 500 ENI file...")

    # The folder...
    if (!file.exists(ENI500_FOLDER)) {
        cli_alert_danger(paste0("Folder '", ENI500_folder, 
            "' was not found!\n"))
        stop("Tip: did you create the necessary folders for the files?")
    }
    
    # ...must contain 1 file with the following name :
    if (!file.exists(file.path(ENI500_FOLDER, ENI500_FILENAME))) {
        cli_alert_danger(paste0(
            "File '", ENI500_FILENAME, 
            "' was not found in '", ENI500_FOLDER,"'!\n"))
        cli_alert_warning(paste0(
            "This dataset can only be obtained from the coordinators ", 
            "of the 500 ENI network. Contact them directly."))
        stop("Tip: this file must be downloaded manually (the version with NAs removed).")
    }
    
    cli_alert_success("Raw 500 ENI file is available.\n\n")
}

check_corine_raster <- function() {
    cli_alert_info("Checking CORINE Land Cover files...")

    # Check if file exists
    if (!file.exists(CORINE_FILEPATH)) {
        cli_alert_danger(
            paste0("File '", CORINE_FILEPATH, "' was not found!\n"))
        cli_alert_warning("If needed, download file from https://www.data.gouv.fr/datasets/corine-land-cover-edition-2018-france-metropolitaine")
        stop("Tip: did you create the necessary folders for the files?")
    }

    cli_alert_success("Raw CLC files are available.\n\n")
}

check_chelsa_rasters <- function() {
    cli_alert_info("Checking CHELSA files...")

    # The folder...
    if (!file.exists(CHELSA_FOLDER)) {
        cli_alert_danger(paste0("Folder '", CHELSA_FOLDER, 
            "' was not found! Creating it...\n"))
        dir.create(CHELSA_FOLDER)
    }

    # ... must contain one subfolder per dataset
    for (chelsa_dataset in CHELSA_DATASETS) {
        if (!file.exists(file.path(CHELSA_FOLDER, chelsa_dataset))) {
            cli_alert_danger(paste0(
                "Folder '", 
                file.path(CHELSA_FOLDER, chelsa_dataset), 
                "' was not found! Creating it...\n"))
            dir.create(file.path(CHELSA_FOLDER, chelsa_dataset))
        }
        
        files <- list.files(file.path(CHELSA_FOLDER, chelsa_dataset))
        if (length(files) != 12) {
            stop(paste0("Folder ", file.path(CHELSA_FOLDER, chelsa_dataset), 
            "contains ", length(files), 
            " files. There should be 12 files (one per month).\n",
            "If needed, files can be downloaded at https://appeears.earthdatacloud.nasa.gov"))
        }

        for (month in 1:12) {
            filename <- paste0(
                "CHELSA_", chelsa_dataset, 
                "_", sprintf("%02d", month),
                "_", OBS_YEAR, ".tif")
            filepath <- file.path(CHELSA_FOLDER, chelsa_dataset, filename)
            if (!file.exists(filepath)) {
                cli_alert_info(paste0(
                    "Missing CHELSA file, auto-download of file ", 
                    filename, "..."))
                date <- as.Date(paste0(OBS_YEAR, "-", sprintf('%02d', month), "-01"))
                chelsa_raster <- getChelsa(
                    chelsa_dataset, 
                    extent = c(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX), 
                    startdate = date, enddate = date,
                    dataset = "chelsa-monthly")

                writeRaster(chelsa_raster, filepath)
            }
        }
    }

    cli_alert_success("Raw CHELSA files are available.\n\n")
}

check_appears_rasters <- function() {
    cli_alert_info("Checking NASA NDVI files...")

    # The folder...
    if (!file.exists(NASA_NDVI_FOLDER)) {
        cli_alert_danger(paste0("Folder '", NASA_NDVI_FOLDER, 
            "' was not found! Creating it...\n"))
        dir.create(NASA_NDVI_FOLDER)
    }

    # ... must contain 12 files following a pattern
    files <- list.files(NASA_NDVI_FOLDER)
    if (length(files) != 12) {
        stop(paste0("Folder ", NASA_NDVI_FOLDER, "contains ", length(files), 
        " files. There should be 12 files (one per month).\n",
        "If needed, files can be downloaded at https://appeears.earthdatacloud.nasa.gov"))
    }
    pattern <- paste0("MOD13A3.061__1_km_monthly_NDVI_doy", OBS_YEAR, "*.tif")
    matches <- grepl(glob2rx(pattern), files)
    if (!all(matches)) {
        stop(paste0(
            "Files in ", NASA_NDVI_FOLDER, 
            " do not follow the expected pattern: ", pattern))
    }
    cli_alert_success("Raw NASA NDVI files are available.\n\n")
    
}

check_elevation_raster <- function(){
    cli_alert_info("Checking elevation files...")

    # The folder...
    if (!file.exists(ELEVATION_FOLDER)) {
        cli_alert_danger(paste0("Folder '", ELEVATION_FOLDER, 
            "' was not found! Creating it...\n"))
        dir.create(ELEVATION_FOLDER)
    }

    # Check if file exists
    if (!file.exists(ELEVATION_FILEPATH)) {
        cli_alert_danger(
            paste0("File '", ELEVATION_FILEPATH, "' was not found!\n"))
        cli_alert_warning("Auto-download...")
        extent_df <- data.frame(
            x = c(LON_MIN, LON_MAX), y = c(LAT_MIN, LAT_MAX))
        raster_elevation <- rast(get_elev_raster(
            locations = extent_df,
            z = ELEVATION_ZOOM_LEVEL,       
            prj = "EPSG:4326",
            clip = "bbox"))
        writeRaster(raster_elevation, ELEVATION_FILEPATH)
    }
    cli_alert_success("Raw elevation file is available.\n\n")
}

check_light_pollution_raster <- function() {
    cli_alert_info("Checking light pollution file...")

    # The folder...
    if (!file.exists(LIGHT_POLLUTION_FOLDER)) {
        cli_alert_danger(paste0("Folder '", LIGHT_POLLUTION_FOLDER, 
            "' was not found! Creating it...\n"))
        dir.create(LIGHT_POLLUTION_FOLDER)
    }

    # Check if file exists
    if (!file.exists(LIGHT_POLLUTION_FILEPATH)) {
        cli_alert_danger(
            paste0("File '", LIGHT_POLLUTION_FILEPATH, "' was not found!\n"))
        stop("Tip: you can download this file at https://doi.org/10.6084/m9.figshare.9828827.v10")
    }

    cli_alert_success("Raw light pollution file is available.\n\n")
}

check_species_data <- function() {
    cli_alert_info("Checking STOC data file...")

    # Check if file exists
    if (!file.exists(STOC_FILEPATH)) {
        cli_alert_danger(paste0("File '", STOC_FILEPATH, "' was not found!\n"))
        cli_alert_warning("If needed, download file from https://doi.org/10.5061/dryad.bnzs7h4g3.")
        stop("Tip: did you create the necessary folders for the files?")
    }

    cli_alert_success("Raw STOC dataframe is available.\n\n")
}


##### Verify each dataset ##### -----------------------------------------------
if (length(list.files(RAW_DATA_PATH)) == 0) {
    cli_alert_danger(paste0(
        "No files were found in", RAW_DATA_PATH, "."))
    cli_alert_warning(paste0(
        "↳ If you cloned this repository from GitHub, ",
        "ignore this script and begin usage with `1-pre_processing.R`\n\n"))
} else {
    # 1. BIODICAPT dataset
    . <- check_biodicapt_files()

    # 2. 500 ENI dataset
    . <- check_eni500_files()

    # 3. CORINE Land Cover dataset
    . <- check_corine_raster()

    # 4. CHELSA datasets
    . <- check_chelsa_rasters()

    # 5. NDVI dataset
    . <- check_appears_rasters()

    # 6. Elevation dataset
    . <- check_elevation_raster()

    # 7. Light pollution dataset
    . <- check_light_pollution_raster()

    # 8. Species absence-presence dataset
    . <- check_species_data()
}


