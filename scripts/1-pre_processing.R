# This scripts:
#   - groups datasets splitted in different files
#   - formats datasets into similar csvs/rasters
#   - anonymizes gps coordinates

##### Libraries ##### ---------------------------------------------------------
library(tidyterra)
library(terra)
library(dplyr)
library(readr)
library(tools)
library(cli)
library(sf)

suppressPackageStartupMessages(source(here::here(file.path("R", "utils_data.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_figures.R"))))
suppressPackageStartupMessages(source(here::here(file.path("R", "utils_analysis.R"))))


##### Parameters ##### --------------------------------------------------------
source(here::here(file.path("data", "config", "config.R"))) # Global parameters
if (file.exists(file.path("data", "config", "seed.R"))) {
    # seed is hidden for confidentiality of data points
    source(here::here(file.path("data", "config", "seed.R"))) 
}

RES_MODE <- "raw"


##### Local functions ##### ---------------------------------------------------
quick_standardised_projections <- function(path_raster, basename_path, mode = "mean") {
    . <- standardized_raster_projection(
        raster_path = path_raster,
        save_to_basename = basename_path)
    . <- standardized_hexagonal_projection(
        raster_path = path_raster,
        save_to_basename = basename_path,
        extraction_mode = mode)    
}

preprocess_stoc_dataset <- function(
    cols_to_remove = c(
        "NDVI", "light_pollution", "precip_spring", "tmp_spring", "p_type", 
        "p_milieu", "altitude", "departement", "qualite_inventaire_stoc", 
        "habitat_principal", "habitat_secondaire", "foret_p", "agricole_p", 
        "urbain_p", "ouvert_p", "Habitat_V2", "test$shape_df.short_name", 
        "Total", "Landscape.PCA1", "Landscape.PCA2", "Landscape.PCA3", 
        "Landscape.PCA4", "annee.x")   
    ) {
    cli_alert_info("Pre-processing of the STOC dataset.")

    if (!authorise_overwrite(STOC_PATH_PREPROCESSED)) {
        # Check if file exists
        cli_alert_warning("Skipping STOC pre-processing.\n\n")

    } else {
        # preprocess file if authorized
        cli_alert_info("Loading STOC dataframe...")

        # loading .Rdata --> variable name: STOC
        load(STOC_FILEPATH)
        
        # select year (STOC contains years 2015-2018)
        if (OBS_YEAR %in% unique(STOC$annee.x)){
            stoc_subset <- subset(STOC, annee.x == OBS_YEAR)
        } else {
            stop(paste0(
                "Selected year '", OBS_YEAR, 
                "' is not in the loaded STOC dataset"))
        }
        
        # convert abundance to occurence
        stoc_subset <- stoc_subset |> 
            mutate(across(23:64, ~ 1L * (. > 0))) 
        # remove useless columns      
        stoc_subset <- stoc_subset |> select(-any_of(cols_to_remove))

        # Rename coordinates columns
        stoc_subset <- stoc_subset |> 
            rename(LON = longitude_wgs84, LAT = latitude_wgs84)
          
        # No blurring necessary (this is a public dataset)
        # Save resulting dataframe
        cli_alert_info("Saving file...")
        write_csv(stoc_subset, STOC_PATH_PREPROCESSED)
        cli_alert_success("Modified STOC file saved!\n\n")
    }
}

save_elevation_over_zero <- function() {
    # Get original file
    base_elevation_raster <- rast(ELEVATION_FILEPATH)

    # apply hard threshold 
    # (we're not interested about whatever is under sea level)
    over_zero_raster <- clamp(base_elevation_raster, lower = 0, values = TRUE)

    # save new raster
    writeRaster(over_zero_raster, ELEVATION_PATH_PREPROCESSED, overwrite = TRUE)
}

show_save_results <- function() {
    cli_alert_info("Showing plots to the user.")
    cli_alert_warning(paste0(
        "All figures are automatically saved as .pdf inside '", 
        MAPS_PATH, "' ."))
    dir.create(MAPS_PATH, recursive = TRUE)

    # I'm only showing hexagonal projections by default
    # If needed, here is an example on how to adapt the function for raster
    # files. (careful though, it is old, I did not check paths in a while)
    # corine_raster <- rast(
    #     paste0(CORINE_BASENAME, "_projection_france_res", 
    #            RES_KM, "km-WGS84.tif")
    # )
    # corine_plot <- ggplot_categorical_raster_on_background_map(
    #     background_map = ggplot_get_france_base_map("national"), 
    #     raster = corine_raster,
    #     layer_name = "NEW_LABEL3")
    # print(corine_plot)
    # standardised_ggplot_save(
    #     figure = corine_plot, 
    #     save_path = file.path(MAPS_PATH, "corine_raster.pdf"))
    
    # 1. Corine Land Cover Dataset
    # shapefile  
    corine_shapefile <- st_read(
        paste0(CORINE_PROJECTION_BASEPATH, 
            "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
        quiet = TRUE)
    corine_shapefile_colors <- read_csv(
        paste0(CORINE_PROJECTION_BASEPATH, 
            "_projection_france_hexagons_res", RES_KM, "km-WGS84_NEW_LABEL3.csv"),
        show_col_types = FALSE)
    corine_plot_bis <- ggplot_categorical_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        shapefile = corine_shapefile,
        layer_name = "NEW_LABEL3",
        color_df = corine_shapefile_colors)
    print(corine_plot_bis)
    standardised_ggplot_save(
        figure = corine_plot_bis, 
        save_path = file.path(MAPS_PATH, "corine_hexagons.pdf"))          
    
    # 2. chelsa datasets
    for (i in 1:length(CHELSA_DATASETS)) {
        # shapefile  
        chelsa_shapefile <- st_read(
            paste0(CHELSA_PROJECTION_BASEPATHS[i], 
                "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
            quiet = TRUE)
        chelsa_plot_bis <- ggplot_quantitative_shapefile_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            shapefile = chelsa_shapefile,
            layer_name = "mean",
            unit=CHELSA_UNITS[i],
            limits=NULL,
            precision_auto_limits = 0.1)
        print(chelsa_plot_bis)
        standardised_ggplot_save(
            figure = chelsa_plot_bis, 
            save_path = file.path(MAPS_PATH, 
                paste0("chelsa_", CHELSA_DATASETS[i], "_hexagons.pdf")))    
    }
    
    # 3. stoc dataset
    stoc_df <- read_csv(STOC_PATH_PREPROCESSED, show_col_types = FALSE)
    cli_alert_warning(paste("Stoc dataset: sampling locations"))

    stoc_plot <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        df = stoc_df, 
        lon_c = "LON",
        lat_c = "LAT",
        legend_title = "Sampling locations of STOC ")
    print(stoc_plot)
    standardised_ggplot_save(
        figure = stoc_plot, 
        save_path = file.path(
            MAPS_PATH, 
            "stoc_sampling_locations_.pdf"))

    # 4. NDVI dataset
    # shapefile  
    ndvi_shapefile <- st_read(
        paste0(NDVI_PROJECTION_BASEPATH, 
            "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
        quiet = TRUE)
    ndvi_plot_bis <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        shapefile = ndvi_shapefile,
        layer_name = "mean",
        unit="Index",
        limits=NULL,
        precision_auto_limits = 0.1)
    print(ndvi_plot_bis)
    standardised_ggplot_save(
        figure = ndvi_plot_bis, 
        save_path = file.path(MAPS_PATH, paste0("NDVI_hexagons.pdf"))) 


    # 5. Elevation dataset
    # shapefile  
    elev_shapefile <- st_read(
        paste0(ELEVATION_PROJECTION_BASEPATH, 
            "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
        quiet = TRUE)
    elev_plot_bis <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        shapefile = elev_shapefile,
        layer_name = "file1254937095365",
        unit="m",
        limits=NULL,
        precision_auto_limits = 0.1)
    print(elev_plot_bis)
    standardised_ggplot_save(
        figure = elev_plot_bis, 
        save_path = file.path(MAPS_PATH, paste0("elevation_hexagons.pdf"))) 

    # 6. Light pollution dataset
    # shapefile  
    light_shapefile <- st_read(
        paste0(LIGHT_POLLUTION_PROJECTION_BASEPATH, 
            "_projection_france_hexagons_res", RES_KM, "km-WGS84.gpkg"),
        quiet = TRUE)
    light_plot_bis <- ggplot_quantitative_shapefile_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        shapefile = light_shapefile,
        layer_name = "Harmonized_DN_NTL_2018_simVIIRS",
        unit="m",
        limits=NULL,
        precision_auto_limits = 0.1)
    print(light_plot_bis)
    standardised_ggplot_save(
        figure = light_plot_bis, 
        save_path = file.path(
            MAPS_PATH, paste0("light_pollution_hexagons.pdf")))    
    
    cli_alert_success("Plots and PDFs are ready!")
}

import_biodicapt_csv <- function() {
    cli_alert_info("Pre-processing of BIODICAPT dataset.")
    # load file and remove useless columns
    biodicapt_df <- read_csv(
        file.path(BIODICAPT_FOLDER, BIODICAPT_FILENAME),
        show_col_types = FALSE)
    
    # Rename coordinates columns
    biodicapt_df <- biodicapt_df |> rename(LON = Longitude, LAT = Latitude)
    cli_alert_success("Dataset is ready!\n\n")

    return(biodicapt_df)
}

import_eni500_csv <- function(
        cols_to_remove = c(
        "lieu_dit", "code_postal", "pourcent_pente", "id_parcelle",
        "commentaire_parcelle", "derniere_modif_parcelle_par", "commune",
        "derniere_modif_parcelle_le", "derniere_modif_donnees_agro_par",
        "derniere_modif_donnees_agro_le", "derniere_modif_pratiques_par",
        "derniere_modif_pratiques_le", "code_parcelle", "nom_parcelle")
    ) {
    cli_alert_info("Pre-processing of 500 ENI dataset.")

    # preprocess file if authorized
    cli_alert_info("Loading 500 ENI files...")

    # load file and remove useless columns
    eni500_df <- read_csv(
        file.path(ENI500_FOLDER, ENI500_FILENAME),
        show_col_types = FALSE)
    eni500_df <- eni500_df |> select(-any_of(cols_to_remove))

    # Rename coordinates columns
    eni500_df <- eni500_df |> rename(LON = X, LAT = Y)
    cli_alert_success("Dataset is ready!\n\n")

    return(eni500_df)
}

save_biodicapt_and_500ENI_map_examples <- function(feature, legend, limits) {
    # with the name of a layer in BIODICAPT_OBS_FULL and ENI500_OBS_FULL
    # shows the associated data points on a map of france.
    see_datasets <- typeline(
        prompt = paste0(paste0(
            "Show processed datasets & Save figures ",
            "(overwrites by default)? [Y/n]: ")))
    cleaned_answer <- tolower(trimws(see_datasets))

    if (cleaned_answer %in% c("y", "yes")) {
        cli_alert_warning(paste0(
            "All figures are automatically saved as .pdf inside '", 
            MAPS_PATH, "' ."))
        cli_alert_info("Showing plots to the user.")

        # BIODICAPT plot
        biodicapt_plot <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map(),
            df = vect(BIODICAPT_OBS_FULL),
            column = feature,
            unit = legend,
            limits = limits)
        print(biodicapt_plot)
        standardised_ggplot_save(
            figure = biodicapt_plot, 
            save_path = file.path(
                MAPS_PATH, paste0("biodicapt_", feature, ".pdf"))
            ) 
        
        # 500 ENI plot
        eni500_plot <- ggplot_quantitative_df_on_background_map(
            background_map = ggplot_get_france_base_map(),
            df = vect(ENI500_OBS_FULL),
            column = feature,
            unit = legend,
            limits = limits)
        print(eni500_plot)
        standardised_ggplot_save(
            figure = eni500_plot, 
            save_path = file.path(
                MAPS_PATH, paste0("eni500_", feature, ".pdf"))
            ) 
        
        cli_alert_success("Plots and PDFs are ready!")
    } else {
        cli_alert_warning("Skipping.\n\n")
    }
}


##### Transformation of raw datasets ##### ------------------------------------
if (exists("BLUR_SEED")) {
    # 1. Corine Land Cover dataset
    . <- save_simplified_clc(
        raster_path = CORINE_FILEPATH,
        save_to = CORINE_PATH_PREPROCESSED)
    . <- quick_standardised_projections(
        path_raster = CORINE_PATH_PREPROCESSED,
        basename_path = CORINE_PROJECTION_BASEPATH)

    # 2. Chelsa dataset 
    for (i in 1:length(CHELSA_DATASETS)) {
        . <- average_monthly_quantitative_rasters(
            folder_path = CHELSA_DATASETS_PATHS[i],
            save_to = CHELSA_ANNUAL_PATHS[i])
        . <- quick_standardised_projections(
            path_raster = CHELSA_ANNUAL_PATHS[i],
            basename_path = CHELSA_PROJECTION_BASEPATHS[i],
            mode = "mean")
    }

    # 3. Stoc dataset
    . <- preprocess_stoc_dataset()

    # 4. NDVI dataset
    . <- average_monthly_quantitative_rasters(
        folder_path = NASA_NDVI_FOLDER,
        save_to = NDVI_ANNUAL_PATH)
    . <- quick_standardised_projections(
        path_raster = NDVI_ANNUAL_PATH,
        basename_path = NDVI_PROJECTION_BASEPATH,
        mode = "mean")    

    # 5. Elevation dataset
    . <- save_elevation_over_zero()
    . <- quick_standardised_projections(
        path_raster = ELEVATION_PATH_PREPROCESSED,
        basename_path = ELEVATION_PROJECTION_BASEPATH,
        mode = "mean")    
    
    # 6. Light pollution dataset
    . <- quick_standardised_projections(
        path_raster = LIGHT_POLLUTION_FILEPATH,
        basename_path = LIGHT_POLLUTION_PROJECTION_BASEPATH,
        mode = "mean")       
} else {
    cli_alert_warning(paste0(
        "Variable 'BLUR_SEED' is not available, ",
        "you are likely running this script without raw datasets available...",
        "In this case, pre-processing is not necessary. Skipping."))
}

### Check results
write_maps <- authorise_overwrite(MAPS_PATH)
if (write_maps) {
    unlink(MAPS_PATH, recursive = TRUE)
    dir.create(MAPS_PATH)
    . <- suppressWarnings(show_save_results())
} else {
    cli_alert_warning("Skipping.\n\n")
}


##### Assemble rasters ##### --------------------------------------------------
cli_alert_info("Assembling rasters/shp together")
. <- create_master_files_features(FEATURES_RAW_FULL, mode = "raw")
. <- create_master_files_features(FEATURES_RES_FULL, mode = "raster")
. <- create_master_files_features(FEATURES_HEX_FULL, mode = "hexagonal")
cli_alert_success("Multi-layer rasters/shp ready!")

##### Add features to df ##### ------------------------------------------------
cli_alert_info(paste0("Extracting features (res_mode = ", RES_MODE, ")."))
. <- save_features_from_obs(
    file = import_biodicapt_csv(), 
    save_to = BIODICAPT_OBS_FULL,
    mode = RES_MODE,
    blur = TRUE)

. <- save_features_from_obs(
    file = import_eni500_csv(), 
    save_to = ENI500_OBS_FULL,
    mode = RES_MODE,
    blur = TRUE)

. <- save_features_from_obs(
    file = STOC_PATH_PREPROCESSED, 
    save_to = STOC_OBS_FULL,
    mode = RES_MODE,
    blur = FALSE)

### Check some results
if (write_maps) {
    . <- save_biodicapt_and_500ENI_map_examples(
    "elevation", "Altitude (m)", limits = c(0, 4800))
    . <- save_biodicapt_and_500ENI_map_examples(
        "chelsa_tas", "T°K", limits = c(270, 290))
}

cli_alert_success("File '1-pre_processing.R' finished running!\n\n")