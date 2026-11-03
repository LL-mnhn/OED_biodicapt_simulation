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


##### Local functions ##### ---------------------------------------------------
check_raster_standardised_raster_projection <- function(
        raster_path, save_to_basename) {
    if (!authorise_overwrite(save_to_basename)) {
        # Check if file exists
        cli_alert_warning("Skipping STOC pre-processing.\n\n")
    } else {
        raster_standardised_raster_projection(
            raster_path = raster_path,
            save_to_basename = save_to_basename)
    }
}
         
preprocess_stoc_dataset <- function(
    cols_to_remove = c(
        "NDVI", "light_pollution", "precip_spring", "tmp_spring", "p_type", 
        "p_milieu", "altitude", "departement", "qualite_inventaire_stoc", 
        "habitat_principal", "habitat_secondaire", "foret_p", "agricole_p", 
        "urbain_p", "ouvert_p", "Habitat_V2", "test$shape_df.short_name", 
        "Total", "Landscape.PCA1", "Landscape.PCA2", "Landscape.PCA3", 
        "Landscape.PCA4", "annee.x", "longitude_wgs84", "latitude_wgs84",
        "longitude_grid_wgs84", "latitude_grid_wgs84")   
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

        # Rename coordinates columns
        stoc_subset <- stoc_subset |>
            dplyr::rename(LON = longitude_wgs84, LAT = latitude_wgs84)

        # remove useless columns      
        stoc_subset <- stoc_subset |> select(-any_of(cols_to_remove))
          
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

    # Save with buffer 
    clip_raster_to_france_shape(
        over_zero_raster, save_to = ELEVATION_PATH_PREPROCESSED)
}

preprocess_light_pollution  <- function() {
    # Pre-cropping
    cropped_pollution <- clip_raster_to_france_shape(
        rast(LIGHT_POLLUTION_FILEPATH),  
        buffer = BUFFER_SIZE*1.25,  # because not reprojected yet
        fill_missing = FALSE,       # wjole world, nothing is missing.
        save_to = NULL)

    # Project to PROJECT_CRS
    cropped_projected_pollution <- project(cropped_pollution, PROJECT_CRS)
    
    # Saving 
    clip_raster_to_france_shape(
        cropped_projected_pollution, 
        fill_missing = FALSE,
        save_to = LIGHT_POLLUTION_PATH_PREPROCESSED)
}

show_save_results <- function() {
    cli_alert_info("Showing plots to the user.")
    cli_alert_warning(paste0(
        "All figures are automatically saved as .pdf inside '", 
        MAPS_PATH, "' ."))
    dir.create(MAPS_PATH, recursive = TRUE)

    # 1. Corine Land Cover Dataset
    # shapefile  
    corine_plot_bis <- ggplot_categorical_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        rast = clip_raster_to_france_shape(
            rast(CORINE_PATH_PREPROCESSED), 
            buffer = 0),
        layer_name = "NEW_LABEL3")
    print(corine_plot_bis)
    standardised_ggplot_save(
        figure = corine_plot_bis, 
        save_path = file.path(MAPS_PATH, "corine.pdf"))   
    
    # 2. chelsa datasets
    for (i in 1:length(CHELSA_DATASETS)) {
        r <- rast(file.path(CHELSA_PATHS_PREPROCESSED[i]))
        chelsa_plot_bis <- ggplot_quantitative_raster_on_background_map(
            background_map = ggplot_get_france_base_map("national"), 
            raster = clip_raster_to_france_shape(r, buffer = 0),
            layer_name = "mean",
            transform_scale = list(
                                    "identity", 
                                    "sqrt", 
                                    mirrored_log_trans(offset_high = 294)
                              )[[i]],
            unit = CHELSA_UNITS[i],
            limits = NULL,
            precision_auto_limits = 0.1)
        print(chelsa_plot_bis)
        standardised_ggplot_save(
            figure = chelsa_plot_bis, 
            save_path = file.path(MAPS_PATH, 
                paste0("chelsa_", CHELSA_DATASETS[i], ".pdf")))    
    }

    # 3. SWF
    r <- rast(file.path(SWF_PATH_PREPROCESSED))
    swf_plot_bis <- ggplot_quantitative_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        raster = clip_raster_to_france_shape(r, buffer = 0),
        transform_scale = "identity",
        layer_name = "swf_2018_100m_03035_v1_0",
        unit = "%",
        limits = NULL,
        precision_auto_limits = 0.1)
    print(swf_plot_bis)
    standardised_ggplot_save(
        figure = swf_plot_bis, 
        save_path = file.path(MAPS_PATH, "swf.pdf"))   

    # 4. NDVI
    r <- rast(file.path(NDVI_PATH_PREPROCESSED))
    ndvi_plot_bis <- ggplot_quantitative_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        raster = clip_raster_to_france_shape(r*100, buffer = 0),
        transform_scale = "identity",
        layer_name = "mean",
        unit = "%",
        limits = NULL,
        precision_auto_limits = 0.1)
    print(ndvi_plot_bis)
    standardised_ggplot_save(
        figure = ndvi_plot_bis, 
        save_path = file.path(MAPS_PATH, "ndvi.pdf"))   

    # 5. Elevation
    r <- rast(file.path(ELEVATION_PATH_PREPROCESSED))
    elev_plot_bis <- ggplot_quantitative_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        raster = clip_raster_to_france_shape(r, buffer = 0),
        transform_scale = "sqrt",
        layer_name = "fileef5d34834802",
        unit = "m",
        limits = NULL,
        precision_auto_limits = 0.1)
    print(elev_plot_bis)
    standardised_ggplot_save(
        figure = elev_plot_bis, 
        save_path = file.path(MAPS_PATH, "elevation.pdf"))   

    # light_pollution
    r <- rast(file.path(LIGHT_POLLUTION_PATH_PREPROCESSED))
    light_plot_bis <- ggplot_quantitative_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        raster = clip_raster_to_france_shape(r, buffer = 0),
        transform_scale = "sqrt",
        layer_name = "harmonized_dn_ntl_2018_simviirs",
        unit = "DN",
        limits = NULL,
        precision_auto_limits = 0.1)
    print(light_plot_bis)
    standardised_ggplot_save(
        figure = light_plot_bis, 
        save_path = file.path(MAPS_PATH, "light_pollution.pdf"))   
    
    # pesticides
    r <- rast(file.path(PESTICIDE_PATH_PREPROCESSED))
    pest_plot_bis <- ggplot_quantitative_raster_on_background_map(
        background_map = ggplot_get_france_base_map("national"), 
        raster = clip_raster_to_france_shape(r, buffer = 0),
        transform_scale = "sqrt",
        layer_name = "mean_tii",
        unit = "eTII",
        limits = NULL,
        precision_auto_limits = 0.1)
    print(pest_plot_bis)
    standardised_ggplot_save(
        figure = pest_plot_bis, 
        save_path = file.path(MAPS_PATH, "pesticides.pdf"))   

    # X. stoc dataset
    stoc_df <- vect(STOC_OBS_FULL)
    cli_alert_warning(paste("Stoc dataset: sampling locations"))

    stoc_plot <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("regional"), 
        df = as.data.frame(stoc_df), 
        lon_c = "lon",
        lat_c = "lat",
        col = PALETTE[3],
        legend_title = "Sampling locations of STOC ")
    print(stoc_plot)
    standardised_ggplot_save(
        figure = stoc_plot, 
        save_path = file.path(
            MAPS_PATH, 
            "stoc_sampling_locations_.pdf"))

    # XX. biodicapt dataset
    biodicapt_df <- vect(BIODICAPT_OBS_FULL)
    cli_alert_warning(paste("Biodicapt dataset: sampling locations"))

    biodicapt_plot <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("regional"), 
        df = as.data.frame(biodicapt_df), 
        lon_c = "lon",
        lat_c = "lat",
        col = PALETTE[1],
        legend_title = "Sampling locations of BIODICAPT")
    print(biodicapt_plot)
    standardised_ggplot_save(
        figure = biodicapt_plot, 
        save_path = file.path(
            MAPS_PATH, 
            "biodicapt_sampling_locations_.pdf"))

    # XXX. 500 ENI dataset
    eni500_df <- vect(ENI500_OBS_FULL)
    cli_alert_warning(paste("500 ENI dataset: sampling locations"))

    eni500_plot <- ggplot_categorical_df_on_background_map(
        background_map = ggplot_get_france_base_map("regional"), 
        df = as.data.frame(eni500_df), 
        lon_c = "lon",
        lat_c = "lat",
        col = PALETTE[2],
        legend_title = "Sampling locations of 500 ENI")
    print(eni500_plot)
    standardised_ggplot_save(
        figure = eni500_plot, 
        save_path = file.path(
            MAPS_PATH, 
            "eni500_sampling_locations_.pdf"))
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


##### Transformation of raw datasets ##### ------------------------------------
if (exists("BLUR_SEED")) {
    # 1. Corine Land Cover dataset
    # WARNING: using "near" and not "mode" to gain some time, might not be as 
    # precise but categorical data will not be re-used in model, it is just for
    # visualisation purposes.
    cli_alert_info("Preprocessing CLC dataset...")
    . <- save_simplified_clc(
        raster_path = CORINE_FILEPATH,
        save_to = CORINE_PATH_PREPROCESSED)
    . <- check_raster_standardised_raster_projection(
        raster_path = CORINE_PATH_PREPROCESSED,
        save_to_basename = CORINE_PROJECTION_BASEPATH)

    # 2. Chelsa dataset 
    cli_alert_info("Preprocessing CHELSA dataset...")
    for (i in 1:length(CHELSA_DATASETS)) {
        cli_alert_info(paste("Preprocessing", CHELSA_DATASETS[i], "dataset..."))
        . <- average_monthly_quantitative_rasters(
            folder_path = CHELSA_DATASETS_PATHS[i],
            save_to = CHELSA_PATHS_PREPROCESSED[i])
        . <- check_raster_standardised_raster_projection(
            raster_path = CHELSA_PREPROCESSED_PATHS[i],
            save_to_basename = CHELSA_PROJECTION_BASEPATHS[i])
    }

    # 3. Stoc dataset
    cli_alert_info("Preprocessing STOC dataset...")
    . <- preprocess_stoc_dataset()

    # 4. Small Woody Features dataset
    cli_alert_info("Preprocessing SWF dataset...")
    # Add buffer
    . <- clip_raster_to_france_shape(
            subst(rast(SWF_FILEPATH), NA, 0), 
            fill_missing = TRUE,
            save_to = SWF_PATH_PREPROCESSED)
    . <- check_raster_standardised_raster_projection(
        raster_path = SWF_PATH_PREPROCESSED,
        save_to_basename = SWF_PROJECTION_BASEPATH)    

    # 5. NDVI dataset
    cli_alert_info("Preprocessing NDVI dataset...")
    . <- average_monthly_quantitative_rasters(
        folder_path = NASA_NDVI_FOLDER,
        save_to = NDVI_PATH_PREPROCESSED)
    . <- check_raster_standardised_raster_projection(
        raster_path = NDVI_PATH_PREPROCESSED,
        save_to_basename = NDVI_PROJECTION_BASEPATH)    

    # 6. Elevation dataset
    cli_alert_info("Preprocessing Elevation dataset...")
    . <- save_elevation_over_zero()
    . <- check_raster_standardised_raster_projection(
        raster_path = ELEVATION_PATH_PREPROCESSED,
        save_to_basename = ELEVATION_PROJECTION_BASEPATH)    
    
    # 7. Light pollution dataset
    cli_alert_info("Preprocessing Light pollution dataset...")
    . <- preprocess_light_pollution()
    . <- check_raster_standardised_raster_projection(
        raster_path = LIGHT_POLLUTION_PATH_PREPROCESSED,
        save_to_basename = LIGHT_POLLUTION_PROJECTION_BASEPATH)    

    # 8. Pesticides use dataset
    cli_alert_info("Preprocessing pesticides use dataset...")
    . <- rasterise_to_france_grid(
            vector = vect(PESTICIDE_FILEPATH)["mean_tii"], 
            layer_name = "mean_tii",
            save_to = PESTICIDE_PATH_PREPROCESSED,  
            res_km = 0.1) 
    . <- check_raster_standardised_raster_projection(
        raster_path = PESTICIDE_PATH_PREPROCESSED,
        save_to_basename = PESTICIDE_PROJECTION_BASEPATH)


##### Assemble rasters ##### --------------------------------------------------
    cli_alert_info("Assembling rasters/shp together")
    . <- create_master_files_features(FEATURES_RAW_FULL, mode = "raw")
    . <- create_master_files_features(FEATURES_RES_FULL, mode = "raster")
    cli_alert_success("Multi-layer rasters/shp ready!")

    
##### Add features to df ##### ------------------------------------------------
    cli_alert_info(paste0("Extracting features (res_mode = ", RES_MODE, ")."))
    . <- save_features_from_obs(
        file = import_biodicapt_csv(),
        raster_path = FEATURES_RAW_FULL,
        save_to = BIODICAPT_OBS_FULL,
        blur = TRUE)

    . <- save_features_from_obs(
        file = import_eni500_csv(),
        raster_path = FEATURES_RAW_FULL,
        save_to = ENI500_OBS_FULL,
        blur = TRUE)

    . <- save_features_from_obs(
        file = STOC_PATH_PREPROCESSED,
        raster_path = FEATURES_RAW_FULL,
        save_to = STOC_OBS_FULL,
        blur = FALSE)

} else {
    cli_alert_warning(paste0(
        "Variable 'BLUR_SEED' is not available, ",
        "you are likely running this script without raw datasets available...",
        "In this case, pre-processing is not necessary. Skipping."))
}

##### Check results ##### -----------------------------------------------------
write_maps <- authorise_overwrite(MAPS_PATH)
if (write_maps) {
    unlink(MAPS_PATH, recursive = TRUE)
    dir.create(MAPS_PATH)
    . <- suppressWarnings(show_save_results())
} else {
    cli_alert_warning("Skipping.\n\n")
}

cli_alert_success("File '1-pre_processing.R' finished running!\n\n")