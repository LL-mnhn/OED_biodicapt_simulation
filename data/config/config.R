##### Global Parameters #####
# These parameters are 'hidden': they can be modified, but you should not 
# change them without a *very* good reason. 
#
# The parameters grouped in this file are re-used by several scripts, having 
# them all in the same place allows for easier access, control and 
# harmonization.
library(ggplot2)
library(colorspace)


### Global variables ----------------------------------------------------------
REMOTE_USERNAME <- "lehnhofl"
OBS_YEAR <- "2018"  # Most recent year for CLC dataset
RES_KM <- 5         # best compromise I think
# Values for map extent (France limits with buffer, in CRS 4326)
LON_MIN <- -5
LON_MAX <- 9.55
LAT_MIN <- 41.35
LAT_MAX <- 51.05


### Global paths --------------------------------------------------------------
WORK_DIR <- file.path(".")
RAW_DATA_PATH <- file.path(WORK_DIR, "data", "raw_data")
PROCESSED_DATA_PATH <- file.path(WORK_DIR, "data", "preprocessed_data")
MAPS_PATH <- file.path(WORK_DIR, "outputs", "maps")
EXPLORE_PATH <- file.path(WORK_DIR, "outputs", "data_exploration")
FIGURES_PATH <- file.path(WORK_DIR, "outputs", "figures")
SIMULATE_PATH <- file.path(WORK_DIR, "outputs", "simulations")
RESULTS_PATH <- file.path(WORK_DIR, "outputs", "results")
ANALYSES_PATH <- file.path(WORK_DIR, "outputs", "results_analyses")


### File paths ----------------------------------------------------------------
## Raw data
# BIODICAPT
BIODICAPT_FOLDER <- file.path(RAW_DATA_PATH, "BIODICAPT_pos")
BIODICAPT_FILENAME <- "parcelles_selectionnees.csv"

# 500 ENI
ENI500_FOLDER <- file.path(RAW_DATA_PATH, "500ENI_pos")
ENI500_FILENAME <- "table_parcelle_NA_excluded.csv"

# STOC
STOC_FILEPATH <- file.path(RAW_DATA_PATH, "STOC", "Dataframe_STOC_input.RData")
STOC_PATH_PREPROCESSED <- file.path(
    PROCESSED_DATA_PATH, "STOC_data_cleaned.csv") 

# CORINE
CORINE_FILEPATH <- file.path(RAW_DATA_PATH, 
    paste0("u", OBS_YEAR, "_clc", OBS_YEAR, "_v2020_20u1_raster100m"), "DATA", 
    paste0("U", OBS_YEAR, "_CLC", OBS_YEAR, "_V2020_20u1.tif"))
CORINE_PROJECTION_BASEPATH <- file.path(
    PROCESSED_DATA_PATH, paste0("CLC", OBS_YEAR))
CORINE_PATH_PREPROCESSED <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("original_CLC", OBS_YEAR, "_simplified_categories.tif"))

# CHELSA
CHELSA_DATASETS <- c("hurs", "pr", "tas") # see https://www.chelsa-climate.org/datasets/chelsa_monthly
CHELSA_UNITS <- c("%", expression(paste("kg m"^{-2}, " month"^{-1})), "K")
CHELSA_FOLDER <- file.path(RAW_DATA_PATH, "CHELSA-monthly")
CHELSA_DATASETS_PATHS <- file.path(CHELSA_FOLDER, CHELSA_DATASETS)
CHELSA_PROJECTION_BASEPATHS <- file.path(
    PROCESSED_DATA_PATH, paste0("CHELSA_", OBS_YEAR, "_", CHELSA_DATASETS))
CHELSA_ANNUAL_PATHS <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("CHELSA_", OBS_YEAR, "_", CHELSA_DATASETS,  "_annual-WGS84.tif"))

# NDVI
NASA_NDVI_FOLDER <- file.path(RAW_DATA_PATH, "NASA-NDVI")
NDVI_PROJECTION_BASEPATH <- file.path(
        PROCESSED_DATA_PATH,
        paste0("NDVI_", OBS_YEAR))
NDVI_ANNUAL_PATH <- file.path(
        PROCESSED_DATA_PATH,
        paste0("NDVI_", OBS_YEAR, "_annual-WGS84.tif"))

# Elevation
ELEVATION_ZOOM_LEVEL <- 7
ELEVATION_FOLDER <- file.path(RAW_DATA_PATH, "ELVATR")
ELEVATION_FILEPATH <- file.path(
    ELEVATION_FOLDER, paste0("elevation_zoom", ELEVATION_ZOOM_LEVEL, ".tif"))
ELEVATION_PROJECTION_BASEPATH <- file.path(
    PROCESSED_DATA_PATH, "Elevation_min-sea-level")
ELEVATION_PATH_PREPROCESSED <- file.path(
    PROCESSED_DATA_PATH, "Elevation_min-sea-level.tif") 

# Light polution
LIGHT_POLLUTION_FOLDER <- file.path(RAW_DATA_PATH, "light_pollution")
LIGHT_POLLUTION_FILEPATH <- file.path(
    LIGHT_POLLUTION_FOLDER, 
    paste0("Harmonized_DN_NTL_", OBS_YEAR, "_simVIIRS.tif"))
LIGHT_POLLUTION_PROJECTION_BASEPATH <- file.path(
    PROCESSED_DATA_PATH,  paste0("light_pollution_", OBS_YEAR))


## Final data files to use
BIODICAPT_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("anonymised_BIODICAPT_", OBS_YEAR, "_obs_features_res", RES_KM, 
    "km-WGS84.gpkg"))
ENI500_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("anonymised_ENI500_", OBS_YEAR, "_obs_features_res", RES_KM, 
    "km-WGS84.gpkg"))
STOC_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("STOC_", OBS_YEAR, "_obs_features_res", RES_KM, "km-WGS84.gpkg"))
FEATURES_RAW_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("RAW_", OBS_YEAR, "_env-data_nokm-WGS84.tif"))
FEATURES_RES_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("RASTER_", OBS_YEAR, "_env-data_", RES_KM, "km-WGS84.tif"))
FEATURES_HEX_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("HEX_", OBS_YEAR, "_env-data_", RES_KM, "km-WGS84.gpkg"))


# Species names
NAMES_SPECIES <- c(
    "Alauda_arvensis", "Anthus_trivialis", "Carduelis_cannabina",
    "Carduelis_carduelis", "Carduelis_chloris", "Certhia_brachydactyla",  
    "Columba_palumbus", "Corvus_corone", "Cyanistes_caeruleus",
    "Dendrocopos_major", "Emberiza_cirlus", "Emberiza_citrinella",
    "Erithacus_rubecula", "Fringilla_coelebs", "Garrulus_glandarius",
    "Hippolais_polyglotta", "Hirundo_rustica", "Luscinia_megarhynchos",
    "Motacilla_alba", "Parus_major", "Passer_domesticus", "Periparus_ater",
    "Phoenicurus_ochruros", "Phylloscopus_collybita", "Pica_pica",
    "picus_viridis", "prunella_modularis", "Regulus_ignicapilla", 
    "Saxicola_rubicola", "Serinus_serinus", "Sitta_europaea", 
    "Streptopelia_decaocto", "Streptopelia_turtur", "Sturnus_vulgaris", 
    "Sylvia_atricapilla", "Sylvia_communis", "Troglodytes_troglodytes", 
    "Turdus_merula", "Turdus_philomelos", "Turdus_viscivorus"     
)
                  

### Simulations and Models ----------------------------------------------------
K_FOLDS <- 10

                       
### Plot styling --------------------------------------------------------------
FONT <- ifelse((Sys.info()["user"] == "lehnhofl")[[1]], "sans", "Lexend")
PALETTE <- c("#D9054E", "#28A349", "#246CBC", 
             "#5D7B84", "#C2562F", "#FFB703",
             "#7B2CBF", "#333333", "#0FA3B1", "#890f40")

SHAPES <- c(21, 22, 23, 24, 25, 8, 13, 7)
SIZES <- c(1.66, 1.85, 1.66, 1.5, 1.5, 1, 1, 1)
STROKES <- c(0, 0, 0, 0, 0, 0.5, 0.5, 0.5)
CUSTOM_SCALES <- list(
    scale_color_manual(values = darken(PALETTE, amount = 0.66)),
    scale_fill_manual(values = PALETTE),
    scale_shape_manual(values = SHAPES),
    scale_size_manual(values = SIZES),
    scale_discrete_manual(aesthetics = "stroke", values = STROKES)
)
LIGHT_CUSTOM_SCALES <- list(
    scale_color_manual(values = PALETTE),
    scale_fill_manual(values = PALETTE),
    scale_shape_manual(values = SHAPES),
    scale_size_manual(values = SIZES),
    scale_discrete_manual(aesthetics = "stroke", values = STROKES)
)