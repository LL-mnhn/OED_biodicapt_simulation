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
# Values for map extent (France limits with buffer, in PROJECT_CRS)
LON_MIN <- 3200000
LON_MAX <- 4300000
LAT_MIN <- 2000000
LAT_MAX <- 3150000

PROJECT_CRS <- "EPSG:3035"  # Must be a metric projection (either 3035 or 2154)
PROJECT_CRS_TAG  <- gsub(":", "", PROJECT_CRS)   # for file names
FRANCE_BBOX_4326 <- c(xmin = -5.5, ymin = 41.2, xmax = 9.7, ymax = 51.2)


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

# Small Woody Features
SWF_FILEPATH <- file.path(RAW_DATA_PATH, "woody_features", 
                          paste0("SWF_", OBS_YEAR, "_100m_03035_V1_0.tif"))
SWF_PATH_PREPROCESSED <- file.path(
    PROCESSED_DATA_PATH, paste0("SWF_", OBS_YEAR, "_100m_with_buffer.tif"))
SWF_PROJECTION_BASEPATH <- file.path(
    PROCESSED_DATA_PATH, paste0("SWF_", OBS_YEAR, "_100m_with_buffer")) 

# CHELSA
CHELSA_DATASETS <- c("hurs", "pr", "tas") # see https://www.chelsa-climate.org/datasets/chelsa_monthly
CHELSA_UNITS <- c("%", expression(paste("kg m"^{-2}, " month"^{-1})), "K")
CHELSA_FOLDER <- file.path(RAW_DATA_PATH, "CHELSA-monthly")
CHELSA_DATASETS_PATHS <- file.path(CHELSA_FOLDER, CHELSA_DATASETS)
CHELSA_PROJECTION_BASEPATHS <- file.path(
    PROCESSED_DATA_PATH, paste0("CHELSA_", OBS_YEAR, "_", CHELSA_DATASETS))
CHELSA_ANNUAL_PATHS <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("CHELSA_", OBS_YEAR, "_", CHELSA_DATASETS,  "_annual-EPSG4326.tif")) # 4326 is the default CRS for CHELSA data

# NDVI
NASA_NDVI_FOLDER <- file.path(RAW_DATA_PATH, "NASA-NDVI")
NDVI_PROJECTION_BASEPATH <- file.path(
        PROCESSED_DATA_PATH,
        paste0("NDVI_", OBS_YEAR))
NDVI_ANNUAL_PATH <- file.path(
        PROCESSED_DATA_PATH,
        paste0("NDVI_", OBS_YEAR, "_annual-EPSG4326.tif"))

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
    "km-", PROJECT_CRS_TAG, ".gpkg"))
BIODICAPT_OBS_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_anonymised_BIODICAPT_", OBS_YEAR, "_obs_features_res", RES_KM, 
    "km-", PROJECT_CRS_TAG, ".gpkg"))

ENI500_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("anonymised_ENI500_", OBS_YEAR, "_obs_features_res", RES_KM, 
    "km-", PROJECT_CRS_TAG, ".gpkg"))
ENI500_OBS_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_anonymised_ENI500_", OBS_YEAR, "_obs_features_res", RES_KM, 
    "km-", PROJECT_CRS_TAG, ".gpkg"))

STOC_OBS_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("STOC_", OBS_YEAR, "_obs_features_res", RES_KM, "km-", PROJECT_CRS_TAG, ".gpkg"))
STOC_OBS_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_STOC_", OBS_YEAR, "_obs_features_res", RES_KM, "km-", PROJECT_CRS_TAG, ".gpkg"))

FEATURES_RAW_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("RAW_", OBS_YEAR, "_env-data_nokm-", PROJECT_CRS_TAG, ".tif"))
FEATURES_RAW_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_RAW_", OBS_YEAR, "_env-data_nokm-", PROJECT_CRS_TAG, ".tif"))

FEATURES_RES_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("RASTER_", OBS_YEAR, "_env-data_", RES_KM, "km-", PROJECT_CRS_TAG, ".tif"))
FEATURES_RES_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_RASTER_", OBS_YEAR, "_env-data_", RES_KM, "km-", PROJECT_CRS_TAG, ".tif"))

FEATURES_HEX_FULL <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("HEX_", OBS_YEAR, "_env-data_", RES_KM, "km-", PROJECT_CRS_TAG, ".gpkg"))
FEATURES_HEX_TRF <- file.path(
    PROCESSED_DATA_PATH, 
    paste0("vtransform_HEX_", OBS_YEAR, "_env-data_", RES_KM, "km-", PROJECT_CRS_TAG, ".gpkg"))


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
K_FOLDS <- 3

# FYI: 
#   - "clc" is removed because it is a categorical variable
#   - "distance to NODATA" is removed, not informative
#   - "peramanent.crops" is removed as it easily heavily biased against brittany
#   - Layers that are correlated together are removed:
#       - "CHELSA_tas" (correlated with elevation)
#       - "prop_artificial.surfaces_within_7071m" (correlated with light_pollution)
#       - "prop_water.bodies_within_7071m" (correlated with distance.to.water.bodies)
#       - "prop_arable.land_within_7071m" (correlated with distance.to.arable.land)
#       - "prop_pastures_within_7071m" (correlated with distance.to.pastures)
X_VARIABLES <- c(
    "chelsa_hurs", "chelsa_pr", "elevation", "ndvi", "light_pollution", 

    "distance.to.artificial.surfaces", "distance.to.arable.land", 
    "distance.to.pastures", "distance.to.heterogeneous.agricultural.areas", 
    "distance.to.wetlands", "distance.to.forest.and.semi.natural.areas", 
    "distance.to.water.bodies",

    "prop_heterogeneous.agricultural.areas_within_7071m", 
    "prop_wetlands_within_7071m",
    "prop_forest.and.semi.natural.areas_within_7071m",
    
    "swf")
X_UNITS <- list(
    "%", expression(paste("kg m"^{-2}, " month"^{-1})), "m", "Index", "DN",
    "m", "m", "m", "m", "m", "m", "m",
    "%", "%", "%",
    "%"
)

# Parameters of interest for Ds-optimality 
# (other parameters are considered as noise)
ACTIONABLE_VARIABLES <- c("light_pollution") 


NAME_SP_SIMUL <- paste0("sp_", seq(1:length(NAMES_SPECIES)))

                       
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

GG_TERRAIN_PALETTE <- colorRampPalette(c(
    "#f7f4f9", "#f7f4f9", "#e8a0a0", "#e8935a", "#e8c93a", 
    "#8fc93a", "#2ca02c", "#22b0b0", "#22b0b0"))(10)