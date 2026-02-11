#' @title Apply Random Forest Classification
#' @description Apply trained random forest model to reference data and assess certainty.
#' @author Megan Coffer
#'
#' @details
#' This script is part of the Chesapeake Bay Marsh Mapping workflow.
#' Run this script third in the Model Development workflow.
#'
#' Processing steps:
#' 1. Apply trained random forest model to formatted imagery
#' 2. Generate class probability predictions
#' 3. Apply certainty thresholding to flag low-certainty predictions
#'
#' @inputs
#' - Formatted Planet imagery rasters
#' - Formatted DEM rasters
#' - Training data shapefiles
#' - Trained Random Forest models (.RData)
#' - Data preprocessing objects (.RData)
#'
#' @outputs
#' - Classified rasters with predictions and associated probabilities
#' - Classification key text files


# ==============================================================================
# SETUP
# ==============================================================================

# Load configuration - check current directory first, then relative paths
if (file.exists("config_local.R")) {
  source("config_local.R")
} else if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists(file.path(dirname(getwd()), "config_local.R"))) {
  source(file.path(dirname(getwd()), "config_local.R"))
} else {
  source(file.path(dirname(getwd()), "config.R"))
}

# Load required packages
library(terra)
library(sf)
library(tidyterra)
library(tidymodels)

# Set working directory
main.dir <- MODEL_DEV_DIR

# Set certainty level for random forest probabilities
prob.level <- CERTAINTY_THRESHOLD


# ==============================================================================
# READ INPUT DATA
# ==============================================================================

# Read in formatted satellite data
sr.formatted <- lapply(
  list.files(
    file.path(main.dir, "Data", "Input_Data", "Planet_Multispectral_Data", "3_Level3B_Data_Formatted"),
    "*.tif$",
    full.names = TRUE
  ),
  rast
)

dem.formatted <- lapply(
  list.files(
    file.path(main.dir, "Data", "Input_Data", "USGS_DEM_Data"),
    "*.tif$",
    full.names = TRUE
  ),
  rast
)

# Read in training data coordinates
training.data.coords <- lapply(
  list.files(
    file.path(main.dir, "Data", "Output_Data", "Training_Data", "Training_Data_Shapefile"),
    "*_Training_Data.shp$",
    full.names = TRUE
  ),
  st_read
)
training.data.coords <- lapply(training.data.coords, function(x) { st_transform(x, crs(sr.formatted[[1]])) })

# Read in training data filenames and model files
training.data.filenames <- list.files(
  file.path(main.dir, "Data", "Output_Data", "Training_Data"),
  "*_Training_Data.csv$",
  full.names = TRUE
)

data.preps <- list.files(
  file.path(main.dir, "Data", "Output_Data", "Random_Forest_Model"),
  "*_data_prep.RData$",
  full.names = TRUE
)

rf.models <- list.files(
  file.path(main.dir, "Data", "Output_Data", "Random_Forest_Model"),
  "*_random_forest_model.RData$",
  full.names = TRUE
)


# ==============================================================================
# APPLY RANDOM FOREST CLASSIFICATION
# ==============================================================================

# Loop through each image and apply classification
for (image in seq_along(sr.formatted)) {

  # Subset to the image of interest and add layer for training data
  training.data.rast <- rasterize(training.data.coords[[image]], sr.formatted[[image]], background = 0)
  names(training.data.rast) <- "training_data"
  formatted.image <- c(sr.formatted[[image]], crop(dem.formatted[[image]], ext(sr.formatted[[image]])), training.data.rast)

  # Read in corresponding random forest data
  data.preps.image <- readRDS(data.preps[[image]])
  rf.models.image <- readRDS(rf.models[[image]])

  # Extract coordinates from image
  sr.coords <- xyFromCell(formatted.image, 1:ncell(formatted.image))
  formatted.image$x <- sr.coords[, 1]
  formatted.image$y <- sr.coords[, 2]

  # Convert to a tibble and compute vegetation indices
  tibble.image <- tidyterra::as_tibble(formatted.image) %>% filter(complete.cases(.))
  tibble.image$NDVI  <- (tibble.image$nir - tibble.image$red) / (tibble.image$nir + tibble.image$red)
  tibble.image$EVI   <- 2.5 * (tibble.image$nir - tibble.image$red) / ((tibble.image$nir + 6 * tibble.image$red - 7.5 * tibble.image$blue) + 1)
  tibble.image$NDWI  <- (tibble.image$green - tibble.image$nir) / (tibble.image$green + tibble.image$nir)
  tibble.image$SAVI  <- ((tibble.image$nir - tibble.image$red) * 1.5) / (tibble.image$nir + tibble.image$red + 0.5)
  tibble.image$WDRVI <- (0.2 * tibble.image$nir - tibble.image$red) / (0.2 * tibble.image$nir + tibble.image$red)
  tibble.image$BG    <- tibble.image$blue / tibble.image$green
  tibble.image$GR    <- tibble.image$green / tibble.image$red
  tibble.image$NIRR  <- tibble.image$nir / tibble.image$red
  tibble.image$CLASS <- tibble.image$VEG_CLASS
  tibble.image$VEG_CLASS <- NULL

  # Apply random forest model
  prep.data.image <- data.preps.image %>% bake(tibble.image)
  predict.data.image <- rf.models.image %>% predict(prep.data.image) %>% bind_cols(prep.data.image)

  # Compute probabilities for each class and append to random forest results
  predict.data.image <- predict(rf.models.image, prep.data.image, type = "prob") %>% bind_cols(predict.data.image)

  # Create certainty-based columns
  predict.data.image$MAX_PROB <- pmax(
    predict.data.image$.pred_Juncus_roemerianus,
    predict.data.image$.pred_NoVeg,
    predict.data.image$.pred_Phragmites_australis,
    predict.data.image$.pred_Spartina_alterniflora,
    predict.data.image$.pred_Spartina_patens,
    na.rm = TRUE
  )

  # Add band based on given certainty threshold
  predict.data.image$PRED_CLASS_PROB <- ifelse(
    predict.data.image$MAX_PROB >= prob.level,
    as.character(predict.data.image$.pred_class),
    "Low_certainty"
  )

  # Convert classification to raster
  predict.data.image$x <- tibble.image$x
  predict.data.image$y <- tibble.image$y
  classified.image <- as_spatraster(
    predict.data.image,
    xycols = c(which(colnames(predict.data.image) == "x"), which(colnames(predict.data.image) == "y")),
    crs = crs(sr.formatted[[image]])
  )

  # Export raster
  writeRaster(
    classified.image,
    file.path(
      main.dir, "Data", "Output_Data", "Classified_Data",
      paste0(substr(basename(sources(sr.formatted[[image]])), 1, 8), "_classification_reference_data.tif")
    ),
    overwrite = TRUE, 
    datatype = "FLT4S"  # Explicitly preserve float precision
  )

  # Export classification key
  sink(file = file.path(
    main.dir, "Data", "Output_Data", "Classified_Data",
    paste0(substr(basename(sources(sr.formatted[[image]])), 1, 8), "_classification_key.txt")
  ))
  cat("Classification key for PRED_CLASS_PROB:\n")
  print(unique(classified.image[["PRED_CLASS_PROB"]]))
  sink(file = NULL)

}
