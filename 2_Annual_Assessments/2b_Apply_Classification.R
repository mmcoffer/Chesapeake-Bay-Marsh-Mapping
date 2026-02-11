#' @title Apply Classification for Annual Assessments
#' @description Apply trained Random Forest model to satellite imagery for annual monitoring.
#' @author Megan Coffer
#'
#' @details
#' This script is part of the Chesapeake Bay Marsh Mapping workflow.
#' Run this script second in the Annual Assessments workflow.
#'
#' Processing steps:
#' 1. Load trained Random Forest model
#' 2. Apply model to formatted annual imagery
#' 3. Generate class probability predictions
#' 4. Apply certainty thresholding
#' 5. Apply 3x3 modal smoothing filter
#'
#' @inputs
#' - Formatted Planet imagery rasters
#' - Formatted DEM rasters
#' - Trained Random Forest model (.RData)
#' - Data preprocessing object (.RData)
#'
#' @outputs
#' - Classified rasters with predictions and smoothed output (per image/date)


# ==============================================================================
# SETUP
# ==============================================================================

# Load configuration - check current directory first, then relative paths
if (file.exists("config_local.R")) {
  source("config_local.R")
} else if (file.exists("config.R")) {
  source("config.R")
} else if (file.exists(file.path(dirname(dirname(getwd())), "config_local.R"))) {
  source(file.path(dirname(dirname(getwd())), "config_local.R"))
} else {
  source(file.path(dirname(dirname(getwd())), "config.R"))
}

# Load required packages
library(terra)
library(sf)
library(tidymodels)
library(tidyterra)

# Set working directory
main.dir <- BASE_DIR


# ==============================================================================
# READ INPUT DATA
# ==============================================================================

# List formatted satellite data
sr.filenames <- list.files(
  file.path(main.dir, "2_Annual_Assessments", "Data", "Input_Data", "Planet_Multispectral_Data", "3_Level3B_Data_Formatted"),
  "*.tif$",
  full.names = TRUE
)

dem.filenames <- list.files(
  file.path(main.dir, "2_Annual_Assessments", "Data", "Input_Data", "USGS_DEM_Data"),
  pattern = "processed.*\\.tif$",
  full.names = TRUE
)


# ==============================================================================
# SELECT BEST MODEL
# ==============================================================================

# Find all agreement statistics files
agreement.dir <- file.path(main.dir, "1_Model_Development", "Data", "Output_Data", "Agreement_Statistics")
agreement.files <- list.files(agreement.dir, "*_agreement_reference_data.csv$", full.names = TRUE)

# Read each file and extract overall balanced agreement
agreement.results <- lapply(agreement.files, function(f) {
  agreement.data <- read.csv(f)
  # Extract date from filename (first 8 characters)
  file.date <- substr(basename(f), 1, 8)
  # Get overall balanced agreement (same value in all rows)
  overall.bal.agreement <- agreement.data$OVERALL_BAL_AGREEMENT[1]
  data.frame(date = file.date, file = f, overall_bal_agreement = overall.bal.agreement)
})

# Combine results and find best date
agreement.summary <- do.call(rbind, agreement.results)
best.date <- agreement.summary$date[which.max(agreement.summary$overall_bal_agreement)]

message("Model selection based on overall balanced agreement:")
message(paste("  Best date:", best.date, "with agreement:", max(agreement.summary$overall_bal_agreement)))

# Read in trained Random Forest model for best date
data.prep <- readRDS(
  list.files(
    file.path(main.dir, "1_Model_Development", "Data", "Output_Data", "Random_Forest_Model"),
    paste0(best.date, "_data_prep.RData$"),
    full.names = TRUE
  )
)

rf.model <- readRDS(
  list.files(
    file.path(main.dir, "1_Model_Development", "Data", "Output_Data", "Random_Forest_Model"),
    paste0(best.date, "_random_forest_model.RData$"),
    full.names = TRUE
  )
)


# ==============================================================================
# APPLY CLASSIFICATION
# ==============================================================================

# Loop through each unique date and apply classification
for (unique.date in seq_along(sr.filenames)) {
  # List files for the given date
  sr.files <- list.files(sr.filenames[[unique.date]], "*.tif$", full.names = TRUE)
  dem.file <- rast(dem.filenames[[unique.date]])

  # Loop through each file and process
  for (sr.file in seq_along(sr.files)) {
    # Read in the raster of interest and append DEM
    sr.rast <- rast(sr.files[[sr.file]])
    formatted.image <- c(sr.rast, crop(dem.file, ext(sr.rast)))

    # Extract coordinates from image
    sr.coords <- xyFromCell(formatted.image, 1:ncell(sr.rast))
    formatted.image$x <- sr.coords[, 1]
    formatted.image$y <- sr.coords[, 2]

    # Convert to a tibble and compute vegetation indices
    tibble.image <- tidyterra::as_tibble(formatted.image) %>%
      filter(complete.cases(.)) %>%
      filter_all(all_vars(!is.infinite(.)))

    tibble.image$NDVI <- (tibble.image$nir - tibble.image$red) / (tibble.image$nir + tibble.image$red)
    tibble.image$EVI <- 2.5 * (tibble.image$nir - tibble.image$red) / ((tibble.image$nir + 6 * tibble.image$red - 7.5 * tibble.image$blue) + 1)
    tibble.image$NDWI <- (tibble.image$green - tibble.image$nir) / (tibble.image$green + tibble.image$nir)
    tibble.image$SAVI <- ((tibble.image$nir - tibble.image$red) * 1.5) / (tibble.image$nir + tibble.image$red + 0.5)
    tibble.image$WDRVI <- (0.2 * tibble.image$nir - tibble.image$red) / (0.2 * tibble.image$nir + tibble.image$red)
    tibble.image$BG <- tibble.image$blue / tibble.image$green
    tibble.image$GR <- tibble.image$green / tibble.image$red
    tibble.image$NIRR <- tibble.image$nir / tibble.image$red
    tibble.image$VEG_CLASS <- NULL

    # Apply random forest model
    prep.data.image <- data.prep %>%
      bake(tibble.image) %>%
      filter(complete.cases(.)) %>%
      filter_all(all_vars(!is.infinite(.)))

    predict.data.image <- rf.model %>%
      predict(prep.data.image) %>%
      bind_cols(prep.data.image)

    # Compute probabilities for each class and append to random forest results
    predict.data.image <- predict(rf.model, prep.data.image, type = "prob") %>%
      bind_cols(predict.data.image)

    # Create probability-based columns
    cert.data.image <- predict.data.image %>%
      mutate(MAX_PROB = pmax(
        .pred_Juncus_roemerianus,
        .pred_NoVeg,
        .pred_Phragmites_australis,
        .pred_Spartina_alterniflora,
        .pred_Spartina_patens,
        na.rm = TRUE
      ))

    pred.data.image$.pred_class_certainty <- ifelse(
      pred.data.image$MAX_PROB >= 0.5,
      levels(pred.data.image$.pred_class)[pred.data.image$.pred_class],
      "Low certainty"
    )

    # Convert classification to raster
    predict.data.image$x <- tibble.image$x
    predict.data.image$y <- tibble.image$y
    classified.image <- as_spatraster(
      predict.data.image,
      xycols = c(which(colnames(predict.data.image) == "x"), which(colnames(predict.data.image) == "y")),
      crs = crs(sr.rast)
    )

    # Apply smoothing via mode with 3x3 moving window
    smoothed.image <- terra::focal(classified.image[[".pred_class"]], w = matrix(1, nc = 3, nr = 3), fun = "modal")
    levels(smoothed.image) <- data.frame(id = 1:5, .pred_class_smooth = unique(classified.image$.pred_class)$.pred_class)
    stacked.image <- c(classified.image, smoothed.image)

    # Export rasters -- create folder for date if it doesn't already exist
    output.dir <- file.path(
      main.dir,
      "2_Annual_Assessments", "Data", "Output_Data", "Classified_Data_Per_Date",
      basename(sr.directories)[unique.date]
    )

    if (!dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE)
    }

    writeRaster(
      stacked.image,
      file.path(
        output.dir,
        paste0(basename(sr.directories)[unique.date], "_classified_image_", sr.file, ".tif")
      ),
      overwrite = TRUE
    )
  }
}
