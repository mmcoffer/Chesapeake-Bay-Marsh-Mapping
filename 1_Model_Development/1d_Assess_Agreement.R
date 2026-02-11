#' @title Assess Classification Agreement
#' @description Validate classification accuracy against reference data.
#' @author Megan Coffer
#'
#' @details
#' This script is part of the Chesapeake Bay Marsh Mapping workflow.
#' Run this script last in the Model Development workflow.
#'
#' Processing steps:
#' 1. Compare classified images against reference vegetation data
#' 2. Compute confusion matrices and agreement statistics
#' 3. Calculate per-class and overall accuracy metrics
#' 4. Assess regional variation in classification performance
#'
#' @inputs
#' - Classified rasters from Random Forest model
#' - VIMS reference vegetation shapefile
#' - Marsh region boundary shapefile
#'
#' @outputs
#' - Overall agreement statistics CSV (per image)
#' - Regional agreement statistics CSV - main (per image)
#' - Regional agreement statistics CSV - supplemental (per image)


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


# ==============================================================================
# READ INPUT DATA
# ==============================================================================

# Read in classified satellite data
classified.images <- lapply(
  list.files(
    file.path(main.dir, "Data", "Output_Data", "Classified_Data"),
    "*.tif$",
    full.names = TRUE
  ),
  rast
)

# Read in reference data
ref <- st_read(file.path(main.dir, "Data", "Input_Data", "VIMS_Reference_Data", "Vegetation_shapefile_May_2021.shp"))


# ==============================================================================
# ASSESS AGREEMENT
# ==============================================================================

# Loop through each image and assess agreement
for (image in 1:length(classified.images)) {
  # Subset to the image of interest
  classified.image <- classified.images[[image]]

  # Format data for assessing agreement
  classified.image$PRED_CLASS_PROB[classified.image$PRED_CLASS_PROB == "Low_certainty"] <- NA
  agreement.tibble <- tidyterra::as_tibble(classified.image) %>%
    filter(if_all(c(CLASS, PRED_CLASS_PROB), complete.cases))

  agreement.tibble$PRED_CLASS_PROB <- factor(
    agreement.tibble$PRED_CLASS_PROB,
    levels = c("Spartina_alterniflora", "Phragmites_australis", "Juncus_roemerianus", "Spartina_patens", "NoVeg")
  )
  agreement.tibble$CLASS <- factor(
    agreement.tibble$CLASS,
    levels = c("Spartina_alterniflora", "Phragmites_australis", "Juncus_roemerianus", "Spartina_patens", "NoVeg")
  )

  multiclass.metrics <- metric_set(sensitivity, specificity, bal_accuracy, accuracy)

  # ============================================================================
  # Agreement across all reference data
  # ============================================================================

  # Create dataframe to populate with results
  agreement.mat <- conf_mat(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_PROB)
  agreement.df <- as.data.frame(matrix(nrow = 5, ncol = 5, as.data.frame(tidy(agreement.mat))$value))
  agreement.df <- cbind(names(agreement.mat[[1]][, 1]), agreement.df)
  colnames(agreement.df) <- c("PREDICTED", names(agreement.mat[[1]][, 1]))

  # Class-by-class agreement statistics
  per.class.metrics <- agreement.tibble %>%
    group_by(CLASS) %>%
    multiclass.metrics(truth = CLASS, estimate = PRED_CLASS_PROB, estimator = "micro")

  agreement.df$CLASS_SENSITIVITY <- round(subset(per.class.metrics, .metric == "sensitivity")$.estimate * 100)
  agreement.df$CLASS_SPECIFICITY <- round(subset(per.class.metrics, .metric == "specificity")$.estimate * 100)
  agreement.df$CLASS_BAL_AGREEMENT <- round(subset(per.class.metrics, .metric == "bal_accuracy")$.estimate * 100)
  agreement.df$CLASS_AGREEMENT <- round(subset(per.class.metrics, .metric == "accuracy")$.estimate * 100)

  # Overall agreement statistics
  agreement.df$OVERALL_SENSITIVITY <- rep(
    round(sensitivity(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_PROB, estimator = "micro")$.estimate * 100),
    each = 5
  )
  agreement.df$OVERALL_SPECIFICITY <- rep(
    round(specificity(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_PROB, estimator = "micro")$.estimate * 100),
    each = 5
  )
  agreement.df$OVERALL_BAL_AGREEMENT <- rep(
    round(bal_accuracy(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_PROB, estimator = "micro")$.estimate * 100),
    each = 5
  )
  agreement.df$OVERALL_AGREEMENT <- rep(
    round(accuracy(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_PROB, estimator = "micro")$.estimate * 100),
    each = 5
  )

  # Export as CSV file
  write.csv(
    agreement.df,
    file.path(
      main.dir, "Data", "Output_Data", "Agreement_Statistics",
      paste0(substr(basename(sources(classified.images[[image]])), 1, 8), "_agreement_reference_data.csv")
    ),
    row.names = FALSE
  )
  
}
