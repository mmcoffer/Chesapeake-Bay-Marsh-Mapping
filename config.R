#' Configuration File for Chesapeake Bay Marsh Mapping
#'
#' INSTRUCTIONS:
#' 1. Copy this file to config_local.R
#' 2. Edit the BASE_DIR path in config_local.R to match your system
#' 3. Do not commit config_local.R to version control
#'
#' The scripts will automatically use config_local.R if it exists,
#' otherwise they will use this file (config.R).


# ==============================================================================
# BASE DIRECTORY - EDIT THIS PATH
# ==============================================================================

# Set this to the root directory of the project on your system
BASE_DIR <- file.path(Sys.getenv("HOME"), "Chesapeake-Bay-Marsh-Mapping")


# ==============================================================================
# DERIVED PATHS (typically no need to edit)
# ==============================================================================

# Model Development directories
MODEL_DEV_DIR <- file.path(BASE_DIR, "1_Model_Development")
MODEL_DEV_INPUT_DIR <- file.path(MODEL_DEV_DIR, "Data", "Input_Data")
MODEL_DEV_OUTPUT_DIR <- file.path(MODEL_DEV_DIR, "Data", "Output_Data")

# Annual Assessments directories
ANNUAL_DIR <- file.path(BASE_DIR, "2_Annual_Assessments")
ANNUAL_INPUT_DIR <- file.path(ANNUAL_DIR, "Data", "Input_Data")
ANNUAL_OUTPUT_DIR <- file.path(ANNUAL_DIR, "Data", "Output_Data")


# ==============================================================================
# MODEL PARAMETERS
# ==============================================================================

# Certainty threshold for classification (predictions below this are flagged)
CERTAINTY_THRESHOLD <- 0.7

# Fraction of smallest class to use for training data
TRAINING_FRACTION <- 0.8

# Correlation threshold for removing highly correlated predictors
CORRELATION_THRESHOLD <- 0.95


# ==============================================================================
# VEGETATION CLASSES
# ==============================================================================

VEG_CLASSES <- c(
  "Spartina_alterniflora",
  "Spartina_patens",
  "Phragmites_australis",
  "Juncus_roemerianus",
  "NoVeg"
)


# ==============================================================================
# HELPER FUNCTION
# ==============================================================================

#' Validate that required directories exist
#' @param dirs Character vector of directory paths to check
#' @param create Logical; if TRUE, create missing directories
validate_dirs <- function(dirs, create = FALSE) {
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      if (create) {
        dir.create(dir, recursive = TRUE)
        message("Created directory: ", dir)
      } else {
        warning("Directory does not exist: ", dir)
      }
    }
  }
}


# ==============================================================================
# CREATE OUTPUT DIRECTORIES
# ==============================================================================

# Create required output directories if they don't exist
required_dirs <- c(
  # Model Development Data directories
  file.path(MODEL_DEV_INPUT_DIR, "Planet_Multispectral_Data", "2_Level_3B_Data_Unzipped"),
  file.path(MODEL_DEV_INPUT_DIR, "Planet_Multispectral_Data", "3_Level3B_Data_Formatted"),
  file.path(MODEL_DEV_INPUT_DIR, "USGS_DEM_Data"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Training_Data"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Training_Data", "Training_Data_Shapefile"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Random_Forest_Model"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Variable_Importance"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Agreement_Statistics"),
  file.path(MODEL_DEV_OUTPUT_DIR, "Classified_Data"),
  # Annual Assessments directories
  file.path(ANNUAL_INPUT_DIR, "Planet_Multispectral_Data", "2_Level_3B_Data_Unzipped"),
  file.path(ANNUAL_INPUT_DIR, "Planet_Multispectral_Data", "3_Level3B_Data_Formatted"),
  file.path(ANNUAL_INPUT_DIR, "USGS_DEM_Data"),
  file.path(ANNUAL_OUTPUT_DIR),
  file.path(ANNUAL_OUTPUT_DIR, "Classified_Data")
)

validate_dirs(required_dirs, create = TRUE)
