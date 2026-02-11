#' Install Required Packages for Chesapeake Bay Marsh Mapping
#'
#' Run this script once before executing the workflow to ensure
#' all required packages are installed.
#'
#' Usage:
#'   source("requirements.R")


# ==============================================================================
# REQUIRED PACKAGES
# ==============================================================================

packages <- c(
  "terra",        # Raster data processing
  "sf",           # Vector data processing
  "tidymodels",   # Machine learning framework (includes recipes, parsnip, etc.)
  "tidyterra",    # Integration of terra with tidyverse
  "FedData",      # Download USGS elevation data
  "randomForest", # Random Forest implementation
  "tools"         # File path utilities (included with base R)
)


# ==============================================================================
# INSTALLATION FUNCTION
# ==============================================================================

install_if_missing <- function(pkg) {
  # Use require() to actually load package (checks all dependencies)
  if (!suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))) {
    message("Installing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
    # Verify it loads after installation
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop("Failed to install package: ", pkg)
    }
  } else {
    message("Package already installed: ", pkg)
  }
}


# ==============================================================================
# INSTALL PACKAGES
# ==============================================================================

message("Checking and installing required packages...\n")

invisible(lapply(packages, install_if_missing))

message("\n")
message("========================================")
message("All required packages are installed.")
message("========================================")


# ==============================================================================
# VERIFY INSTALLATION
# ==============================================================================

message("\nVerifying package versions:\n")

for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    message(sprintf("  %-15s %s", pkg, version))
  }
}

message("\nSetup complete. You can now run the analysis scripts.")
