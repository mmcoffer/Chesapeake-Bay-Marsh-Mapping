#' @title Format Satellite Data
#' @description Format Planet multispectral imagery and USGS DEM data for model development.
#' @author Megan Coffer
#'
#' @details
#' This script is part of the Chesapeake Bay Marsh Mapping workflow.
#' Run this script first in the Model Development workflow.
#'
#' Processing steps:
#' 1. Unzip Planet Level 3B multispectral imagery
#' 2. Apply cloud masking using UDM2 quality band included in Planet Level 3B multispectral imagery 
#' 3. Merge overlapping images from the same date
#' 4. Download and resample USGS DEM to match imagery
#'
#' @inputs
#' - Planet Level 3B zipped imagery files
#' - VIMS reference vegetation shapefile
#'
#' @outputs
#' - Formatted Planet imagery rasters (per date)
#' - Resampled DEM rasters (per date)


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
require(terra)
require(sf)
require(FedData)
require(tools)


# Set working directory
main.dir <- MODEL_DEV_INPUT_DIR


# ==============================================================================
# READ INPUT DATA
# ==============================================================================

# Planet imagery - unzip if needed
zipped.list <- list.files(
  file.path(main.dir, "Planet_Multispectral_Data", "1_Level_3B_Data_Zipped"),
  "*.zip",
  full.names = TRUE
)

zipped.path <- lapply(zipped.list, function(x) {
  file.path(
    main.dir,
    "Planet_Multispectral_Data",
    "2_Level_3B_Data_Unzipped",
    substr(basename(x), 1, nchar(basename(x)) - 4)
  )
})

for (zipped.file in seq_along(zipped.list)) {
  if (!file.exists(zipped.path[[zipped.file]])) {
    unzip(zipped.list[[zipped.file]], exdir = zipped.path[[zipped.file]])
  }
}

# List unzipped data
unzipped.list <- list.files(
  file.path(main.dir, "Planet_Multispectral_Data", "2_Level_3B_Data_Unzipped"),
  full.names = TRUE
)

# Reference data
ref <- st_read(file.path(main.dir, "VIMS_Reference_Data", "Vegetation_shapefile_May_2021.shp"))

# DEM data
feddata.get <- get_ned(
  template = ref,
  label = "VIMS_CCRM",
  res = 13,
  extraction.dir = file.path(main.dir, "USGS_DEM_Data", "Extracted_DEM_Raster"),
  force.redo = TRUE
)
feddata.dem <- rast(file.path(main.dir, "USGS_DEM_Data", "Extracted_DEM_Raster", "VIMS_CCRM_NED_13.tif"))


# ==============================================================================
# 1. PROCESS PLANET MULTISPECTRAL DATA
# ==============================================================================

# List all surface reflectance files
sr.filenames <- list.files(
  file.path(unzipped.list, "PSScene"),
  "*SR_8b.tif$",
  full.names = TRUE
)

sr.masked <- vector(mode = "list", length = length(sr.filenames))

# Loop through each Planet file and process
for (image in seq_along(sr.filenames)) {

  # Read in raster, corresponding UDM raster, and remove invalid data
  sr.image <- rast(sr.filenames[[image]])
  udm2.image <- rast(paste0(
    dirname(sr.filenames[[image]]), "/",
    substr(basename(sr.filenames[[image]]), 1, 27),
    "udm2.tif"
  ))
  sr.image[udm2.image$clear == 0] <- NA

  # Mask to corresponding reference data
  sr.masked[[image]] <- mask(sr.image, ref)
  names(sr.masked)[[image]] <- file_path_sans_ext(basename(sr.filenames[[image]]))
}

# Find unique dates for Planet data
unique.dates <- unique(substr(basename(sr.filenames), 1, 8))
sr.formatted <- vector(mode = "list", length = length(unique.dates))

# Combine Planet data from the same dates
for (unique.date in seq_along(unique.dates)) {

  # Take the mean of overlapping pixels for images collected on the same day
  i.dates <- sr.masked[which(substr(names(sr.masked), 1, 8) == unique.dates[unique.date])]
  unique.date.merge <- merge(sprc(i.dates))
  sr.formatted[[unique.date]] <- crop(unique.date.merge, ext(ref))
  names(sr.formatted)[[unique.date]] <- unique.dates[unique.date]

  # Rasterize reference data and add layer for corresponding reference data
  ref.rast <- rasterize(ref, sr.formatted[[unique.date]], "VEG_CLASS")
  sr.formatted[[unique.date]] <- c(sr.formatted[[unique.date]], ref.rast)

  # Export raster
  writeRaster(
    sr.formatted[[unique.date]],
    file.path(
      main.dir,
      "Planet_Multispectral_Data",
      "3_Level3B_Data_Formatted",
      paste0(unique.dates[[unique.date]], "_3B_AnalyticMS_SR_8b_clip_processed.tif")
    ),
    overwrite = TRUE
  )
}


# ==============================================================================
# 2. PROCESS DEM DATA
# ==============================================================================

# Create corresponding DEM dataset for each unique date
for (unique.date in seq_along(unique.dates)) {

  # Project and resample to match Planet data
  feddata.dem.proj <- project(feddata.dem, crs(sr.formatted[[unique.date]]))
  feddata.dem.resample <- resample(feddata.dem.proj, sr.formatted[[unique.date]], "bilinear")

  # Crop and mask to reference data
  feddata.dem.ref <- crop(feddata.dem.resample, ref, mask = TRUE)

  # Rename raster data
  names(feddata.dem.ref) <- "DEM"

  # Export processed raster
  writeRaster(
    feddata.dem.ref,
    file.path(main.dir, "USGS_DEM_Data", paste0("USGS_DEM_Data_processed_", unique.dates[[unique.date]], ".tif")),
    overwrite = TRUE
  )
}
