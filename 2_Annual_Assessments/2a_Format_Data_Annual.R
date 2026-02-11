#' @title Format Data for Annual Assessments within NOAA's Middle Peninsula Habitat Focus Area (HFA)
#' @description Define the Middle Peninsula study region and format Planet multispectral imagery and DEM data for annual assessments.
#' @author Megan Coffer
#'
#' @details
#' This script is part of the Chesapeake Bay Marsh Mapping workflow.
#' Run this script first in the Annual Assessments workflow (but after the entire model development workflow)
#'
#' Processing steps:
#' 1. Define Middle Peninsula study region using TMI, HFA, and salinity data
#' 3. Unzip Planet Level 3B multispectral imagery
#' 4. Apply cloud masking using UDM2 quality band
#' 5. Subset to NOAA's Habitat Focus Area (HFA) with salinity constraint
#' 6. Merge overlapping images from the same date
#' 7. Resample USGS DEM to match imagery
#'
#' @inputs
#' - VIMS CCRM TMI data shapefile
#' - Middle Peninsula HFA boundary shapefile
#' - Chesapeake Bay mean surface salinity shapefile
#' - Planet Level 3B zipped imagery files
#' - UDM2 cloud mask files
#'
#' @outputs
#' - Simplified HFA boundary shapefile
#' - Formatted salinity shapefile
#' - Formatted Planet imagery rasters (per date)
#' - Resampled DEM rasters (per date)
#' - PlanetScope image identifier lists (per date)


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
library(alphahull)
library(FedData)
library(tools)

# Set working directories
main.dir <- ANNUAL_INPUT_DIR


# ==============================================================================
# DEFINE MIDDLE PENINSULA STUDY AREA
# ==============================================================================

# Read in VIMS CCRM TMI data, HFA boundary, and salinity data; reproject to WGS84
tmi <- st_transform(st_read(file.path(main.dir, "VIMS_Tidal_Marsh_Inventory","VA_TMI_2011_2019_utm18.shp")), crs = "EPSG:4326")
hfa <- st_transform(st_read(file.path(main.dir, "Middle_Peninsula_HFA", "Middle_Peninsula_HFA_Boundary.shp")), crs = "EPSG:4326")
salinity <- st_transform(st_read(file.path(main.dir, "Mean_Surface_Salinity", "Chesapeake_Bay_Mean_Surface_Salinity_1985-2018.shp")), crs = "EPSG:4326")

# Crop TMI by HFA
tmi.hfa.crop <- st_intersection(st_make_valid(tmi), st_make_valid(st_geometry(hfa)))

# Simplify geometry of tmi.hfa.crop object for cropping classification datasets
# Retrieve centroid of each TMI polygon
tmi.hfa.centroid <- st_centroid(tmi.hfa.crop)

# Apply alpha hull on centroids to create boundary of points
unique.centroids <- unique(round(st_coordinates(tmi.hfa.centroid), digits = 4))
tmi.hfa.ahull <- ahull(unique.centroids, alpha = 0.03)

# Convert alpha hull object to sf object
tmi.hfa.ahull.edges <- data.frame(tmi.hfa.ahull$ashape.obj$edges)[, c("x1", "y1", "x2", "y2")]
tmi.hfa.ahull.linestring <- st_linestring(matrix(as.numeric(tmi.hfa.ahull.edges[1, ]), ncol = 2, byrow = TRUE))
for (tmi.hfa.ahull.edge in 2:nrow(tmi.hfa.ahull.edges)) {
  tmi.hfa.ahull.linestring <- c(
    tmi.hfa.ahull.linestring,
    st_linestring(matrix(as.numeric(tmi.hfa.ahull.edges[tmi.hfa.ahull.edge, ]), ncol = 2, byrow = TRUE))
  )
}
tmi.hfa.ahull.sf <- st_sf(geom = st_sfc(tmi.hfa.ahull.linestring), crs = 4326) %>%
  st_polygonize() %>%
  st_collection_extract()

# Dissolve to remove interior polygons
tmi.hfa.union <- tmi.hfa.ahull.sf %>%
  st_make_valid() %>%
  st_union()

# Buffer by 1.25 km to ensure all polygons are covered and simplify
tmi.hfa.buffer <- st_buffer(tmi.hfa.union, dist = 1250)
tmi.hfa.simplify <- st_simplify(tmi.hfa.buffer, preserveTopology = TRUE, dTolerance = 1000)

# Subset to salinities greater than 15 ppt and intersect with HFA
salinity.subset <- subset(salinity, SALINITY > 15)
mid.pen.study.area <- st_intersection(st_union(salinity.subset), tmi.hfa.simplify)

message("TMI successfully cropped to HFA with salinities above 15 ppt")

# ==============================================================================
# READ CLASSIFICATION INPUT DATA
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

if (length(zipped.list) > 0) {
  for (zipped.file in seq_along(zipped.list)) {
    if (!file.exists(zipped.path[[zipped.file]])) {
      unzip(zipped.list[[zipped.file]], exdir = zipped.path[[zipped.file]])
    }
  }
}

unzipped.list <- list.files(
  file.path(main.dir, "Planet_Multispectral_Data", "2_Level_3B_Data_Unzipped"),
  full.names = TRUE
)

# DEM data
feddata.get <- get_ned(
  template = mid.pen.study.area,
  label = "Mid_Pen_Study_Area",
  res = 13,
  extraction.dir = file.path(main.dir, "USGS_DEM_Data"),
  force.redo = FALSE
)
feddata.dem <- rast(file.path(main.dir, "USGS_DEM_Data", "Mid_Pen_Study_Area_NED_13.tif"))


# ==============================================================================
# PROCESS PLANET MULTISPECTRAL DATA
# ==============================================================================

# List all surface reflectance files
sr.filenames <- list.files(
  file.path(unzipped.list, "PSScene"),
  "*SR_8b.tif$",
  full.names = TRUE
)

unique.dates <- unique(substr(basename(sr.filenames), 1, 8))

# Loop through each unique date and process
for (unique.date in unique.dates) {
  # Subset to the date of interest
  sr.filenames.date <- sr.filenames[which(substr(basename(sr.filenames), 1, 8) == unique.date)]
  unique.date.list <- vector(mode = "list", length = length(sr.filenames.date))

  # Loop through each PlanetScope raster and process
  for (sr.filename.date in seq_along(sr.filenames.date)) {
    # Read in raster
    sr.image <- rast(sr.filenames.date[[sr.filename.date]])

    # Check if it intersects mid.pen.study.area
    if (!is.null(intersect(ext(sr.image), st_bbox(mid.pen.study.area)))) {
      # Read in corresponding UDM raster and remove invalid data
      udm2.image <- rast(paste0(
        dirname(sr.filenames.date[[sr.filename.date]]), "/",
        substr(basename(sr.filenames.date[[sr.filename.date]]), 1, 27),
        "udm2.tif"
      ))
      sr.image[udm2.image$clear == 0] <- NA

      # Mask to corresponding reference data
      mid.pen.proj <- project(vect(mid.pen.study.area), crs(sr.image))
      unique.date.list[[sr.filename.date]] <- mask(sr.image, mid.pen.proj)
      names(unique.date.list)[[sr.filename.date]] <- file_path_sans_ext(basename(sr.filenames.date[[sr.filename.date]]))

      # Progress update
      message(paste("Processing raster", sr.filename.date, "of", length(sr.filenames.date), "within region of interest"))
    }
  }

  # Aggregate overlapping pixels for images collected on the same day
  unique.date.list.nonnull <- Filter(Negate(is.null), unique.date.list)
  unique.date.merge <- terra::merge(sprc(unique.date.list.nonnull))
  unique.date.crop <- crop(unique.date.merge, mid.pen.proj, mask = TRUE)

  # Create corresponding DEM dataset for each unique date
  feddata.dem.proj <- project(feddata.dem, crs(unique.date.crop), mask = TRUE, threads = TRUE, gdal = TRUE)
  feddata.dem.resample <- resample(feddata.dem.proj, unique.date.crop, "bilinear")
  names(feddata.dem.resample) <- "DEM"

  # Export processed rasters
  writeRaster(
    unique.date.crop,
    file.path(
      main.dir,
      "Planet_Multispectral_Data",
      "3_Level3B_Data_Formatted",
      paste0(unique.date, "_3B_AnalyticMS_SR_8b_clip_processed.tif")
    ),
    overwrite = TRUE
  )

  writeRaster(
    feddata.dem.resample,
    file.path(main.dir, "USGS_DEM_Data", paste0("USGS_DEM_Data_processed_", unique.date, ".tif")),
    overwrite = TRUE
  )

  # Export PlanetScope image identifiers used in annual assessments
  image.identifiers.file <- file.path(
    main.dir,
    "Planet_Multispectral_Data",
    "3_Level3B_Data_Formatted",
    paste0(unique.date, "_PlanetScope_image_identifiers.txt")
  )
  writeLines(substr(names(unique.date.list.nonnull), 1, 23), image.identifiers.file)

  # Progress update
  message(paste("Processing complete for", unique.date))
}
