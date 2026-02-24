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
#' 2. Apply model to all dates, grouped by year
#' 3. Generate per-class probability predictions for each tile and date
#' 4. Mosaic tiles within each date
#' 5. Compute pixel-wise median class probabilities across all dates in a year
#' 6. Assign final classification as the class with the highest median probability
#' 7. Export a single multi-band annual raster per year
#'
#' @inputs
#' - Formatted Planet imagery rasters (filenames starting with YYYYMMDD)
#' - Formatted DEM rasters (filenames containing YYYYMMDD)
#' - Trained Random Forest model (.RData)
#' - Data preprocessing object (.RData)
#'
#' @outputs
#' - One GeoTIFF per year: classification band (VEG_CLASS) + 5 median probability bands


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

# Probability class column names (must match tidymodels .pred_* output)
prop.class.names <- c(
  ".pred_Juncus_roemerianus",
  ".pred_NoVeg",
  ".pred_Phragmites_australis",
  ".pred_Spartina_alterniflora",
  ".pred_Spartina_patens"
)

# Determine unique dates (YYYYMMDD) and years (YYYY) from sr.filenames
unique.dates <- unique(substr(basename(sr.filenames), 1, 8))
unique.years <- unique(substr(unique.dates, 1, 4))

# Create annual output directory
annual.output.dir <- file.path(
  main.dir,
  "2_Annual_Assessments", "Data", "Output_Data", "Classified_Data_Annual"
)
if (!dir.exists(annual.output.dir)) {
  dir.create(annual.output.dir, recursive = TRUE)
}

# ==== Loop over years =========================================================
for (unique.year in unique.years) {

  message(paste("Processing year:", unique.year))

  # Dates belonging to this year
  year.dates <- unique.dates[substr(unique.dates, 1, 4) == unique.year]

  # Paths to temp files holding one mosaicked probability raster per date
  date.prob.raster.files <- vector("character", length(year.dates))

  # ==== Loop over dates within the year =======================================
  for (date.idx in seq_along(year.dates)) {
    current.date <- year.dates[date.idx]

    message(paste("  Processing date:", current.date))

    # All tiles for this date
    date.sr.files <- sr.filenames[substr(basename(sr.filenames), 1, 8) == current.date]

    # Match DEM by 8-char date prefix in DEM filename
    dem.match <- grep(current.date, basename(dem.filenames), value = FALSE)
    if (length(dem.match) == 0) {
      warning(paste("No DEM found for date", current.date, "-- skipping."))
      next
    }
    dem.file <- rast(dem.filenames[dem.match[1]])

    # Collect probability rasters for each tile
    tile.prob.rasters <- vector("list", length(date.sr.files))

    # ==== Loop over tiles within the date =====================================
    for (tile.idx in seq_along(date.sr.files)) {

      # Read raster and append DEM
      sr.rast <- rast(date.sr.files[tile.idx])
      formatted.image <- c(sr.rast, crop(dem.file, ext(sr.rast)))

      # Split tile into row-blocks so each block fits in RAM individually.
      # terra::blocks() returns row offsets sized to available memory.
      blk <- blocks(sr.rast, n = 4)
      block.prob.list <- vector("list", blk$n)

      # ==== Loop over row-blocks within the tile ==============================
      for (b in seq_len(blk$n)) {

        # Crop formatted image to this block's row extent
        row.start <- blk$row[b]
        row.end   <- blk$row[b] + blk$nrows[b] - 1
        blk.ext   <- ext(
          xmin(formatted.image),
          xmax(formatted.image),
          yFromRow(formatted.image, row.end),
          yFromRow(formatted.image, row.start)
        )
        blk.rast <- crop(formatted.image, blk.ext)

        # Extract coordinates for this block
        blk.coords   <- xyFromCell(blk.rast, 1:ncell(blk.rast))
        blk.rast$x   <- blk.coords[, 1]
        blk.rast$y   <- blk.coords[, 2]
        rm(blk.coords)

        # Convert block to data frame and filter
        tibble.block <- as.data.frame(blk.rast, na.rm = FALSE) %>%
          filter(complete.cases(.)) %>%
          filter_all(all_vars(!is.infinite(.)))
        rm(blk.rast)

        if (nrow(tibble.block) == 0) {
          rm(tibble.block); gc(); next
        }

        # Compute vegetation indices
        tibble.block$NDVI  <- (tibble.block$nir - tibble.block$red) / (tibble.block$nir + tibble.block$red)
        tibble.block$EVI   <- 2.5 * (tibble.block$nir - tibble.block$red) / ((tibble.block$nir + 6 * tibble.block$red - 7.5 * tibble.block$blue) + 1)
        tibble.block$NDWI  <- (tibble.block$green - tibble.block$nir) / (tibble.block$green + tibble.block$nir)
        tibble.block$SAVI  <- ((tibble.block$nir - tibble.block$red) * 1.5) / (tibble.block$nir + tibble.block$red + 0.5)
        tibble.block$WDRVI <- (0.2 * tibble.block$nir - tibble.block$red) / (0.2 * tibble.block$nir + tibble.block$red)
        tibble.block$BG    <- tibble.block$blue / tibble.block$green
        tibble.block$GR    <- tibble.block$green / tibble.block$red
        tibble.block$NIRR  <- tibble.block$nir / tibble.block$red
        tibble.block$VEG_CLASS <- NULL

        # Retain x/y before baking -- the recipe may drop them as non-features.
        # bake() preserves row count and order, so we can safely re-attach
        # coordinates onto the baked result before filtering. That way x/y
        # survive the post-bake complete.cases/!is.infinite filter in lockstep
        # with the predictor columns, with no row-alignment risk.
        block.xy <- tibble.block %>% select(x, y)

        # Bake preprocessing recipe
        baked.block <- data.prep %>% bake(tibble.block)
        rm(tibble.block)

        # Re-attach x/y if the recipe dropped them
        if (!"x" %in% names(baked.block)) baked.block$x <- block.xy$x
        if (!"y" %in% names(baked.block)) baked.block$y <- block.xy$y
        rm(block.xy)

        prep.block <- baked.block %>%
          filter(complete.cases(.)) %>%
          filter_all(all_vars(!is.infinite(.)))
        rm(baked.block)

        if (nrow(prep.block) == 0) {
          rm(prep.block); gc(); next
        }

        # Compute class probabilities and retain x/y coordinates
        prob.block <- predict(rf.model, prep.block, type = "prob") %>%
          bind_cols(prep.block %>% select(x, y))
        rm(prep.block)

        # Convert block probabilities to SpatRaster
        block.prob.rast <- as_spatraster(
          prob.block,
          xycols = c(which(colnames(prob.block) == "x"),
                     which(colnames(prob.block) == "y")),
          crs = crs(sr.rast)
        )
        names(block.prob.rast) <- prop.class.names
        block.prob.list[[b]] <- block.prob.rast
        rm(prob.block, block.prob.rast)
        gc()
      }
      # ==== End block loop ====================================================

      # Mosaic row-blocks back into a single tile probability raster
      block.prob.list <- Filter(Negate(is.null), block.prob.list)
      if (length(block.prob.list) == 0) {
        rm(block.prob.list); next
      } else if (length(block.prob.list) == 1) {
        tile.prob.rasters[[tile.idx]] <- block.prob.list[[1]]
      } else {
        tile.prob.rasters[[tile.idx]] <- mosaic(sprc(block.prob.list), fun = "mean")
      }
      rm(block.prob.list)
      gc()
    }
    # ==== End tile loop =======================================================

    # Mosaic tiles for this date into a single probability raster
    tile.prob.rasters <- Filter(Negate(is.null), tile.prob.rasters)
    if (length(tile.prob.rasters) == 1) {
      date.mosaic <- tile.prob.rasters[[1]]
    } else {
      date.mosaic <- mosaic(sprc(tile.prob.rasters), fun = "mean")
    }
    names(date.mosaic) <- prop.class.names

    # Write mosaicked date raster to a temp file and release from RAM
    tmp.file <- tempfile(fileext = ".tif")
    writeRaster(date.mosaic, tmp.file, overwrite = TRUE)
    date.prob.raster.files[date.idx] <- tmp.file
    rm(date.mosaic, tile.prob.rasters)
    gc()
  }
  # ==== End date loop =========================================================

  # Remove empty entries (dates skipped due to missing DEM)
  date.prob.raster.files <- date.prob.raster.files[nchar(date.prob.raster.files) > 0]

  if (length(date.prob.raster.files) == 0) {
    warning(paste("No valid dates processed for year", unique.year, "-- skipping."))
    next
  }

  # ==== Compute pixel-wise median class probabilities across all dates ========
  # terra::union(SpatExtent, SpatExtent) takes exactly two arguments, so use
  # Reduce() to accumulate the union pairwise across all date extents.
  union.ext <- Reduce(terra::union, lapply(date.prob.raster.files, function(f) ext(rast(f))))

  # For each class, load one band at a time from each date file, extend to
  # union extent, stack, and compute the median. This keeps only one
  # single-class stack in RAM at a time rather than all dates x all bands.
  class.medians <- lapply(seq_along(prop.class.names), function(ci) {
    class.stack <- rast(lapply(date.prob.raster.files, function(f) {
      extend(rast(f)[[ci]], union.ext)
    }))
    app(class.stack, function(x) median(x, na.rm = TRUE))
  })

  year.median.probs <- rast(class.medians)
  rm(class.medians)
  gc()

  names(year.median.probs) <- prop.class.names

  # ==== Assign final classification: class with highest median probability ====
  year.classification <- which.lyr(year.median.probs)
  class.labels <- gsub("^\\.pred_", "", prop.class.names)
  levels(year.classification) <- list(data.frame(id = 1:5, label = class.labels))
  names(year.classification) <- "VEG_CLASS"

  # ==== Combine classification and probability bands and export ===============
  annual.raster <- c(year.classification, year.median.probs)

  writeRaster(
    annual.raster,
    file.path(annual.output.dir, paste0(unique.year, "_annual_classification.tif")),
    overwrite = TRUE
  )

  message(paste("  Annual raster written for year", unique.year))

  # Clean up temp files and free memory before the next year
  unlink(date.prob.raster.files)
  rm(annual.raster, year.median.probs, year.classification, date.prob.raster.files)
  gc()
}
# ==== End year loop ===========================================================
