

## Script to format input data for annual assessments 


## Load required packages 
require(terra)
require(sf)
require(FedData)
require(tools)


## Set main working directory 
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/2_Annual_Assessments/Data/Input_Data/"
sat.dir <- "C:/Users/mmama/Documents/"


## Read in input shapefiles 
tmi.hfa <- st_transform(st_read(paste0(main.dir,"Format_VIMS_Data_For_Planet/TMI_HFA_For_Planet/TMI_HFA_For_Planet.shp")), crs = "EPSG:32618")
salinity <- st_transform(st_read(paste0(main.dir,"Mean_Surface_Salinity/HFA_Mean_Surface_Salinity_1985-2018_formatted.shp")), crs = "EPSG:32618")
tmi.hfa.15ppt <- st_crop(tmi.hfa, ext(salinity))


## Read in all datasets 
## Planet imagery 
zipped.list <- list.files(paste0(sat.dir,"Input_Data/Planet_Multispectral_Data/1_Level_3B_Data_Zipped"), "*.zip", full.names = TRUE)
zipped.path <- lapply(zipped.list, function(x){paste0(sat.dir,"Input_Data/Planet_Multispectral_Data/2_Level_3B_Data_Unzipped/", substr(basename(x), 1, nchar(basename(x)) - 4))})
if(length(zipped.list) > 0){
  for(zipped.file in 1:length(zipped.list)){
    if(!(file.exists(zipped.path[[zipped.file]]))){
      unzip(zipped.list[[zipped.file]], exdir = zipped.path[[zipped.file]])
    }
  } 
}
unzipped.list <- list.files(paste0(sat.dir,"Input_Data/Planet_Multispectral_Data/2_Level_3B_Data_Unzipped/"), full.names = TRUE)
## DEM data 
feddata.get <- get_ned(template = tmi.hfa.15ppt, label = "TMI_HFA", res = 13, extraction.dir = paste0(main.dir,"USGS_DEM_Data/Extracted_DEM_Raster"), force.redo = FALSE)
feddata.dem <- rast(paste0(main.dir,"USGS_DEM_Data/Extracted_DEM_Raster/TMI_HFA_NED_13.tif"))


## Process Planet Multispectral data 
sr.filenames <- list.files(paste0(unzipped.list, "/PSScene"), "*SR_8b.tif$", full.names = T)
unique.dates <- unique(substr(basename(sr.filenames), 1, 8))
## Loop through each unique date and process
for(unique.date in unique.dates){
  
  ## Subset to the date of interest 
  sr.filenames.date <- sr.filenames[which(substr(basename(sr.filenames), 1, 8) == unique.date)]
  unique.date.list <- vector(mode = "list", length = length(sr.filenames.date))
  
  ## Loop through each PlanetScope raster and process 
  for(sr.filename.date in 1:length(sr.filenames.date)){
    
    ## Read in raster
    sr.image <- rast(sr.filenames.date[[sr.filename.date]])
    
    ## Check if it intersects tmi.hfa.15ppt
    if(!(is.null(intersect(ext(sr.image), ext(tmi.hfa.15ppt))))){
      
      ## Read in corresponding UDM raster and remove invalid data
      udm2.image <- rast(paste0(dirname(sr.filenames.date[[sr.filename.date]]), "/",substr(basename(sr.filenames.date[[sr.filename.date]]), 1, 27), "udm2.tif"))
      sr.image[udm2.image$clear == 0] <- NA
      ## Mask to corresponding reference data 
      unique.date.list[[sr.filename.date]] <- mask(sr.image, tmi.hfa.15ppt)
      names(unique.date.list)[[sr.filename.date]] <- file_path_sans_ext(basename(sr.filenames.date[[sr.filename.date]]))
      ## Update progress of the loop
      print(paste("Processing raster", sr.filename.date, "of", length(sr.filenames.date), "within region of interest"))
    
    }
  }
  
  ## Aggregate overlapping pixels for images collected on the same day
  unique.date.list.nonnull <- Filter(Negate(is.null), unique.date.list)
  unique.date.merge <- terra::merge(sprc(unique.date.list.nonnull))
  unique.date.crop <- crop(unique.date.merge, tmi.hfa.15ppt, mask = TRUE)
  
  ## Create corresponding DEM dataset for each unique date
  feddata.dem.proj <- project(feddata.dem, crs(unique.date.crop), mask = TRUE, threads = TRUE, gdal = TRUE)
  feddata.dem.resample <- resample(feddata.dem.proj, unique.date.crop, "bilinear")
  names(feddata.dem.resample) <- "DEM"
  
  ## Export processed rasters
  writeRaster(unique.date.crop, paste0(sat.dir, "Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/", unique.date, "_3B_AnalyticMS_SR_8b_clip_processed.tif"), overwrite = TRUE)
  writeRaster(feddata.dem.resample, paste0(main.dir, "USGS_DEM_Data/", "USGS_DEM_Data_processed_", unique.date, ".tif"), overwrite = TRUE)
  ## Export PlanetScope image identifiers used in annual assessments 
  image.identifiers <- paste0(sat.dir,"Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/", unique.date,"_PlanetScope_image_identifiers.txt")
  sink(image.identifiers)
  print(substr(names(unique.date.list.nonnull),1,23))
  sink()

  ## Update progress of the loop
  print(paste("Processing complete for", unique.date))
  
}

