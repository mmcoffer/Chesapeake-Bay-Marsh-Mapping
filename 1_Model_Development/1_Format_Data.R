

## Script to format satellite data:
##    1. Planet multispectral data
##    2. USGS DEM data  


## Load required packages 
require(terra)
require(sf)
require(FedData)
require(tools)


## Set main working directory 
#main.dir <- "C:/Users/mmama/Documents/Input_Data/"
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/1_Model_Development/Data/Input_Data/"


## Read in all datasets 
## Planet imagery 
zipped.list <- list.files(paste0(main.dir,"Planet_Multispectral_Data/1_Level_3B_Data_Zipped"), "*.zip", full.names = TRUE)
zipped.path <- lapply(zipped.list, function(x){paste0(main.dir,"Planet_Multispectral_Data/2_Level_3B_Data_Unzipped/", substr(basename(x), 1, nchar(basename(x)) - 4))})
for(zipped.file in 1:length(zipped.list)){
  if(!(file.exists(zipped.path[[zipped.file]]))){
    unzip(zipped.list[[zipped.file]], exdir = zipped.path[[zipped.file]])
  }
}
## List unzipped data 
unzipped.list <- list.files(paste0(main.dir,"Planet_Multispectral_Data/2_Level_3B_Data_Unzipped/"), full.names = TRUE)
## Reference data 
ref <- st_read(paste0(main.dir,"VIMS_Reference_Data/Vegetation_shapefile_May_2021.shp"))
## DEM data 
feddata.get <- get_ned(template = ref, label = "VIMS_CCRM", res = 13, extraction.dir = paste0(main.dir,"USGS_DEM_Data/Extracted_DEM_Raster"), force.redo = TRUE)
feddata.dem <- rast(paste0(main.dir,"USGS_DEM_Data/Extracted_DEM_Raster/VIMS_CCRM_NED_13.tif"))
  

###############################################
######## 1. Planet multispectral data #########
###############################################


## Process Planet data for May 2021 
sr.filenames <- list.files(paste0(unzipped.list, "/PSScene"), "*SR_8b.tif$", full.names = T)
sr.masked <- vector(mode = "list", length = length(sr.filenames))
## Loop through each Planet file and process 
for(image in 1:length(sr.filenames)){
  
  ## Read in raster, corresponding UDM raster, and remove invalid data
  sr.image <- rast(sr.filenames[[image]])
  udm2.image <- rast(paste0(dirname(sr.filenames[[image]]), "/",substr(basename(sr.filenames[[image]]), 1, 27), "udm2.tif"))
  sr.image[udm2.image$clear == 0] <- NA
  
  ## Mask to corresponding reference data 
  sr.masked[[image]] <- mask(sr.image, ref)
  names(sr.masked)[[image]] <- file_path_sans_ext(basename(sr.filenames[[image]]))
  
}

## Find unique dates for Planet data
unique.dates <- unique(substr(basename(sr.filenames), 1, 8))
sr.formatted <- vector(mode = "list", length = length(unique.dates))
## Combine Planet data from the same dates 
for(unique.date in 1:length(unique.dates)){
  
  ## Take the mean of overlapping pixels for images collected on the same day
  i.dates <- sr.masked[which(substr(names(sr.masked), 1, 8) == unique.dates[unique.date])]
  unique.date.merge <- merge(sprc(i.dates))
  sr.formatted[[unique.date]] <- crop(unique.date.merge, ext(ref))
  names(sr.formatted)[[unique.date]] <- unique.dates[unique.date]
  ## Rasterize reference data and add layer for corresponding reference data
  ref.rast <- rasterize(ref, sr.formatted[[unique.date]], "VEG_CLASS")
  sr.formatted[[unique.date]] <- c(sr.formatted[[unique.date]], ref.rast)
  ## Export raster 
  writeRaster(sr.formatted[[unique.date]], paste0(main.dir, "Planet_Multispectral_Data/3_Level3B_Data_Formatted/", unique.dates[[unique.date]], "_3B_AnalyticMS_SR_8b_clip_processed.tif"), overwrite = TRUE)
  
}


###############################################
################ 2. DEM data ##################
###############################################


## Create corresponding DEM dataset for each unique date
for(unique.date in 1:length(unique.dates)){
  
  ## Project and resample to match Planet data 
  feddata.dem.proj <- project(feddata.dem, crs(sr.formatted[[unique.date]]))
  feddata.dem.resample <- resample(feddata.dem.proj, sr.formatted[[unique.date]], "bilinear")
  ## Crop and mask to reference data 
  feddata.dem.ref <- crop(feddata.dem.resample, ref, mask = TRUE)
  ## Rename raster data 
  names(feddata.dem.ref) <- "DEM"
  ## Export processed raster 
  writeRaster(feddata.dem.ref, paste0(main.dir, "USGS_DEM_Data/", "USGS_DEM_Data_processed_", unique.dates[[unique.date]], ".tif"), overwrite = TRUE)
  
}

