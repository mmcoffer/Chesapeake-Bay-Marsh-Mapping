

## Apply random forest classifier to satellite imagery for annual assessments 


## Load required packages 
require(terra)
require(sf)
require(tidymodels)
require(tidyterra)


## Set main working directory 
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/"


## List formatted satellite data and read in random forest model 
sr.directories <- list.files(paste0(main.dir,"2_Annual_Assessments/Data/Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/"), "*", full.names = TRUE)
dem.filenames <- list.files(paste0(main.dir,"2_Annual_Assessments/Data/Input_Data/USGS_DEM_Data/"), "*.tif$", full.names = TRUE)
data.prep <- readRDS(list.files(paste0(main.dir, "1_Model_Development/Code/Random_Forest_Model/"), "20210506_data_prep.RData$", full.names = TRUE))
rf.model <- readRDS(list.files(paste0(main.dir, "1_Model_Development/Code/Random_Forest_Model/"), "20210506_random_forest_model.RData$", full.names = TRUE))


## Loop through each unique date and apply classification 
for(unique.date in 1:length(basename(sr.directories))){
  
  ## List files for the given date 
  sr.files <- list.files(sr.directories[[unique.date]], "*.tif$", full.names = TRUE)
  dem.file <- rast(dem.filenames[[unique.date]])
  
  ## Loop through each file and process 
  for(sr.file in 1:length(sr.files)){
    
    ## Read in the raster of interest and append DEM 
    sr.rast <- rast(sr.files[[sr.file]])
    formatted.image <- c(sr.rast, crop(dem.file, ext(sr.rast)))
    
    ## Extract coordinates from image 
    sr.coords <- xyFromCell(formatted.image, 1:ncell(sr.rast))
    formatted.image$x <- sr.coords[,1]
    formatted.image$y <- sr.coords[,2]
    ## Convert to a tibble and compute vegetation indices 
    tibble.image <- tidyterra::as_tibble(formatted.image) %>% filter(complete.cases(.)) %>% filter_all(all_vars(!is.infinite(.)))
    tibble.image$NDVI  <- (tibble.image$nir - tibble.image$red) / (tibble.image$nir + tibble.image$red)
    tibble.image$EVI   <-  2.5 * (tibble.image$nir - tibble.image$red) / ((tibble.image$nir + 6 * tibble.image$red - 7.5 * tibble.image$blue) + 1)
    tibble.image$NDWI  <- (tibble.image$green - tibble.image$nir) / (tibble.image$green + tibble.image$nir)
    tibble.image$SAVI  <- ((tibble.image$nir - tibble.image$red) * 1.5) / (tibble.image$nir + tibble.image$red + 0.5)
    tibble.image$WDRVI5 <- (0.2 * tibble.image$nir - tibble.image$red) / (0.2 * tibble.image$nir + tibble.image$red)
    tibble.image$BG    <- tibble.image$blue / tibble.image$green
    tibble.image$GR    <- tibble.image$green / tibble.image$red 
    tibble.image$NIRR  <- tibble.image$nir / tibble.image$red 
    tibble.image$VEG_CLASS <- NULL
    
    ## Apply random forest model 
    prep.data.image <- data.prep %>% bake(tibble.image) %>% filter(complete.cases(.)) %>% filter_all(all_vars(!is.infinite(.)))
    predict.data.image <- rf.model %>% predict(prep.data.image) %>% bind_cols(prep.data.image)
    ## Compute probabilities for each class and append to random forest results 
    predict.data.image <- predict(rf.model, prep.data.image, type = "prob") %>% bind_cols(predict.data.image)
    ## Create confidence-based columns 
    conf.data.image <- predict.data.image %>% mutate(MAX_CONFIDENCE = pmax(predict.data.image$.pred_Juncus_roemerianus, predict.data.image$.pred_NoVeg, predict.data.image$.pred_Phragmites_australis, predict.data.image$.pred_Spartina_alterniflora, predict.data.image$.pred_Spartina_patens, na.rm = TRUE))
    conf.data.image$.pred_class_confidence <- ifelse(conf.data.image$MAX_CONFIDENCE >= 0.5, levels(conf.data.image$.pred_class)[conf.data.image$.pred_class], "Low confidence")
    
    ## Convert classification to raster 
    predict.data.image$x <- tibble.image$x
    predict.data.image$y <- tibble.image$y
    classified.image <- as_spatraster(predict.data.image, xycols = c(which(colnames(predict.data.image) == "x"), which(colnames(predict.data.image) == "y")), crs = crs(sr.rast))
    ## Apply smoothing via mode with 3 x 3 moving window  
    smoothed.image <- terra::focal(classified.image[[".pred_class"]], w = matrix(1, nc = 3, nr = 3), fun = "modal")
    levels(smoothed.image) <- data.frame(id = 1:5, .pred_class_smooth = unique(classified.image$.pred_class)$.pred_class)
    stacked.image <- c(classified.image, smoothed.image)
    ## Export rasters -- create folder for date if it doesn't already exist
    if(!(dir.exists(paste0(main.dir,"2_Annual_Assessments/Data/Output_Data/Classified_Data_Per_Date/", basename(sr.directories)[unique.date])))){
      dir.create(paste0(main.dir,"2_Annual_Assessments/Data/Output_Data/Classified_Data_Per_Date/", basename(sr.directories)[unique.date]))
    }
    writeRaster(stacked.image, paste0(main.dir,"2_Annual_Assessments/Data/Output_Data/Classified_Data_Per_Date/", basename(sr.directories)[unique.date], "/", basename(sr.directories)[unique.date], "_classified_image_", sr.file,".tif"), overwrite = TRUE)
    
    
  }
  
}



