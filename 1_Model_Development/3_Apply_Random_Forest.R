

## Script to apply random forest to all reference data 


## Set main working directory 
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/1_Model_Development/"


## Load required packages 
require(terra)
require(sf)
require(tidyterra)
require(tidymodels)


## Set confidence level for random forest probabilities 
conf.level <- 0.7
## Read in formatted satellite data and training data 
sr.formatted <- lapply(list.files(paste0(main.dir,"Data/Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/"), "*.tif$", full.names = TRUE), rast)
dem.formatted <- lapply(list.files(paste0(main.dir,"Data/Input_Data/USGS_DEM_Data/"), "*.tif$", full.names = TRUE), rast)
training.data.coords <- lapply(list.files(paste0(main.dir, "Code/Training_Data/Training_Data_Shapefile/"), "*_Training_Data.shp$", full.names = TRUE), st_read)
training.data.coords <- lapply(training.data.coords, function(x){st_transform(x, crs(sr.formatted[[1]]))})
training.data.filenames <- list.files(paste0(main.dir,"Code/Training_Data/"), "*_Training_Data.csv$", full.names = TRUE)
data.preps <- list.files(paste0(main.dir, "Code/Random_Forest_Model/"), "*_data_prep.RData$", full.names = TRUE)
rf.models <- list.files(paste0(main.dir, "Code/Random_Forest_Model/"), "*_random_forest_model.RData$", full.names = TRUE)


## Loop through each image and build random forest model 
for(image in 1:length(sr.formatted)){
  
  ## Subset to the image of interest and add layer for training data 
  training.data.rast <- rasterize(training.data.coords[[image]], sr.formatted[[image]], background = 0)
  names(training.data.rast) <- "training_data"
  formatted.image <- c(sr.formatted[[image]], crop(dem.formatted[[image]], ext(sr.formatted[[image]])), training.data.rast)
  ## Read in corresponding random forest data 
  data.preps.image <- readRDS(data.preps[[image]])
  rf.models.image <- readRDS(rf.models[[image]])
  
  ## Extract coordinates from image 
  sr.coords <- xyFromCell(formatted.image, 1:ncell(formatted.image))
  formatted.image$x <- sr.coords[,1]
  formatted.image$y <- sr.coords[,2]
  ## Convert to a tibble and compute vegetation indices 
  tibble.image <- tidyterra::as_tibble(formatted.image) %>% filter(complete.cases(.))
  tibble.image$NDVI  <- (tibble.image$nir - tibble.image$red) / (tibble.image$nir + tibble.image$red)
  tibble.image$EVI   <-  2.5 * (tibble.image$nir - tibble.image$red) / ((tibble.image$nir + 6 * tibble.image$red - 7.5 * tibble.image$blue) + 1)
  tibble.image$NDWI  <- (tibble.image$green - tibble.image$nir) / (tibble.image$green + tibble.image$nir)
  tibble.image$SAVI  <- ((tibble.image$nir - tibble.image$red) * 1.5) / (tibble.image$nir + tibble.image$red + 0.5)
  tibble.image$WDRVI <- (0.2 * tibble.image$nir - tibble.image$red) / (0.2 * tibble.image$nir + tibble.image$red)
  tibble.image$BG    <- tibble.image$blue / tibble.image$green
  tibble.image$GR    <- tibble.image$green / tibble.image$red 
  tibble.image$NIRR  <- tibble.image$nir / tibble.image$red 
  tibble.image$CLASS <- tibble.image$VEG_CLASS
  tibble.image$VEG_CLASS <- NULL
  
  ## Apply random forest model 
  prep.data.image <- data.preps.image %>% bake(tibble.image) 
  predict.data.image <- rf.models.image %>% predict(prep.data.image) %>% bind_cols(prep.data.image)
  ## Compute probabilities for each class and append to random forest results 
  predict.data.image <- predict(rf.models.image, prep.data.image, type = "prob") %>% bind_cols(predict.data.image)
  ## Create confidence-based columns 
  predict.data.image$MAX_CONF <- pmax(predict.data.image$.pred_Juncus_roemerianus, predict.data.image$.pred_NoVeg, predict.data.image$.pred_Phragmites_australis, predict.data.image$.pred_Spartina_alterniflora, predict.data.image$.pred_Spartina_patens, na.rm = TRUE)
  ## Add band based on given confidence threshold 
  predict.data.image$PRED_CLASS_CONF <- ifelse(predict.data.image$MAX_CONF >= conf.level, as.character(predict.data.image$.pred_class), "Low_confidence")
  
  ## Convert classification to raster 
  predict.data.image$x <- tibble.image$x
  predict.data.image$y <- tibble.image$y
  classified.image <- as_spatraster(predict.data.image, xycols = c(which(colnames(predict.data.image) == "x"), which(colnames(predict.data.image) == "y")), crs = crs(sr.formatted[[image]]))
  ## Export raster and classification key
  writeRaster(classified.image, paste0(main.dir, "Data/Output_Data/Classified_Data/", substr(basename(sources(sr.formatted[[image]])), 1, 8), "_classification_reference_data.tif"), overwrite = TRUE)
  sink(file = paste0(main.dir,"Data/Output_Data/Classified_Data/",substr(basename(sources(sr.formatted[[image]])), 1, 8), "_classification_key.txt"))
  paste("Classification key for PRED_CLASS_CONF:")
  unique(classified.image[["PRED_CLASS_CONF"]])
  sink(file = NULL)
  
  ## Perform sensitivity analysis on confidence threshold 
  conf.df <- as.data.frame(matrix(nrow = 5, ncol = 3))
  colnames(conf.df) <- c("CONF_LEVEL","AREA_LOW_CONF","PCT_LOW_CONF")
  conf.df$CONF_LEVEL <- seq(0.5, 0.9, by = 0.1)
  for(conf in conf.df$CONF_LEVEL){
    classified.conf <- ifelse(predict.data.image$MAX_CONF >= conf, levels(predict.data.image$.pred_class)[predict.data.image$.pred_class], "Low_confidence")
    conf.df$AREA_LOW_CONF[which(conf.df$CONF_LEVEL == conf)] <- round((length(which(classified.conf == "Low_confidence")) * 9))
    conf.df$PCT_LOW_CONF[which(conf.df$CONF_LEVEL == conf)] <- round((length(which(classified.conf == "Low_confidence")) / nrow(predict.data.image)) * 100, digits = 2)
  }
  write.csv(conf.df, paste0(main.dir,"Code/Agreement_Statistics/",substr(basename(sources(sr.formatted[[image]])), 1, 8),"_conf_threshold_iterations.csv"), row.names = F)

}


