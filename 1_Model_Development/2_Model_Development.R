

## Script to generate random forest model: 
##    1. Create training dataset 
##    2. Hyperparameter tuning
##    3. Create Random Forest model 


## Load required packages 
require(terra)
require(sf)
require(tools)
require(tidymodels)
require(tidyterra)
require(randomForest)


## Set main working directory 
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/1_Model_Development/"


###############################################
######### 1. Create training dataset ##########
###############################################


## Read in formatted satellite data 
sr.formatted <- lapply(list.files(paste0(main.dir,"Data/Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/"), "*.tif$", full.names = TRUE), rast)
dem.formatted <- lapply(list.files(paste0(main.dir,"Data/Input_Data/USGS_DEM_Data/"), "*.tif$", full.names = TRUE), rast)
## Read in reference data 
ref <- st_read(paste0(main.dir,"Data/Input_Data/VIMS_Reference_Data/Vegetation_shapefile_May_2021.shp"))

## Find the number of pixels that overlap with each vegetation class (images follow the same grid, so just need to assess one raster)
i.class.mask <- lapply(as.list(unique(ref$VEG_CLASS)), function(x){mask(sr.formatted[[1]][[1]], subset(ref, VEG_CLASS == x), touches = FALSE)})
i.pixel.count <- lapply(as.list(c(1:length(i.class.mask))), function(y){lapply(i.class.mask[[y]], function(x){length(which(!is.na(values((x)))))})})
names(i.pixel.count) <- as.list(unique(ref$VEG_CLASS))
## Set the number of pixels per class as 80% of the smallest class 
n.pixels.per.class <- floor(min(unlist(i.pixel.count)) * 0.8) 
## Export the number of pixels per class 
sink(file = paste0(main.dir,"Code/Training_Data/n_pixels_per_class.txt"))
print("n pixels selected per class =")
n.pixels.per.class
print("count of pixels per class =")
unlist(i.pixel.count)
sink(file = NULL)

## Loop through each image and extract training data 
for(image in 1:length(sr.formatted)){
  
  ## Subset to the raster of interest 
  sr.image <- sr.formatted[[image]]
  
  ## Define a dataframe to populate with training data  
  veg.classes <- c("Spartina_alterniflora","Spartina_patens","Phragmites_australis","Juncus_roemerianus","NoVeg")
  training.data <- as.data.frame(matrix(nrow = n.pixels.per.class * length(veg.classes), ncol = 21)) 
  colnames(training.data) <- c("DATE","CLASS","PT_LAT","PT_LON","coastal_blue","blue","green_i","green","yellow","red","rededge","nir","NDVI","EVI","NDWI","SAVI","WDRVI","BG","GR","NIRR","DEM")
  training.data$DATE <- rep(substr(basename(sources(sr.image)), 1, 8), times = n.pixels.per.class * length(veg.classes))
  training.data$CLASS <- rep(veg.classes, times = n.pixels.per.class) 
  
  ## Define a list to populate with random sample coordinates 
  random.sample.list <- vector(mode = "list", length = length(veg.classes))
    
  ## Loop through each vegetation class
  for(veg.class in veg.classes){
    
    ## Subset to the class of interest
    ref.class <- subset(ref, VEG_CLASS == veg.class)
    ## Crop and mask Planet Multispectral data to vims.subset
    sr.crop <- mask(sr.image, ref.class, touches = FALSE)
    sr.extract <- terra::extract(sr.crop, ref.class, cells = T)
    ## Crop and mask DEM data to vims.subset
    dem.crop <- mask(dem.formatted[[image]], ref.class)
    dem.extract <- terra::extract(dem.crop, ref.class, cells = T)
    ## Randomly select n.pixels.per.class pixels
    random.sample <- sample(nrow(sr.extract), size = n.pixels.per.class, replace = T) 
    sr.sample <- sr.extract[random.sample,]
    dem.sample <- dem.extract[random.sample,]
    ## Extract Planet spectral imagery at each of these locations; populate r.df 
    training.data[which(training.data$CLASS == veg.class),"PT_LAT"] <- xyFromCell(sr.crop, sr.sample$cell)[,2]
    training.data[which(training.data$CLASS == veg.class),"PT_LON"] <- xyFromCell(sr.crop, sr.sample$cell)[,1]
    training.data[which(training.data$CLASS == veg.class),c("coastal_blue","blue","green_i","green","yellow","red","rededge","nir")] <- sr.sample[,c("coastal_blue","blue","green_i","green","yellow","red","rededge","nir")]
    training.data[which(training.data$CLASS == veg.class),"DEM"] <- dem.sample$DEM
  }
  
  ## Populate vegetation indices 
  training.data$NDVI  <- (training.data$nir - training.data$red) / (training.data$nir + training.data$red)
  training.data$EVI   <-  2.5 * (training.data$nir - training.data$red) / ((training.data$nir + 6 * training.data$red - 7.5 * training.data$blue) + 1)
  training.data$NDWI  <- (training.data$green - training.data$nir) / (training.data$green + training.data$nir)
  training.data$SAVI  <- ((training.data$nir - training.data$red) * 1.5) / (training.data$nir + training.data$red + 0.5)
  training.data$WDRVI <- (0.2 * training.data$nir - training.data$red) / (0.2 * training.data$nir + training.data$red)
  training.data$BG    <- training.data$blue / training.data$green
  training.data$GR    <- training.data$green / training.data$red
  training.data$NIRR  <- training.data$nir / training.data$red
  ## Remove any rows containing NA values
  training.data.cc <- training.data[complete.cases(training.data),]
  ## Export training dataset 
  write.csv(training.data.cc, paste0(main.dir, "Code/Training_Data/", substr(basename(sources(sr.image)), 1, 8), "_Training_Data.csv"), row.names = FALSE)
  
  ## Export as spatial object 
  random.sample.coords <- st_as_sf(training.data.cc, coords = c("PT_LON","PT_LAT"), crs = crs(sr.crop))
  write_sf(random.sample.coords, paste0(main.dir, "Code/Training_Data/Training_Data_Shapefile/", substr(basename(sources(sr.image)), 1, 8), "_Training_Data.shp"), append = FALSE)
}

## Remove variables 
rm(sr.formatted, dem.formatted, i.class.mask, i.pixel.count, n.pixels.per.class, image, 
   sr.image, veg.classes, training.data, veg.class, ref.class, sr.crop, sr.extract, 
   dem.crop, dem.extract, random.sample, sr.sample, dem.sample, training.data.cc,
   random.sample.list, random.sample.coords, sr.sample.xy)


###############################################
########## 2. Hyperparameter tuning ###########
###############################################

## Read in training data 
training.data.filenames <- list.files(paste0(main.dir,"Code/Training_Data/"), "*.csv$", full.names = TRUE)
## Define a dataframe to populate with results 
hyperparameters <- as.data.frame(matrix(nrow = length(training.data.filenames), ncol = 4))
colnames(hyperparameters) <- c("IMAGE","OPTIMAL_TREES","OPTIMAL_MTRY","OPTIMAL_MIN_N")
hyperparameters$IMAGE <- substr(basename(training.data.filenames), 1, 8)

## Loop through each training dataset and tune hyperparameters
for(image in 1:length(training.data.filenames)){
  
  ## Read in training data CSV
  training.data.image <- read.csv(training.data.filenames[[image]], header = T)
  ## Split data 75/25 for training and testing 
  data.split <- initial_split(training.data.image[,c("CLASS","coastal_blue","blue","green_i","green","yellow","red","rededge","nir","NDVI","EVI","NDWI","SAVI","WDRVI","BG","GR","NIRR","DEM")])
  
  ## Apply dataset pre-processing
  data.recipe <- training(data.split) %>% recipe(CLASS ~.) %>%
    step_corr(all_predictors(), threshold = 0.95, method = "spearman") %>%
    step_center(all_predictors(), -all_outcomes()) %>% step_scale(all_predictors(), -all_outcomes()) 
  ## Set up object to tune based on mtry and min_n 
  data.rf.tune <- rand_forest(trees = tune(), mode = "classification", mtry = tune(), min_n = tune()) %>% set_engine("randomForest")
  ## Set up workflow object 
  tuning.workflow <- workflow() %>% add_recipe(data.recipe) %>% add_model(data.rf.tune)
  ## Define resampling strategy and tune parameters 
  data.vfold <- vfold_cv(training(data.split), v = 5) 
  tune.resample <- tune_grid(tuning.workflow, resamples = data.vfold, grid = 10)
  ## Choose the model with the highest performance 
  best.model <- select_best(tune.resample, metric = "roc_auc")
  ## Populate dataframe with results 
  hyperparameters$OPTIMAL_TREES[image] <- best.model$trees
  hyperparameters$OPTIMAL_MTRY[image]  <- best.model$mtry
  hyperparameters$OPTIMAL_MIN_N[image] <- best.model$min_n
  
}

## Export hyperparameters dataframe as CSV file 
write.csv(hyperparameters, paste0(main.dir,"Code/Random_Forest_Model/Optimal_hyperparameters.csv"), row.names = F)

## Remove variables 
rm(training.data.filenames, hyperparameters, image, training.data.image, data.split, data.recipe, data.rf.tune, 
   tuning.workflow, data.vfold, tune.resample, best.model)


###############################################
######## 3. Create Random Forest model ########
###############################################


## Read in formatted satellite data and training data 
sr.formatted <- lapply(list.files(paste0(main.dir,"Data/Input_Data/Planet_Multispectral_Data/3_Level3B_Data_Formatted/"), "*.tif$", full.names = TRUE), rast)
training.data.filenames <- list.files(paste0(main.dir,"Code/Training_Data/"), "*_Training_Data.csv$", full.names = TRUE)
optimized.hyperparameters <- read.csv(paste0(main.dir,"Code/Random_Forest_Model/Optimal_hyperparameters.csv"), header = T)

## Loop through each image and build random forest model 
for(image in 1:length(sr.formatted)){
  
  ## Read in training data CSV
  training.data.image <- read.csv(training.data.filenames[[image]], header = T)
  ## Split data 75/25 for training and testing and set seed for reproducibility 
  data.split <- initial_split(training.data.image[,c("CLASS","coastal_blue","blue","green_i","green","yellow","red","rededge","nir","NDVI","EVI","NDWI","SAVI","WDRVI","BG","GR","NIRR","DEM")])
  
  ## Apply data pre-processing 
  data.recipe <- training(data.split) %>% recipe(CLASS ~.) %>%
    step_corr(all_predictors(), threshold = 0.95, method = "spearman") %>%
    step_center(all_predictors(), -all_outcomes()) %>% step_scale(all_predictors(), -all_outcomes())
  data.prep <- data.recipe %>% prep()
  data.testing <- data.prep %>% bake(testing(data.split)) 
  data.training <- bake(data.prep, new_data = NULL)
  ## Train a random forest model specifying the hyperparameters selected below 
  optimized.hyperparameters.image <- subset(optimized.hyperparameters, IMAGE == substr(basename(training.data.filenames[[image]]), 1, 8))
  data.rf <- rand_forest(trees = optimized.hyperparameters.image$OPTIMAL_TREES, mode = "classification", mtry = optimized.hyperparameters.image$OPTIMAL_MTRY, min_n = optimized.hyperparameters.image$OPTIMAL_MIN_N) %>% 
    set_engine("randomForest") %>% fit(CLASS ~ ., data = data.training)
  
  ## Compute variable importance as mean decrease in node impurity and export
  var.imp <- as.data.frame(importance(data.rf$fit, type = 2))
  write.csv(var.imp, paste0(main.dir,"Code/Variable_Importance/", substr(basename(training.data.filenames[[image]]), 1, 8), "_variable_importance.csv"))

  ## Export model and data preparation information 
  saveRDS(data.prep, file = paste0(main.dir,"Code/Random_Forest_Model/", substr(basename(training.data.filenames[[image]]), 1, 8), "_data_prep.RData"))
  saveRDS(data.rf, file = paste0(main.dir,"Code/Random_Forest_Model/", substr(basename(training.data.filenames[[image]]), 1, 8), "_random_forest_model.RData"))
  
}

## Remove variables 
rm(sr.formatted, training.data.filenames, optimized.hyperparameters, image, 
   training.data.image, data.split, data.recipe, data.prep, data.testing, data.training, var.imp, 
   optimized.hyperparameters.image, data.rf, ref, main.dir)
