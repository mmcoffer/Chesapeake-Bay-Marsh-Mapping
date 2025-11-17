

## Script to assess agreement with reference data 


## Set main working directory 
main.dir <- "C:/Users/mmama/Documents/COASTAL_WATER/Planet_Marsh_Mapping/1_Model_Development/"


## Load required packages 
require(terra)
require(sf)
require(tidyterra)
require(tidymodels)


## Read in classified satellite data
classified.images <- lapply(list.files(paste0(main.dir,"Data/Output_Data/Classified_Data/"), "*.tif$", full.names = TRUE), rast)
ref <- st_read(paste0(main.dir,"Data/Input_Data/VIMS_Reference_Data/Vegetation_shapefile_May_2021.shp"))
regions <- st_transform(st_read(paste0(main.dir,"Data/Input_Data/VIMS_Reference_Data/Coffer-defined_marsh_regions.shp")), st_crs(ref))
ref$REGION <- st_intersection(regions,ref)$Name


## Loop through each image and assess agreement 
for(image in 2:length(classified.images)){
  
  ## Subset to the image of interest
  classified.image <- classified.images[[image]]
  ## Format data for assessing agreement 
  classified.image$PRED_CLASS_CONF[classified.image$PRED_CLASS_CONF == "Low_confidence"] <- NA
  agreement.tibble <- tidyterra::as_tibble(classified.image) %>% filter(if_all(c(CLASS,PRED_CLASS_CONF), complete.cases))
  agreement.tibble$PRED_CLASS_CONF <- factor(agreement.tibble$PRED_CLASS_CONF, levels = c("Spartina_alterniflora","Phragmites_australis","Juncus_roemerianus","Spartina_patens","NoVeg"))
  agreement.tibble$CLASS <- factor(agreement.tibble$CLASS, levels = c("Spartina_alterniflora","Phragmites_australis","Juncus_roemerianus","Spartina_patens","NoVeg"))
  multiclass.metrics <- metric_set(sensitivity, specificity, bal_accuracy, accuracy)
  
  ## Agreement across all reference data
  ## Create dataframe to populate with results  
  agreement.mat <- conf_mat(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_CONF)
  agreement.df <- as.data.frame(matrix(nrow = 5, ncol = 5, as.data.frame(tidy(agreement.mat))$value))
  agreement.df <- cbind(names(agreement.mat[[1]][,1]), agreement.df)
  colnames(agreement.df) <- c("PREDICTED", names(agreement.mat[[1]][,1]))
  ## Class-by-class agreement statistics 
  per.class.metrics <- agreement.tibble %>% group_by(CLASS) %>% multiclass.metrics(truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")
  agreement.df$CLASS_SENSITIVITY <- round(subset(per.class.metrics, .metric == "sensitivity")$.estimate * 100)
  agreement.df$CLASS_SPECIFICITY <- round(subset(per.class.metrics, .metric == "specificity")$.estimate * 100)
  agreement.df$CLASS_BAL_AGREEMENT <- round(subset(per.class.metrics, .metric == "bal_accuracy")$.estimate * 100)
  agreement.df$CLASS_AGREEMENT <- round(subset(per.class.metrics, .metric == "accuracy")$.estimate * 100)
  ## Overall agreement statistics 
  agreement.df$OVERALL_SENSITIVITY <- rep(round(sensitivity(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100), each = 5) 
  agreement.df$OVERALL_SPECIFICITY <- rep(round(specificity(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100), each = 5)
  agreement.df$OVERALL_BAL_AGREEMENT <- rep(round(bal_accuracy(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100), each = 5)
  agreement.df$OVERALL_AGREEMENT <- rep(round(accuracy(agreement.tibble, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100), each = 5)
  ## Export dfs as CSV files
  write.csv(agreement.df, paste0(main.dir, "Code/Agreement_Statistics/", substr(basename(sources(classified.images[[image]])),1,8), "_agreement_reference_data.csv"), row.names = FALSE)

  
  ## Agreement across regional reference data 
  ## Create dataframes to populate with results 
  region.df.supp <- as.data.frame(matrix(nrow = (length(unique(regions$Name)) + 1) * 4, ncol = 6))
  colnames(region.df.supp) <- c("REGION","SPECIES","AREA","SENSITIVITY","SPECIFICITY","BAL_AGREEMENT")
  region.df.supp$REGION <- rep(c(unique(regions$Name),"ALL_VIMS_CCRM_data"), each = 4)
  region.df.supp$SPECIES <- rep(levels(agreement.tibble$PRED_CLASS_CONF)[1:4], times = length(unique(regions$Name)) + 1)
  region.df.main <- as.data.frame(matrix(nrow = (length(unique(regions$Name)) + 1), ncol = 5))
  colnames(region.df.main) <- c("REGION","AREA","SENSITIVITY","SPECIFICITY","BAL_AGREEMENT")
  region.df.main$REGION <- c(unique(regions$Name),"ALL_VIMS_CCRM_data")
  ## Loop through each region 
  for(region in unique(regions$Name)){
    
    ## Crop to the region of interest 
    ref.region <- subset(ref, REGION == region)
    classified.image.region <- crop(classified.image, ref.region, mask = TRUE)
    agreement.tibble.region <- tidyterra::as_tibble(classified.image.region) %>% filter(if_all(c(CLASS,PRED_CLASS_CONF), complete.cases))
    agreement.tibble.region$PRED_CLASS_CONF <- factor(agreement.tibble.region$PRED_CLASS_CONF, levels = c("Spartina_alterniflora","Phragmites_australis","Juncus_roemerianus","Spartina_patens","NoVeg"))
    agreement.tibble.region$CLASS <- factor(agreement.tibble.region$CLASS, levels = c("Spartina_alterniflora","Phragmites_australis","Juncus_roemerianus","Spartina_patens","NoVeg"))
    
    ## Populate dataframe with results -- region.df.supp
    region.df.supp$AREA[which(region.df.supp$REGION == region)] <- sum(st_area(ref.region)) / 1000000
    supp.metrics <- agreement.tibble.region %>% group_by(CLASS) %>% multiclass.metrics(truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")
    region.df.supp$SENSITIVITY[which(region.df.supp$REGION == region & region.df.supp$SPECIES %in% as.character(unique(supp.metrics$CLASS)))] <- round(subset(supp.metrics, .metric == "sensitivity")$.estimate * 100)
    region.df.supp$SPECIFICITY[which(region.df.supp$REGION == region & region.df.supp$SPECIES %in% as.character(unique(supp.metrics$CLASS)))] <- round(subset(supp.metrics, .metric == "specificity")$.estimate * 100)
    region.df.supp$BAL_AGREEMENT[which(region.df.supp$REGION == region & region.df.supp$SPECIES %in% as.character(unique(supp.metrics$CLASS)))] <- round(subset(supp.metrics, .metric == "bal_accuracy")$.estimate * 100)

    ## Populate dataframe with results -- region.df.main
    region.df.main$AREA[which(region.df.main$REGION == region)] <- sum(st_area(ref.region)) / 1000000
    region.df.main$SENSITIVITY[which(region.df.main$REGION == region)] <- round(sensitivity(agreement.tibble.region, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100)
    region.df.main$SPECIFICITY[which(region.df.main$REGION == region)] <- round(specificity(agreement.tibble.region, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100)
    region.df.main$BAL_AGREEMENT[which(region.df.main$REGION == region)] <- round(bal_accuracy(agreement.tibble.region, truth = CLASS, estimate = PRED_CLASS_CONF, estimator = "micro")$.estimate * 100)
    
  }
  
  ## Populate rows for ALL_VIMS_CCRM_data
  ## region.df.supp
  region.df.supp$AREA[which(region.df.supp$REGION == "ALL_VIMS_CCRM_data")] <- sum(st_area(ref)) / 1000000
  region.df.supp$SENSITIVITY[which(region.df.supp$REGION == "ALL_VIMS_CCRM_data")] <- agreement.df$CLASS_SENSITIVITY[1:4]
  region.df.supp$SPECIFICITY[which(region.df.supp$REGION == "ALL_VIMS_CCRM_data")] <- agreement.df$CLASS_SPECIFICITY[1:4]
  region.df.supp$BAL_AGREEMENT[which(region.df.supp$REGION == "ALL_VIMS_CCRM_data")] <- agreement.df$CLASS_BAL_AGREEMENT[1:4]
  ## region.df.main
  region.df.main$AREA[which(region.df.main$REGION == "ALL_VIMS_CCRM_data")] <- sum(st_area(ref)) / 1000000
  region.df.main$SENSITIVITY[which(region.df.main$REGION == "ALL_VIMS_CCRM_data")] <- unique(agreement.df$OVERALL_SENSITIVITY)
  region.df.main$SPECIFICITY[which(region.df.main$REGION == "ALL_VIMS_CCRM_data")] <- unique(agreement.df$OVERALL_SPECIFICITY)
  region.df.main$BAL_AGREEMENT[which(region.df.main$REGION == "ALL_VIMS_CCRM_data")] <- unique(agreement.df$OVERALL_BAL_AGREEMENT)
  
  ## Export dataframes 
  write.csv(region.df.supp, paste0(main.dir, "Code/Agreement_Statistics/", substr(basename(sources(classified.images[[image]])),1,8), "_regional_agreement_supplemental.csv"), row.names = FALSE)
  write.csv(region.df.main, paste0(main.dir, "Code/Agreement_Statistics/", substr(basename(sources(classified.images[[image]])),1,8), "_regional_agreement_main.csv"), row.names = FALSE)
  
}


