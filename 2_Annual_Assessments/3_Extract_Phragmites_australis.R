

## Extract satellite classification in areas of Phragmites australis 


## Load required packages 
require(terra)
require(sf)


## Set main working directory 
main.dir <- "C:/Users/mcoffer/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/COASTAL_WATER/Planet_Marsh_Mapping/"


## Read in input shapefiles and subset to polygons with presence Phragmites australis 
tmi.15ppt <- st_transform(st_read(paste0(main.dir,"2_Annual_Assessments/Data/Input_Data/VA_TMI_2011_2019_HFA/TMI_15ppt_edited.shp")), crs = "EPSG:32618")
tmi.15ppt.phrag <- subset(tmi.15ppt, !(PrcntPhrag == "undetermined"))
## Format PrcntPhrag data 
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "<1")]     <- 0
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "1-25")]   <- 13
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "25-50")]  <- 37.5
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "75-85")]  <- 80
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "75-100")] <- 87.5
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "1-30")]   <- 15.5
tmi.15ppt.phrag$PrcntPhrag[which(tmi.15ppt.phrag$PrcntPhrag == "50-75")]  <- 62.5
tmi.15ppt.phrag$PrcntPhrag <- as.numeric(tmi.15ppt.phrag$PrcntPhrag)


## Read in output classification for the year of interest and its classification key 
sat <- rast(paste0(main.dir,"2_Annual_Assessments/Data/Output_Data/2021_classified_image.tif"))
sat.crop <- crop(sat, tmi.15ppt.phrag, mask = TRUE, touches = FALSE)
sat.key <- read.table(paste0(main.dir,"2_Annual_Assessments/Data/Output_Data/2021_classification_key.txt"), header = FALSE, skip = 2)
## Create a dataframe to populate with results 
year.df <- as.data.frame(matrix(nrow = nrow(tmi.15ppt.phrag), ncol = 8))
colnames(year.df) <- c("MarshNo","REF_cont_pct","REF_cont_binned","SAT_total_n","SAT_lowconf_n","SAT_lowconf_pct","SAT_Phrag_n","SAT_Phrag_pct")
year.df$MarshNo <- tmi.15ppt.phrag$MarshNo
## Loop through each polygon in tmi.15ppt.phrag and extract satellite-derived results 
for(phrag.polygon in 1:nrow(tmi.15ppt.phrag)){
  
  ## Extract sat.crop within phrag.polygon  
  sat.phrag.polygon <- extract(sat.crop, tmi.15ppt.phrag[phrag.polygon,], touches = FALSE)
  sat.phrag.polygon$PRED_CLASS <- as.character(sat.phrag.polygon$PRED_CLASS)
  sat.phrag.polygon$PRED_CLASS[which(sat.phrag.polygon$MAX_CONFIDENCE < 0.7)] <- "Low confidence"
  sat.phrag.polygon$PRED_CLASS[which(sat.phrag.polygon$PRED_CLASS == sat.key$V1[which(sat.key$V2 == "MED_phragmites_australis")])] <- "MED_phragmites_australis"
  ## Populate year.df with results -- reference data 
  year.df$REF_binary_pct[phrag.polygon]  <- ifelse(tmi.15ppt.phrag$PhragPres[phrag.polygon] == "<50% Phragmites australis", tmi.15ppt.phrag$PhragPres, ">=50% Phragmites australis")
  year.df$REF_cont_pct[phrag.polygon]    <- tmi.15ppt.phrag$PrcntPhrag[phrag.polygon]
  year.df$REF_cont_binned[phrag.polygon] <- ifelse(year.df$REF_cont_pct[phrag.polygon] < 75, ifelse(year.df$REF_cont_pct[phrag.polygon] <= 25, "Low_percentage", "Mid_percentage"), "High_percentage")
  ## Populate year.df with results -- satellite data
  year.df$SAT_total_n[phrag.polygon]     <- length(which(!(is.na(sat.phrag.polygon$PRED_CLASS))))
  year.df$SAT_lowconf_n[phrag.polygon]   <- length(which(sat.phrag.polygon$PRED_CLASS == "Low confidence"))
  year.df$SAT_lowconf_pct[phrag.polygon] <- round((year.df$SAT_lowconf_n[phrag.polygon] / year.df$SAT_total_n[phrag.polygon]) * 100, digits = 2)
  year.df$SAT_Phrag_n[phrag.polygon]     <- length(which(sat.phrag.polygon$PRED_CLASS == "MED_phragmites_australis"))
  if(nrow(sat.phrag.polygon) > 10 & year.df$SAT_lowconf_pct[phrag.polygon] < 50){
    year.df$SAT_Phrag_pct[phrag.polygon] <- round((year.df$SAT_Phrag_n[phrag.polygon]/(year.df$SAT_total_n[phrag.polygon] - year.df$SAT_lowconf_n[phrag.polygon])) * 100, digits = 2)
  }
}
  
## Export CSV file 
write.csv(year.df, paste0(main.dir, "Manuscript/Figures/Figure_5_MWU_Boxplots/Extracted_Phragmites_australis/2021_extracted_Phragmites_australis.csv"), row.names = F)

