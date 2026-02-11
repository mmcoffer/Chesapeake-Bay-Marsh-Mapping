# Chesapeake Bay Marsh Vegetation Mapping

Random forest classification of marsh vegetation types using Planet SuperDove 8-band multispectral imagery and USGS DEM data for the Chesapeake Bay region.

## Overview

This repository contains R scripts for developing and applying a machine learning model to classify marsh vegetation species from satellite imagery. The workflow supports both initial model development using reference data and annual assessments for monitoring marsh vegetation changes over time. Much of the input data is not publicly available; see below and [METHODOLOGY.md](METHODOLOGY.md) for options to obtain these datasets and for details on data access for those datasets that are publicly available. 

## Target Vegetation Classes

- ***Spartina alterniflora*** (synonmously *Sporobolus alterniflorus*; saltmarsh cordgrass)
- ***Phragmites australis*** (common reed)
- ***Spartina patens*** (saltmeadow hay)
- ***Juncus roemerianus*** (black needlerush)
- **No vegetation** (primarily water and tidal flats)

## Repository Structure

This repository includes all necessary script files and the folder structure for required input data (most of which is not publicly available). The code will create the remaining folder structure as needed. This code is designed to run on one or multiple PlanetScope images for both model development and annual assessments. If multiple images are used for model development, the random forest model with the highest balanced agreement is used for annual assessments. 

```
├── config.R                                 # Configuration template (copy to config_local.R)
├── requirements.R                           # Package installation script
├── METHODOLOGY.md                           # Scientific methodology documentation
├── LICENSE                                  # MIT License
│
├── 1_Model_Development/                     # Create random forest classification model 
│   ├── 1a_Format_Data.R                     # Format satellite and DEM data
│   ├── 1b_Create_Random_Forest.R            # Create training data and train random forest
│   ├── 1c_Apply_Classification.R            # Apply model to reference data
│   ├── 1d_Assess_Agreement.R                # Validate model agreement
│   └── Data/                                # Location for input and output data folders 
|       └── Input_Data/                      # Subfolder containing data used to generate model
|           ├── VIMS_Reference_Data/         # Place formatted reference data here; not publicly available
|           └── Planet_Multispectral_Data/   # Subfolder containing Planet imagery 
|               └── 1_Level3B_Data_Zipped/   # Place zipped data from Planet Explorer here; not publicly available
│
└── 2_Annual_Assessments/                    # Generate large-scale annual classifications 
    ├── 2a_Format_Data.R                     # Format data for NOAA Habitat Focus Area
    ├── 2b_Apply_Classification.R            # Apply trained model to new imagery
    └── Data/                                # Location for input and output data folders 
        └── Input_Data/                      # Subfolder containing data used to generate model
            ├── VIMS_Tidal_Marsh_Inventory/  # Place VIMS CCRM TMI data here; publicly available
            └── Mean_Surface_Salinity/       # Place unzipped Chesapeake Bay Program surface salinity data here; publicly available
            └── Middle_Peninsula_HFA         # Includes shapefile of Middle Peninsula HFA boundary; not publicly available 
            └── Planet_Multispectral_Data/   # Subfolder containing Planet imagery
                └── 1_Level3B_Data_Zipped/   # Place zipped data from Planet Explorer here; not publicly available 
```

## Requirements

### R Version
- R >= 4.0.0

### Required Packages
Run `requirements.R` to install all dependencies:
```r
source("requirements.R")
```

Packages used:
- `terra` - Raster data processing
- `sf` - Vector data processing
- `tidymodels` - Machine learning framework
- `tidyterra` - Integration of terra with tidyverse
- `FedData` - Download USGS elevation data
- `randomForest` - Random forest implementation
- `tools` - File path utilities
- `alphahull` - Used to simplify complex shapefiles

## Configuration

1. Copy the configuration template:
   ```r
   file.copy("config.R", "config_local.R")
   ```

2. Edit `config_local.R` to set your local paths:
   ```r
   BASE_DIR <- "/path/to/your/project/directory"
   ```

3. The scripts will automatically use `config_local.R` if it exists, otherwise fall back to `config.R`.

## Input Data Requirements

### Reference data

   #### 1. Local, species-specific reference data
   - VIMS (Virginia Institute of Marine Science) CCRM (Center for Coastal Resources Management) aerial vegetation classification shapefile, generated in May 2021 and delineating the five classes listed above
   - **Not publicly available**

   #### 2. Broad, more generalized reference data
   - VIMS CCRM Shoreline and Tidal Marsh Inventory shapefiles for annual assessments
   - **Publicly available** for the state of Virginia through the [VIMS CCRM Data Request Form](https://www.vims.edu/ccrm/research/inventory/virginia/)

### Regional shapefiles 

   #### 1. NOAA Habitat Focus Area boundary 
   - [NOAA Middle Peninsula Habitat Focus area](https://www.habitatblueprint.noaa.gov/habitat-focus-areas/middle-peninsula-virginia/) boundary 
   - **Not publicly available,** but included in this repository

   #### 2. Surface salinities 
   - Chesapeake Bay Program mean surface salinities from 1985-2018 
   - **Publicly available** from the [Chesapeake Bay Program](https://www.chesapeakebay.net/what/publications/chesapeake-bay-mean-surface-salinity-1985-2018)

### Image classification data 

   #### 1. Planet multispectral imagery 
   - Level 3B SuperDove 8-band multispectral imagery
   - Bands: coastal_blue, blue, green_i, green, yellow, red, rededge, nir
   - Associated UDM2 (Usable Data Mask) files for cloud masking
   - **Not publicly available,** but see [METHODOLOGY.md](METHODOLOGY.md) for data access options 

   #### 2. USGS DEM
   - National Elevation Dataset (NED) 1/3 arc-second (~9 m resolution)
   - **Publicly available** and automatically downloaded via the `FedData` package

## Usage

### Random Forest Model Development

Execute scripts in order:

```r
# 1. Format input data
source("1_Model_Development/1a_Format_Data.R")

# 2. Create training data and train random forest model
source("1_Model_Development/1b_Create_Random_Forest.R")

# 3. Apply random forest model to reference data
source("1_Model_Development/1c_Apply_Classification.R")

# 4. Assess classification agreement
source("1_Model_Development/1d_Assess_Agreement.R")
```

### Annual Assessments

After model development is complete:

```r
# 1. Format new imagery
source("2_Annual_Assessments/2a_Format_Data.R")

# 2. Apply classification
source("2_Annual_Assessments/2b_Apply_Classification.R")
```

## Output Data

### Model Development Outputs

Located in `1_Model_Development/`:

- **Classified rasters** (`Data/Output_Data/Classified_Data/`): GeoTIFF files with vegetation class predictions and per-pixel probabilities
- **Trained models** (`Data/Output_Data/Random_Forest_Model/`): Random forest model (.RData) and data preprocessing object (.RData)
- **Training data** (`Data/Output_Data/Training_Data/`): CSV files and shapefiles of sampled training pixels
- **Variable importance** (`Data/Output_Data/Variable_Importance/`): CSV files with random forest feature importance rankings
- **Agreement statistics** (`Data/Output_Data/Agreement_Statistics/`): CSV files with agreement metrics (sensitivity, specificity, balanced agreement)

### Annual Assessment Outputs

Located in `2_Annual_Assessments/Data/Output_Data/`:

- **Classified rasters**: GeoTIFF files with vegetation class predictions and per-pixel probabilities, with multiple images collected in the same year aggregated via the median probability for each individual date

## Methodology

See [METHODOLOGY.md](METHODOLOGY.md) for detailed documentation of:
- Input Data
- Vegetation Indices
- Training Data Selection
- Random Forest Model
- Classification and Certainty Thresholding
- Agreement Assessment
- Annual Assessments 

## Citation

If you use this code in your work, please cite:

> Coffer, M.M., Trinh, R., Mitchell, M., Angstadt, K., Stanhope, D., Lv, Z., Nunez, K., Bartlett, N.D., Wiltsie, D., Sullivan, S., & Schaeffer, B.A. (2026). Tidal marsh species mapping using commercial satellite imagery for enhanced coastal management in Chesapeake Bay. *Remote Sensing Applications: Society and Environment*, 101902. https://doi.org/10.1016/j.rsase.2026.101902

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
