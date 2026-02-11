# Scientific Methodology

This document describes the scientific methodology used in the Chesapeake Bay Marsh Vegetation Mapping workflow.

## Table of Contents

1. [Input Data](#input-data)
2. [Vegetation Indices](#vegetation-indices)
3. [Training Data Selection](#training-data-selection)
4. [Random Forest Model](#random-forest-model)
5. [Classification and Certainty Thresholding](#classification-and-certainty-thresholding)
6. [Agreement Assessment](#agreement-assessment)
7. [Annual Assessments](#annual-assessments)

---

## Input Data

### Reference data

#### 1. Local, species-specific reference data
VIMS (Virginia Institute of Marine Science) CCRM (Center for Coastal Resources Management) conducted aerial surveys in May 2021, which were used to generate a vegetation classification shapefile delineating *S. alterniflora*, *P. australis*, *S. patens*, *J. roemerianus*, and no vegetation, where no vegetation primarily consists of water and tidal flats. This data was provided by project co-authors and is not publicly available. Data that is obtained for this step should be a shapefile (.shp and its required support files). This shapefile will be read in using the `sf` package as a sf object with geometry type MULTIPOLYGON and have the following feature structure:

| VEG_CLASS             | DATE     | geometry        |
|-----------------------|----------|-----------------|
| Spartina_alterniflora | May_2021 | ...             |
| Spartina_patens       | May_2021 | ...             |
| ...                   | ...      | ...             |

#### 2. Broad, more generalized reference data
VIMS CCRM Shoreline and Tidal Marsh Inventory represents Virginia's most comprehensive assessment of coastal marsh and serves as a benchmark for shoreline management. The most recent inventory at the time of article publication was released in 2019 and was generated from 2006-2017 aerial imagery, although the updated Tidal Marsh Inventory is planned for release in 2026. Data are publicly available and can be requested for the state of Virginia through the [VIMS CCRM Data Request Form](https://www.vims.edu/ccrm/research/inventory/virginia/). 

### Regional shapefiles 

#### 1. NOAA Habitat Focus Area boundary  
Virginia's Middle Peninsula is designated by NOAA as a [Habitat Focus Area](https://www.habitatblueprint.noaa.gov/habitat-focus-areas/middle-peninsula-virginia/), a priority region where restoration is targeted to have the greatest impact, increase coastal resiliency, and advance science and conservation. This data is not available online, but this respository includes the required input shapefile in the **Middle_Peninsula_HFA**(`2_Annual_Assessments/Data/Input_Data/Middle_Peninsula_HFA/`) folder.

#### 2. Surface salinities 
Surface salinities are used to subset the Middle Peninsula Habitat Focus Area to terrestrial areas bordering water with salinities above 15 ppt, a threshold selected given salinity preferences of the marsh species considered. Salinity zones were defined using 1985-2018 mean surface salinities from the [Chesapeake Bay Program](https://www.chesapeakebay.net/what/publications/chesapeake-bay-mean-surface-salinity-1985-2018). This data should be downloaded, unzipped, placed in the **Mean_Surface_Salinity** (`2_Annual_Assessments/Data/Input_Data/Mean_Surface_Salinity/`) folder. 

### Image classification data 

#### 1. Planet multispectral imagery  
The classification uses Planet SuperDove Level 3B surface reflectance imagery with 8 spectral bands:

| Band | Name | Wavelength (nm) |
|------|------|-----------------|
| 1 | Coastal Blue | 431-452 |
| 2 | Blue | 465-515 |
| 3 | Green I | 513-549 |
| 4 | Green | 547-583 |
| 5 | Yellow | 600-620 |
| 6 | Red | 650-680 |
| 7 | Red Edge | 697-713 |
| 8 | NIR | 845-885 |

Cloud and cloud shadow masking is performed using the UDM2 (Usable Data Mask 2) quality band, where pixels flagged by Planet as non-clear are set to NA.
The UDM2 is included in Planet SuperDove Level 3B surface reflectance imagery. Unfortunately, Planet imagery cannot be released with the publication, but there are avenues for Planet data access, including:
1. [Planet's Education and Research Program](https://www.planet.com/industries/education-and-research/)
2. [NASA's Commercial Satellite Data Acquisition Program](https://www.earthdata.nasa.gov/about/csda)

Once data access is obtained, [Planet Explorer](https://www.planet.com/explorer/) can be used for requesting and downloading data. This code is designed to work with the ZIP files provided through Planet Explorer; other data formats will require edits to the R code.

#### 2. USGS DEM

Digital Eleveation Model (DEM) elevation data is obtained from the USGS National Elevation Dataset (NED) at 1/3 arc-second resolution (~9 m). The DEM is resampled to match the Planet imagery resolution using bilinear interpolation. No external data is required. 

---

## Vegetation Indices

The following spectral indices are computed from the Planet imagery bands:

### 1. Normalized Difference Vegetation Index (NDVI)
```
NDVI = (NIR - Red) / (NIR + Red)
```
Standard vegetation index sensitive to chlorophyll content and vegetation density.
**Reference:** Rouse, J.W., et al. (1974). Monitoring vegetation systems in the Great Plains with ERTS.

### 2. Enhanced Vegetation Index (EVI)
```
EVI = 2.5 × (NIR - Red) / (NIR + 6 × Red - 7.5 × Blue + 1)
```
Optimized vegetation index that reduces atmospheric and soil background influences.
**Reference:** Huete, A., et al. (2002). Overview of the radiometric and biophysical performance of the MODIS vegetation indices. Remote Sensing of Environment.

### 3. Normalized Difference Water Index (NDWI)

```
NDWI = (Green - NIR) / (Green + NIR)
```
Sensitive to water content in vegetation and standing water.

**Reference:** McFeeters, S.K. (1996). The use of the Normalized Difference Water Index (NDWI) in the delineation of open water features.

### 4. Soil Adjusted Vegetation Index (SAVI)

```
SAVI = ((NIR - Red) × 1.5) / (NIR + Red + 0.5)
```
Minimizes soil brightness influences using a soil adjustment factor (L = 0.5).

**Reference:** Huete, A.R. (1988). A soil-adjusted vegetation index (SAVI). Remote Sensing of Environment.

### 5. Wide Dynamic Range Vegetation Index (WDRVI)

```
WDRVI = (0.2 × NIR - Red) / (0.2 × NIR + Red)
```
Enhances dynamic range in high biomass regions where NDVI saturates.

**Reference:** Gitelson, A.A. (2004). Wide Dynamic Range Vegetation Index for remote quantification of biophysical characteristics of vegetation.

### 6. Band Ratios

```
BG   = Blue / Green
GR   = Green / Red
NIRR = NIR / Red
```
Simple band ratios that capture spectral relationships useful for vegetation discrimination.

---

## Training Data Selection

### Sampling Strategy

Training data is selected using stratified random sampling:

1. **Class balance**: The number of training pixels per class is set to the percentage of the smallest class to ensure balanced representation set by `TRAINING_FRACTION` in `config.R` or `config_local.R` (the accompanying publication uses a percentage of 80%).

2. **Spatial sampling**: Pixels are randomly sampled from within each vegetation class polygon.

3. **Per-image sampling**: Training data is extracted separately for each satellite acquisition date to account for temporal variations.

### Data Preprocessing

The `tidymodels` recipe applies the following preprocessing steps:

1. **Correlation filtering** (`step_corr`): Removes predictor variables with Spearman correlation > `CORRELATION_THRESHOLD` set in `config.R` or `config_local.R` (the accompanying publication uses a threshold of 0.95) to reduce multicollinearity.

2. **Centering** (`step_center`): Centers predictors to a mean of zero.

3. **Scaling** (`step_scale`): Scales predictors to a standard deviation of one.

---

## Random Forest Model

### Algorithm

Random forest is an ensemble learning method that constructs multiple decision trees during training and outputs the mode (classification) of individual tree predictions.

### Hyperparameter Tuning

Hyperparameters are optimized using 5-fold cross-validation with grid search:

| Parameter | Description | Search Range |
|-----------|-------------|--------------|
| `trees` | Number of trees in the forest | Tuned via grid search |
| `mtry` | Number of variables randomly sampled at each split | Tuned via grid search |
| `min_n` | Minimum number of observations in terminal nodes | Tuned via grid search |

The optimal model is selected based on the highest ROC-AUC (Area Under the Receiver Operating Characteristic Curve).

### Variable Importance

Variable importance is computed as the mean decrease in node impurity (Gini importance) across all trees in the forest.

---

## Classification and Certainty Thresholding

### Probability Predictions

For each pixel, the random forest model outputs class probabilities for all five vegetation classes. The predicted class is the one with the highest probability.

### Certainty Threshold

A certainty threshold of set by `CERTAINTY_THRESHOLD` in `config.R` or `config_local.R` is applied to filter low-certainty predictions(the accompanying publication uses a threshold of 0.7), although this is adjustable to fit stakeholder needs:

- Pixels with maximum class probability >= `CERTAINTY_THRESHOLD` (defaults to 0.7) are assigned the predicted class
- Pixels with maximum class probability < `CERTAINTY_THRESHOLD` (defaults to 0.7) are flagged as "low certainty"

---

## Accuracy Assessment

### Metrics

Classification agreement is assessed using the following metrics:

| Metric | Description |
|--------|-------------|
| **Sensitivity** | True positive rate: TP / (TP + FN) |
| **Specificity** | True negative rate: TN / (TN + FP) |
| **Balanced Agreement** | Average of sensitivity and specificity |
| **Overall Agreement** | (TP + TN) / Total |

Balanced agreement is recommended for imbalanced datasets. For comparison with other studies, the more commonly used overall agreement is also reported, despite being less statistically representative for imbalanced datasets. Agreement is evaluated across all reference data and for each vegetation class individually.

## Annual Assessments

Annual assessments can be assessed for multiple years and for multiple images per year. As described in the publication, model development can be generated using multiple satellite images and the image that produces the highest balanced agreement with reference data is selected for application on larger-scale imagery. Larger-scale imagery then applies the selected random forest classification model to each image. Data are subset to the input parameters defined by the random forest model, removing those with high collinearity.
If there are multiple images collected per year, their median probabilities per vegetation class are computed and the final classification assigns each pixel to the class with the highest median probabilities across satellite images from the same year.

---

## References

1. Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5-32.

2. Gitelson, A.A. (2004). Wide Dynamic Range Vegetation Index for remote quantification of biophysical characteristics of vegetation. Journal of Plant Physiology, 161(2), 165-173.

3. Huete, A., et al. (2002). Overview of the radiometric and biophysical performance of the MODIS vegetation indices. Remote Sensing of Environment, 83(1-2), 195-213.

4. Huete, A.R. (1988). A soil-adjusted vegetation index (SAVI). Remote Sensing of Environment, 25(3), 295-309.

5. McFeeters, S.K. (1996). The use of the Normalized Difference Water Index (NDWI) in the delineation of open water features. International Journal of Remote Sensing, 17(7), 1425-1432.

6. Rouse, J.W., et al. (1974). Monitoring vegetation systems in the Great Plains with ERTS. NASA Special Publication, 351, 309-317.
