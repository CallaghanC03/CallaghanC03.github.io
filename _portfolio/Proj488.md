---
title: "GES488 Final Project"
excerpt: "I looked at the effect of high amounts of vegetation on the access to cable internet. I did this using satellite data from Landsat imagery. This was then used to create an NDVI vegetation index. This compares the amount of Red light to Near Infrared light because plants relect a lot on NIR but very little Red light. Then the averages of NDVI in all tracts were then compared to percentages of household cable internet returning a strong correlation. Below is the clustering of high to high(red) and low to low(blue) values between the two. <br/><img src='/images/Moran.png'>"
collection: portfolio
---

# Final Project GES 488

```{r Library_Loader, warning = FALSE}
library(sf) # 
library(tidyverse)
library(tidycensus)
library(purrr)
library(dplyr) 
library(tmap)
library(ggplot2)
library(tmaptools)
library(leaflet)
library(corrplot) 
source("http://www.sthda.com/upload/rquery_cormat.r")
library(ppcor) 
library(olsrr) 
library(lares)
library(jtools)

census_api_key("Enter API KEY") # Enter API "KEY"

```

### Data Pull

Using the tidyverse/tidycensus library with the `get_acs()` function to pull data for tracts across DC and Maryland of boardband internet access

```{r Data_Grab, warning = FALSE}
census_api_key('a91b19846dc2e7972c640b84abc39f74fef6e5b8')

vars <- load_variables(2022, "acs5/subject")  
vars$groupid <- lares::left(vars$name, n=9)
vars_filter <- filter(vars, groupid %in% c('S2801_C01'))

tract <- get_acs(
  geography = 'tract',
  variables = c(internet_cnt = "S2801_C01_012", # Total Households with Internet
                sat_net = 'S2801_C01_018', # Sattelite Internet Estimates
                dial_up = 'S2801_C01_013', # Dial Up Internet Estimates
                no_net = 'S2801_C01_019', # No Internet Esstimates
                cable_net = "S2801_C01_017", # Cable/Fiber Optic Internet Estimates
                hh_total ="S2801_C01_001" # Total number of Households
                ),
  state = c('District of Columbia', 'Maryland'),
  year = 2022,
  output = "wide",
  geometry = TRUE
) %>% filter(!st_is_empty(.)) %>% na.omit() # REMOVES EMPTY GEOMETRY TRACTS AND NA ENTRIES
tract <- separate(tract, 
                   NAME, 
                   into = c("tract", "countyname", "state"),
                   sep = ", "
                   )

```

### Data Frame Edits

Prepares Data for further Manipulation and creates new fields to be used for statistics. Normalizes data across tracts to avoid outliers from population density. Tracts and data then written to shapefile

```{r Data_Frame_Edits, warning = FALSE}

st_transform(tract, crs = 3857)
tract$dial_up_pct <- tract$dial_upE/tract$hh_totalE # % of households that have Dial Up internet
tract$cable_pct <- tract$cable_netE/tract$hh_totalE # % of households that have cable(fiber-optic or other)
tract$no_net_pct <- tract$no_netE/tract$hh_totalE
tract$sat_net_pct <- tract$sat_netE/tract$hh_totalE

# st_write(tract, "I:/Final_Proj/data/tract.shp")

```

#### Landsat Data

Landsat Data gotten from link and code below:

##### <https://code.earthengine.google.com/b010c2b896602c6e95745fd2d1b9b79f> -\>

##### Landsat Collection 2 Level- 2 Surface Reflectance Science Product courtesy of the U.S. Geological Survey

```{r LANDSAT_GEE_Code}
// ** Function to apply mask and add quality bands to Landsat 8 images
var addQualityBands = function(image, bands) {
  return maskClouds(image)
    //  convert reflectance numbers to 0 to 1 (NASA spec document conversion)
    .multiply(0.0000275).add(-0.2)
    // time in days
    .addBands(image.metadata('system:time_start'))
     // NDVI
    .addBands(image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI'));
};


// ** Function to mask clouds in imagery
var maskClouds = function(image) {
  var cloudDilBitMask = 1 << 1;  // this code to make a binary number
  var cirrusBitMask = 1 << 2; 
  var cloudsBitMask = 1 << 3; 
  var cloudShadowBitMask = 1 << 4;   
  var qa = image.select('QA_PIXEL');  //  in this QA band, binary numbers represent mask codes
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
    .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
    .and(qa.bitwiseAnd(cloudDilBitMask).eq(0))
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));  
  return image.updateMask(mask);   //this applies the mask above to the image.
}



// --> Add the Landsat surface reflectance collection
// Filter by date
var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate('2017-06-01', '2017-09-30').filterBounds(MDDCVA)

// Calculate median for each band using filtered Collection
// Map the cloud screening function over the collection
var collection_m = landsat.map(addQualityBands);


// Create a greenest pixel (NDVI) composite
var NDVIPixelComposite = collection_m.qualityMosaic('NDVI');
print(NDVIPixelComposite)
Map.addLayer(NDVIPixelComposite, {bands: ['NDVI'], min: 0, max: 0.5, palette: ['red', 'yellow','green']}, 'Greenest pixel (NDVI) composite');
Map.addLayer(NDVIPixelComposite, {bands: ['SR_B4', 'SR_B5', 'SR_B3'], min: 0, max: 0.5}, '(RGB)Greenest pixel (NDVI) composite');

//Exports data to Google Drive
Export.image.toDrive({
  image: NDVIPixelComposite.float(),
  description: 'NDVI',
  scale: 30,
  region: MDDCVA,
  folder: 'GEE',
  maxPixels: 395764910,
});

```

### Basic Regression

```{r Regression}
tract_fnl <- st_read("I:/Final_Proj/Arc/tract_final2.shp")
variables <- c('cbl_pct','st_nt_p', 'dl_p_pc', 'NDVI_mean')
proj <- tract_fnl[,variables]
str(proj)
summary(proj)


corr_subset <-st_drop_geometry(proj)
var_proj <- var(corr_subset, use = "complete.obs") # variance matrix
cov_proj <-cov(corr_subset, use = "complete.obs") # covariance matrix
cor_proj <- cor(corr_subset, use = "complete.obs") # correlation matrix

var_proj
cov_proj
cor_proj

rquery.cormat(corr_subset)



mod.1 <- lm(cbl_pct ~ NDVI_mean, data = proj)
summary(mod.1)
summ(mod.1, confint = TRUE, digits = 4)

```
