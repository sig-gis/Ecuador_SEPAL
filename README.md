# Ecuador_SEPAL
------
This repository contains the preprocessing code for landsat-8 and and sentinel-2 to create biweekly composites. Both scripts are based upon the Google Earth engine (GEE). 


## sentinel-2
------
Sentinel-2 data are available in Google Earth Engine in the “Sentinel-2 MSI: MultiSpectral Instrument, Level-1C” image collection (ID: “COPERNICUS/S2”). These data are available as top-of-atmosphere radiometrically-corrected imagery. 

Processing is done through:

**Cloud Masking:** The QA60 band is used for cloud removal and a custom build cloud masking algorithm is used to remove clouds from the images.

**Cloud shadow masking:** The Temporal Dark-Outlier Mask algorithm (TDOM) approachis used to remove cloud shadows. The algorithm identifies cloud shadows as pixels that are both dark in the infrared bands and that are dark statistical outliers across the time series.

**BRDF correction:** bidirectional reflectance distribution function is applied to for the correction of view and illumination angle effects.

**illumination correction:** correction for topographic effects on illumation.

**Atmospheric correction:** Atmospheric corresion is done using the 6s emulator. The 6S emulator is an open-source atmospheric correction tool. Follow the steps on here (https://github.com/samsammurphy/6S_emulator ) to install the 6s emulator. 

## Landsat-8
------
Landsat 8 data used to create the cloud free annual composites or biweekly mosaics are available in Google Earth Engine in the “USGS Landsat 8 Surface Reflectance Tier 1” image collection (ID: “LANDSAT/LC08/C01/T1_SR”). These data have been atmospherically corrected from top-of-atmosphere to surface reflectance using the Landsat 8 Surface Reflectance Code (LaSRC), produced by the U.S. Geological Survey (USGS). Surface Reflectance Landsat data is then corrected for topographic, radiometric, and atmospheric distortion; and clouds and cloud-shadows are masked.

**Cloud Masking:** The QA60 band is used for cloud removal and a custom build cloud masking algorithm is used to remove clouds from the images.

**Cloud shadow masking:** The Temporal Dark-Outlier Mask algorithm (TDOM) approachis used to remove cloud shadows. The algorithm identifies cloud shadows as pixels that are both dark in the infrared bands and that are dark statistical outliers across the time series.

**BRDF correction:** bidirectional reflectance distribution function is applied to for the correction of view and illumination angle effects.

**illumination correction:** correction for topographic effects on illumation.
