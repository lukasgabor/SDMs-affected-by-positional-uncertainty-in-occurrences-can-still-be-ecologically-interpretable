# Species distribution models affected by positional uncertainty in species occurrences can still be ecologically interpretable

This repository was created as a supporting material for the article exploring the ecological interpretability of species distribution models negatively affected by positional uncertainty in occurrence data. 

Previous work has documented the effect of positional uncertainty on model predictive performance, but its consequences for inference about species-environment relationships remains largely unknown.  In our study, we use simulated data and two example species to investigate this issue. Specifically, we artificially applied positional error to species occurrence data and assessed the effects on variable importance and coefficient estimates from Maximum Entropy-based SDMs. 

In this repository we provide used R scripts and "vs" object that can be used to recreate our study.

Band-tailed Pigeon dataset can be downloaded [here](https://www.sciencebase.gov/catalog/item/5eb4485182ce25b5135abeea). Environmental data used for Band-tailed Pigeon SDMs are available for download on [modis](www.modis.gsfc.nasa.gov) (Mean summer EVI, Mean winter EVI, ), [chelsea](www.chelsa-climate.org) (Seasonality of precipitation, Mean annual temperature, ), [earthenv](www.earthenv.org) (Cloud cover, Terrain ruggedness, EVI spatial heterogeneity), [soilgrids](www.soilgrids.org) (Soil clay content, Soil silt content), [envidat](www.envidat.ch) (Growing season precipitation). 

Cantabrian Brown bear dataset is available upon request from the authors and all used environmental data are available for download on [Download Center of Autonomous body National Center for Geographic Information web site](https://centrodedescargas.cnig.es/CentroDescargas/locale?request_locale=en#) (Elevation, Amount of forest, Amount of grassland). Topography wetness index and Aspect were derived from the Elevation model.

Code author: [Lukas Gabor](https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs)

Date: 10/14/2022

**Relevant paper:**
Here I will add the paper citation

**File Description:**

virtual.species: Saved "vs" object that can be used for further analysis.

R_script_virtual: Step-by-step guideline for generating virtual species, sampling occurrence data, simulating positional error, developing SDMs and plotting results. 

R_script_bear: Step-by-step guideline for simulating positional error in bear´s data, developing SDMs and plotting results. The same script was also used for Band-tailed pigeon.
