# Species distribution models affected by positional uncertainty in species occurrences can still be ecologically interpretable

This repository was created as a supporting material for the article exploring the ecological interpretability of species distribution models affected by positional uncertainty in occurrence data. 

Previous work has documented the effect of positional uncertainty on model predictive performance, but its consequences for inference about species-environment relationships remains largely unknown.  In our study, we explored the extent to which parameter estimation is affected by positional uncertainty. In total, we simulated 12 560 combinations of virtual and real environmental data and virtual species (Workflows 1, 2) to investigate our assumptions and have fitted over 628 000 models. Additionally, we tested our assumptions using real environmental variables and real species (Band-tailed Pigeon; Workflow 3).

**Relevant paper:**
[Link for published paper](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.06358)

Code authors: [Lukas Gabor](https://scholar.google.cz/citations?user=pLQXY5wAAAAJ&hl=cs),
              [Vojtech Bartak](https://scholar.google.cz/citations?user=p8WAo8oAAAAJ&hl=cs&oi=ao),
              [Matthew S. Rogan](https://scholar.google.cz/citations?user=OTgp4V8AAAAJ&hl=cs&oi=ao)
              

Date: 10/14/2022

In this repository we provide used R scripts that can be used to recreate our study.

**Environmental and species data**

***Workflow 1 - virtual environmental data + virtual species***

As both, environemtal and species data have been simulated there is no need to download any. Just follow attached R script. 

***Workflow 2 - real environmental data + virtual species***

Environmental data are available for download on [Download Center of Autonomous body National Center for Geographic Information web site](https://centrodedescargas.cnig.es/CentroDescargas/locale?request_locale=en#) (amount of forest, amount of grassland, elevation). Aspect and topography wetness index were derived from the elevation.

***Workflow 3 - real environmental data + real species***

Band-tailed Pigeon dataset can be downloaded [here](https://www.sciencebase.gov/catalog/item/5eb4485182ce25b5135abeea). Environmental data used for Band-tailed Pigeon SDMs are available for download on [modis](https://modis.gsfc.nasa.gov/) (Mean winter EVI, ), [chelsea](https://chelsa-climate.org/) (Seasonality of precipitation, Mean annual temperature, ), [earthenv](https://www.earthenv.org/) (Cloud cover, Terrain ruggedness, EVI spatial heterogeneity), [soilgrids](https://www.soilgrids.org/) (Soil clay content, Soil silt content), [envidat](https://www.envidat.ch/#/) (Growing season precipitation).

**File Description:**

R_script_workflow1: Step-by-step guideline which allow replication of Workflow 1 from generating virtual data to plotting results. 

R_script_workflow2: Step-by-step guideline which allow replication of Workflow 2 from generating virtual species to plotting results. 

R_script_workflow3: Step-by-step guideline which allow replication of Workflow 3 from simulating positional error in occurrence data to plotting results.

Cite the code: [![DOI](https://zenodo.org/badge/470330827.svg)](https://zenodo.org/badge/latestdoi/470330827)
