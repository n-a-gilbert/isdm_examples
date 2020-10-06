# Data & code for examples of the 3 models in Gilbert et al., "Integrating harvest and camera trap data in species distribution models"

### NOTE: Using NIMBLE

These models are written in NIMBLE, a language derived from BUGS. Unlike previous programs (WinBUGS, JAGS), NIMBLE models are programmable objects in R, but are compiled in C++ for speed. Therefore, prior to installing NIMBLE, you must have Rtools installed so the code can be compiled. Please see the [NIMBLE website](https://r-nimble.org/download) for instructions.

### Focal species: bobcat (Lynx rufus)

We provide example data and code for bobcat. In Wisconsin, bobcats are known to be associated with areas of greater forest cover and to avoid urban areas--hence, we use canopy cover and impervious cover within 5x5 km grid cells to predict bobcat occurrence. Camera data comes from 6/1/2018-10/14/2018, a period following parturition but before the harvest season (20 October 2018-31 January 2019).

### Data

The data included in this repository consists of three objects:

1. constants : constants (e.g. indices) for the model. In list form.
2. data : data for the model. In list form.
3. cov : two covariates, in a polygon shapefile

More detail about each is as follows:

#### constants

* ncell : the number of grid cells (5x5 km)
* nsurveys: the number of week-long sampling occasions the camera was active for
* cell : the identity of grid cells containing a camera
* ncams : the number of cameras
* ncounty: the number of counties
* low: a lower bounding index used to define which grid cells fall within which county
* high: an upper bounding index used to define which grid cells fall within which county
* neigh: the number of adjacencies between grid cells 

#### data

* y : a matrix with rows representing cameras and columns representing sampling occasions. 1 = bobcat detected during sampling occasion, 0 = bobcat NOT detected during sampling occasion, NA = camera not active during this sampling occasion
* e: an effort bias term used in the harvest submodel. This is, for each county, the percent of hunters statewide who operated in a given county
* w: the number of bobcats harvested in each county in 2018
* num: number of adjacent grid cells to each grid cell
* adj: the identies of adjacent grid cells for each grid cell
* weights: how neighboring grid cells are weighed. All = 1
* forest: mean % canopy cover within each 5x5 km grid cell
* imperv: mean % impervious cover within each 5x5 km grid cell
* cam_can: % canopy over of the 30x30m cell containing the camera coordinates
* yday: matrix with ordinal date of the beginning of each survey occasion for each camera
* yday2: matrix with ordinal date^2

#### cov

Shapefile of covariates (vectorized grid).

