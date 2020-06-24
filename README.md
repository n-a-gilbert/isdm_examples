# Data & code for examples of the 3 models in Gilbert et al., "Integrating harvest and camera trap data to improve species distribution models"

### NOTE: Using NIMBLE

These models are written in NIMBLE, a language derived from BUGS. Unlike previous programs (WinBUGS, JAGS), NIMBLE models are programmable objects in R, but are compiled in C++ for speed. Therefore, prior to installing NIMBLE, you must have Rtools installed so the code can be compiled. Please see the [NIMBLE website](https://r-nimble.org/download) for instructions.

### Data

The data included in this repository consists of four objects:

1. constants : constants (e.g. indices) for the model. In list form.
2. data : data for the model. In list form.
3. cov : two covariates, in a polygon shapefile
4. test : a dataframe containing testing data to evaluate predictive performance

More detail about each is as follows:

#### constants

* ncell : the number of grid cells (8.5x8.5 km)
* cell : the identity of grid cells containing a camera
* nsite : the number of cameras
* ncounty: the number of counties
* low: a lower bounding index used to define which grid cells fall within which county
* high: an upper bounding index used to define which grid cells fall within which county
* k: the number of adjacencies between grid cells 

#### data

* y : the number of week-long sampling occasions bear was detected at a given camera
* nsuryveys: the number of week-long sampling occasions the camera was active for
* e: an effort offset used in the harvest model. This is, for each county, the number of harvest authorizations per square kilometer of the bear management zone in which the county falls
* w: the number of bears harvested in each county in 2018
* num: number of adjacent grid cells to each grid cell
* adj: the identies of adjacent grid cells for each grid cell
* weights: how neighboring grid cells are weighed. All = 1
* forest: mean % canopy cover within each 8.5x8.5 km grid cell
* imperv: mean % impervious cover within each 8.5x8.5 km grid cell

#### cov

Shapefile of covariates (vectorized grid).

#### test

Camera data witheld to evaluate predictive performance of model via AUC. A dataframe with two columns:

* INDEX: the identity of the grid cell that the camera falls inside of
* Y: Indicator of whether (1) or not (0) bear was ever detected at the camera during 2018

