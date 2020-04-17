# Monte Carlo for Redistricting
There have been numerous attempts to apply Monte Carlo techniques to quanitfy the process of political redistricting in the United States. Due to the need for extensive optimisation and the lack of a 
one size fits all approach, previous methods do not translate well to the UK. In the following project, I have implemented a number of popular techniques to this end, making them applicable to the structure of the UK electoral system. 
Building upon this, I have created new algorithms that perform better under our system than existing models.

## Getting started 

### Prerequisites 
To run these algorithms, the following packages are required: 
numpy, fiona, random, geopandas, csv, shapely, itertools, scipy, cvxpy

### Required Data
The data input type is modelled around what is available from the Office of National Statistics. Any data source can be used however, lookup tables and data sets from the ONS will be already be in the correct format. Within the folder data, there is an example for each one of these require data types, pertaining to the wards of the West Midlands. 


#### Main Dataframe - df
Type: Pandas DataFrame
Columns (must be in the following order):  unique ward code, ward population

#### Shape file - shapefile
Type: .shp
Column: a column named ‘geometry’ with ward shape data in, must be type shapely polygons. The data must be in the same order as the data in df.

#### Ordinance Survey Area data - OAframe (optional: only for use with MMI compactness)
Type: Dataframe
Columns (in the following order): unique OA code, population, ward code that contains OA area

#### Ward Adjacency lists - adjlist
Type: csv 
File with each row containing a list of df index numbers corresponding to wards that are adjacent geographically. The file should be structured such that each row features information about a single ward and each column contains a index. This file can be created using the function neighbours in MCMC_functions.py, if provided with the ward shape file. 

#### OA ward relationships - OAcol
type: csv 
File with each row containing a list of unique OA codes corresponding to OAs that are contained inside a specific ward (corresponding to that row). The file should be structured such that each row features information about a single ward and each column contains a unique OA code. This file can be created using the function OAcol in MCMC_functions.py, if provided with the OAdataframe. 


* If data relating to OA areas (OAcol and OAframe) is not available, or the user does not want to use MMI compactness, these files can be omitted. All MCMC algorithms will not require these files by default. 

## Running the algorithms

