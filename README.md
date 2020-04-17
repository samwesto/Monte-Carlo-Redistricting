# Monte Carlo for Redistricting
There have been numerous attempts to apply Monte Carlo techniques to quanitfy the process of political redistricting in the United States. Due to the need for extensive optimisation and the lack of a 
one size fits all approach, previous methods do not translate well to the UK. In the following project I have implemented a number of popular techniques to this end, making them applicable to the structure of the UK electoral system. 
More than this, I have created new algorithms that perform better than those existing when applied to UK redistricting problems. 

## Getting started 

### Prerequisites 
To run these algorithms, the following packages are required: 
numpy, fiona, random, geopandas, csv, shapely, itertools, scipy, cvxpy

### Required Data
The data input type is modelled around what is available from the Office of National Statistics. Any data source can be used however, lookup tables and data sets from the ONS will be already be in the correct format.

#### Required filed
Main Dataframe - df
Type: Pandas DataFrame
Columns (must be in the following order):  unique ward code, ward population, 

OAframe
Type: Dataframe
Columns (in the following order): unique OA code, population, ward that contains OA area

Adjacency lists:
Type: csv files with each row containing a list of ward (OA) index numbers for df (OAframe), such that each column is a separate number. 
ward neighbours,* OAs**. 

Shape file
Type: .shp
Column: a column named ‘geometry’ with ward shape data in, must be shapely polygons.


*Neighbours = all wards that are adjacent to the given ward, measured by touching edges. If you do not have this information, use function neighbours in MC. 
** OAs = list of OA areas contained inside this ward. If not available use argument OAs = None. You will then be restricted to the use of Area compactness only.

