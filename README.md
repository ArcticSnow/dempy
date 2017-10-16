# Dempy: a python package to work with GIS data types

Simon Filhol, October 2017

## Installation

```sh
pip install dempy
```

WARNING: Use at your own risk. 



## Objective

This package is a collection of useful tool to work with geospatial type of data within Python. Tools for working with rasters, shapefile, pointcloud data, and more advanced dem analysis tool such as fourier decomposition. 



Some of the tools included are:

- geoTIFF import
- DEM detrend with linear plane
- DEM filtering (FFT or Convolution)
- Generating random DEM with diamondSquare (implemented in Cython)
- Derive morphometric data
- plotting DEM
- derive Blowing snow terrain parameter: 
  - Tabler's drift profile (Tabler, R.D., 1975. Predicting Profiles of Snowdrifts in Topographic Catchments. In *Proceedings of the 43rd Annual Western Snow Conference*. San Diego, CA, pp. 87–97.)
  - Winstral Sx (Winstral, A., Elder, K. & Davis, R.E., 2002. Spatial Snow Modeling of Wind-Redistributed Snow Using Terrain-Based Parameters. *Journal of Hydrometeorology*, 3(5), pp.524–538.)




## Requirements

Requires a few libraries such as laspy, gdal, cython, geopandas ...
