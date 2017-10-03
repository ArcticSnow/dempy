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
- Generating random DEM
- Derive morphometric data
- plotting DEM




## Requirements

Requires to install third party libraries for using certain functions:

    openCV
    gdal
    pyproj
    numpy
    matplotlib
    scipy
    cv2
    laspy
    
    pyfftw