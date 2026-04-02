# README
S. Filhol, March 2026

This is a fresh start of `dempy` that had become too combersome to track. `dempy` intents to bring together a collection of tools to work with Digital Elevation Models. `dempy` started some years ago under Python 2.7 and could use some refreshing and cleanup.

Many of the current GIS library minimize memory usage by making use of raster files for many steps along any processing pipeline. When working on local projects, computers have nowadays sufficient memory to handle on the fly many processing steps. Also, dempy tries using JIT compiled code running in parallel computing using Numba, and efficient IO storage using Zarr stores.  

The intent here is to convert to a cleaner toolbox, integrated around xarray dataset and prevent using raster files completely. This is particularly interesting to work from memory when working on small projects.

## TODO:
- [ ] create a main class pulling the various method together
- [ ] avoid gdal at almost all cost
- [ ] make the method compatible with xarray dataset
- [ ] write a documentation 
- [ ] integrate the compilation of Diamond_square with the package install
- [ ] keep adding new functionalities :)
