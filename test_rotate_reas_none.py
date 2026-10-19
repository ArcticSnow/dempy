import xarray as xr
import numpy as np
from dempy.data import open_dataset
from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
from dempy.terrain import hillshade
import matplotlib.pyplot as plt 


dem = (open_dataset("dem_alp")
       .isel(band=0).band_data
       .isel(x=slice(1000, 2500), y=slice(2000,3500))
       )

dem.plot(); plt.show()

# Compute hillshade with default azimuth (315°) and altitude (45°)
shade = hillshade(dem)

shade.plot(cmap="gray"); plt.show()

# dem
# <xarray.DataArray 'band_data' (y: 1500, x: 1500)> Size: 9MB
# [2250000 values with dtype=float32]
# Coordinates:
#   * y            (y) float64 12kB 6.541e+06 6.541e+06 ... 6.496e+06 6.496e+06
#   * x            (x) float64 12kB 8.911e+05 8.911e+05 ... 9.36e+05 9.36e+05
#     band         int64 8B 1
#     spatial_ref  int64 8B ...


# deepseek 4 pro res effort none hallucinates a non existing  dempy.terrain
