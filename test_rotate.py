import xarray as xr
import numpy as np
from dempy.data import open_dataset
from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
import matplotlib.pyplot as plt 


dem = (open_dataset("dem_alp")
       .isel(band=0).band_data
       .isel(x=slice(1000, 2500), y=slice(2000,3500))
       )

dem.plot(); plt.show()

# dem
# <xarray.DataArray 'band_data' (y: 1500, x: 1500)> Size: 9MB
# [2250000 values with dtype=float32]
# Coordinates:
#   * y            (y) float64 12kB 6.541e+06 6.541e+06 ... 6.496e+06 6.496e+06
#   * x            (x) float64 12kB 8.911e+05 8.911e+05 ... 9.36e+05 9.36e+05
#     band         int64 8B 1
#     spatial_ref  int64 8B ...

# Hillshade parameters
azimuth = 315.0   # light direction in degrees
altitude = 45.0   # light height in degrees

# Terrain derivatives, using the projected x/y coordinates for true spacing
dy, dx = np.gradient(dem.values, dem.y.values, dem.x.values)

slope = np.arctan(np.hypot(dx, dy))
aspect = np.arctan2(-dy, dx)

az = np.deg2rad(azimuth)
alt = np.deg2rad(altitude)

hillshade = 255 * (
    np.sin(alt) * np.cos(slope)
    + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
)

hillshade = xr.DataArray(
    hillshade, dims=dem.dims, coords=dem.coords
).clip(0, 255)

hillshade.plot(cmap="Greys_r"); plt.show()




