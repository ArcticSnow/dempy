import xarray as xr
import numpy as np
from dempy.data import open_dataset
from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource


dem = (open_dataset("dem_alp")
       .isel(band=0).band_data
       .isel(x=slice(1000, 2500), y=slice(2000,3500))
       )

# dem
# <xarray.DataArray 'band_data' (y: 1500, x: 1500)> Size: 9MB
# [2250000 values with dtype=float32]
# Coordinates:
#   * y            (y) float64 12kB 6.541e+06 6.541e+06 ... 6.496e+06 6.496e+06
#   * x            (x) float64 12kB 8.911e+05 8.911e+05 ... 9.36e+05 9.36e+05
#     band         int64 8B 1
#     spatial_ref  int64 8B ...

dem.plot(); plt.show()

# Create a light source (azimuth in degrees, altitude in degrees)
ls = LightSource(azdeg=315, altdeg=45)

# Compute hillshade
hillshade = ls.hillshade(dem.values, vert_exag=3)

# Plot the hillshade
plt.figure(figsize=(10, 8))
plt.imshow(hillshade, cmap='gray', origin='upper')
plt.colorbar(label='Hillshade intensity')
plt.title('Hillshade')
plt.show()


# Compute the slope (gradient) components in easting and northing
import numpy as np

dzdx = np.gradient(dem.values, dem.x.values, axis=1)   # easting slope (d z/d x)
dzdy = np.gradient(dem.values, dem.y.values, axis=0)   # northing slope (d z/d y)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

im0 = axes[0].imshow(dzdx, cmap='RdBu_r', origin='upper')
axes[0].set_title('Slope (easting) d z/d x')
fig.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(dzdy, cmap='RdBu_r', origin='upper')
axes[1].set_title('Slope (northing) d z/d y')
fig.colorbar(im1, ax=axes[1])

plt.tight_layout()
plt.show()



