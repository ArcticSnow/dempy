import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
from dempy.fast_shadows import compute_shadows, compute_lon_lat_from_lambert_xy


dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"

dem = xr.open_dataset(dem_path).isel(band=0).band_data.isel(x=slice(None, None, 10), y=slice(None, None, 10))


lon, lat = compute_lon_lat_from_lambert_xy(dem)

# Attach as coordinates to the DEM for easy indexing
dem = dem.assign_coords(lon=lon, lat=lat)

shadows = compute_shadows(dem, 2, 0.05)


