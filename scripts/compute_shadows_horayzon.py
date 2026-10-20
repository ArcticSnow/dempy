import xarray as xr
import numpy as np
from dempy.shadow import DEMShadowMapper
import matplotlib.pyplot as plt
from datetime import datetime
from skyfield.api import utc

import matplotlib
from dempy.fast_shadows import HorizonComputer, compute_lon_lat_from_lambert_xy

matplotlib.use('module://matplotlib-backend-kitty')


dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"
dem = xr.open_dataset(dem_path).band_data

# dem.plot()
# plt.show()

# dem_4326 = dem.rio.reproject(4326)
# dem_4326.to_netcdf("dem_4326.nc")

dem_4326 = xr.open_dataset("dem_4326.nc")
# dem_4326.band_data.plot(); plt.show()


# Initialize shadow Mapper
mapper = DEMShadowMapper(
    dem_file="dem_4326.nc",
    domain={"lon_min": 6.95, "lon_max": 7.0, "lat_min": 45.915, "lat_max": 45.935},
    dist_search=0.5, # [km]
)
# date = datetime(2019,3,22,10,32,45, tzinfo=utc)
# mapper.get_solarposition(date)
sw_dir_cor = mapper.compute_sw_dir_cor_map_from_position(np.pi/3, np.pi + np.deg2rad(2.7), 1e6)
# sw_dir_cor2 = mapper.compute_sw_dir_cor_map_from_position(np.pi/4, np.pi, 1e9)
sw_reproj_match = sw_dir_cor.rename(dict(lon="x", lat="y")).transpose('y', 'x').rio.reproject_match(dem)

# (sw_dir_cor - sw_dir_cor2).plot(); plt.show()
# sw_dir_cor2.plot(); plt.show()

shadows1 = (sw_reproj_match > 0)
# (sw_reproj_match > 0).plot(); plt.show()

horizon = HorizonComputer(dem)
for _step_size in 100*[1]:# + 20*[2] + 20*[5] + 10*[20]:
    horizon.step(_step_size)
shadows2 = (np.atan(horizon.horizons[-1]) < np.pi/3)


np.atan(horizon.horizons[-1]).plot(); plt.show()

(shadows2.astype('int') - shadows1.astype('int')).plot(); plt.show()

shadows2.plot(); plt.show()
shadows1.plot(); plt.show()
