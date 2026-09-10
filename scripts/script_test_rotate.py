import xarray as xr
import numpy as np
from dempy.data import open_dataset
from dempy.rotate import Transformer, get_rotated_dem, get_target_grid


dem_o = (open_dataset("dem_alp")
       .isel(band=0).band_data
       .isel(x=slice(1500, 2000), y=slice(3000,3500))
       )





T = Transformer(dem_o)

dem = dem_o.pipe(T.center_dem)

theta = np.pi/8


target_grid = get_target_grid(dem_o)
demrot = get_rotated_dem(dem, target_grid, theta)
demback = get_rotated_dem(demrot, dem, -1*theta)


demrot.plot(); plt.show()
demback.plot(); plt.show()


demback
