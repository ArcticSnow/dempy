import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from .rotate import get_rotated_dem, get_target_grid

class HorizonComputer:

    def __init__(self, dem):
        self.dem = dem
        self.y_spacing = (dem.y.isel(y=0) - dem.y.isel(y=1)).item()
        self.diffs_dem = []
        self.horizons = []
        self.steps = []

    def step(self, step_size):
        shift = int(np.array(self.steps).sum() + step_size)
        if len(self.steps) == 0: # first step
            diff = (self.dem.shift(y=-shift) - self.dem)/(self.y_spacing*shift)
            self.horizons.append(diff.where(~diff.isnull(),0))
            self.diffs_dem.append(diff)
            self.steps.append(step_size)
            return
        diff = (self.dem.shift(y=-1*shift) - self.dem)/(self.y_spacing*shift)
        self.horizons.append(diff.where(diff > self.horizons[-1], self.horizons[-1]))
        self.diffs_dem.append(diff)
        self.steps.append(step_size)




class ShadowIterator():

    def __init__(self, dem, sun_elevation, strategy="max"):
        self.sun_elevation  = sun_elevation
        self.strategy = strategy
        self.dem = dem
        self.diff_dem = dem.shift({'y':1}, fill_value=None) - dem
        self.shadows = xr.full_like(dem, False, dtype=bool)
        self.peaks_to_explore = (
                xr.full_like(dem, True, dtype=bool)
                .where(self.diff_dem/30 <= -1*sun_elevation, False)
                )
        self.steps = []


    def iterate(self):
        if self.strategy == "max":
            max_y_idx = self.dem.where(self.peaks_to_explore).argmax(dim="y").item()
        else:
            max_y_idx = np.nonzero(self.peaks_to_explore.values)[0][-1]
        max_elevation = self.dem.isel(y=max_y_idx).item()
        max_y = self.dem.y.isel(y=max_y_idx).item()
        top_solar_ray = max_elevation-self.sun_elevation*(self.dem.y - max_y)
        shadows = ((self.dem < top_solar_ray) & (self.dem.y >= max_y))
        self.shadows |= shadows
        self.peaks_to_explore &= ~shadows
        self.peaks_to_explore[max_y_idx] = False
        # step related parameters
        self.steps.append(
                dict(
                    max_elevation= max_elevation,
                    max_y_idx=max_y_idx,
                    max_y=max_y,
                    shadows=shadows,
                    top_solar_ray=top_solar_ray,
                    explored=max_y_idx
                    )

                )
        return self


    def compute(self):
        try:
            while True:
                self.iterate()
        except ValueError, IndexError:
            print("iteration finished")
            return


class VectorisedShadowIterator:

    def __init__(self, dem, sun_elevation, name_x="x", name_y="y"):
        """
        aucune diff avec la version non vectorisée!
        sun_elevation: elevation du soleil en radians
        """
        self.sun_elevation  = sun_elevation
        self.dem = dem.rename({name_x: "x", name_y: "y"})
        self.size_x = self.dem.sizes["x"]
        self.size_y = self.dem.sizes["y"]
        self.resolution = (self.dem.x[1] - self.dem.x[0]).item()
        self.diff_dem = self.dem.shift({'y':1}, fill_value=None) - self.dem
        self.shadows = xr.full_like(self.dem, False, dtype=bool)
        self.sun_slope = np.tan(self.sun_elevation)
        # self.current_iteration = 0
        # self.steps = []


    @property
    def peaks_to_explore(self):
        return (
            xr
            .full_like(self.dem, True, dtype=bool)
            .where(self.diff_dem/self.resolution <= -1*np.tan(self.sun_elevation), False)
        )

    # @property
    # def current_y(self):
    #     return self.size_y - self.current_iteration - 1

    def compute(self, idx_y):
         idx_x = np.nonzero(self.peaks_to_explore.isel(y=idx_y).values)[0]
         elevations_explored = self.dem[dict(y=idx_y, x=idx_x)]
         y_explored = elevations_explored.y.item()
         top_solar_ray = elevations_explored - self.sun_slope*(self.dem.y - y_explored)
         new_shadows = ((self.dem <= top_solar_ray) & (self.dem.y >= y_explored))
         self.idx_x = idx_x
         self.y_explored = y_explored
         self.top_solar_ray = top_solar_ray
         self.new_shadows = new_shadows


    # def iterate(self):
    #     all_idx_x =  np.nonzero(self.peaks_to_explore.isel(y=self.current_y).values)[0]
    #     elevations_explored = self.dem[dict(y=self.current_y, x=all_idx_x)]
    #     y_explored = elevations_explored.y.item()
    #     sun_slope = np.tan(self.sun_elevation)
    #     top_solar_ray = elevations_explored - sun_slope*(self.dem.y - y_explored)
    #     # (-np.tan(0.15)*(demtest.y - demtest.y[-1])+last_z).plot(ax=_ax1)
    #     self.shadows[dict(x=all_idx_x)] = (
    #             self.shadows[dict(x=all_idx_x)] 
    #             | ((self.dem <= top_solar_ray) & (self.dem.y >= y_explored))
    #              )
    #     self.peaks_to_explore &= ~self.shadows
    #     self.peaks_to_explore[dict(y=self.current_y, x=all_idx_x)] = False
    #     # documenting the step
    #     return self
    #
    #
    # def compute(self):
    #     for y in range(self.size_y):
    #         self.current_iteration = y
    #         self.iterate()


class ShadowerStep:

    def __init__(self, dem, sun_elevation, name_x="x", name_y="y"):
        self.sun_elevation  = sun_elevation # elevation solaire en radians
        self.dem = dem.rename({name_x: "x", name_y: "y"})
        self.size_x = self.dem.sizes["x"]
        self.size_y = self.dem.sizes["y"]
        self.resolution = (self.dem.x[1] - self.dem.x[0]).item()
        self.diff_dem = self.dem.shift({'y':1}, fill_value=None) - self.dem
        # self.diff_dem = self.dem - self.dem.shift({'y':-1}, fill_value=None)
        self.shadows = xr.full_like(self.dem, False, dtype=bool)
        self.peaks_to_explore = (
                xr.full_like(self.dem, True, dtype=bool)
                .where(self.diff_dem/self.resolution <= -1*sun_elevation, False)
                )
        self.current_iteration = 0
        self.steps = []

# j0, i0 = 500, 350
# lat0 = lat.isel(y=j0,x=i0).item()
# lon0 = lon.isel(y=j0,x=i0).item()
# R = 13000
# # distance de la sphère (terrestre) plan tangent au point choisi
# h_tangent = 1000*R*(1-np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)-np.sin(lat0)*np.sin(lat))
# h_tangent.plot()

from pyproj import Transformer as ProjTransformer

def compute_lon_lat_from_lambert_xy(dem):
    """
    """
    # Get 1‑D coordinate vectors
    x_vals = dem.x.values          # shape (x)
    y_vals = dem.y.values          # shape (y)
    # Create 2‑D grid matching the DEM layout
    X, Y = np.meshgrid(x_vals, y_vals)   # both shape (y, x)

    # Set up transformer: Lambert‑93 (EPSG:2154) → WGS‑84 (EPSG:4326)
    lambert_to_wgs84 = ProjTransformer.from_crs(
        crs_from="EPSG:2154",
        crs_to="EPSG:4326",
        always_xy=True,
    )
    # 4️⃣ Transform the whole grid at once
    lon_vals, lat_vals = lambert_to_wgs84.transform(X, Y)   # both (y, x)
    # 5️⃣ Wrap as DataArrays, preserving DEM dimensions & coordinates
    lon = xr.DataArray(
        lon_vals,
        dims=dem.dims,
        coords=dem.coords,
        name="longitude",
    )
    lat = xr.DataArray(
        lat_vals,
        dims=dem.dims,
        coords=dem.coords,
        name="latitude",
    )
    return lon, lat



def compute_shadows(dem, sun_azimuth, sun_elevation):
    x_center = dem.x.mean().item()
    y_center = dem.y.mean().item()
    dem_centered = dem.copy()
    dem_centered.coords['x'] = dem.x - x_center
    dem_centered.coords['y'] = dem.y - y_center
    target_grid = get_target_grid(dem_centered)
    dem_rotation_angle = sun_azimuth - np.pi
    demrot = get_rotated_dem(dem_centered, target_grid, dem_rotation_angle)
    SI = VectorisedShadowIterator(demrot, sun_elevation)
    SI.compute()
    shadowsrot = SI.shadows
    shadows_centered = get_rotated_dem(SI.shadows, dem_centered, -1*dem_rotation_angle)
    shadows = shadows_centered.copy()
    shadows.coords['x'] = shadows.x + x_center
    shadows.coords['y'] = shadows.y + y_center
    return shadows

