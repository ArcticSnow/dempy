"""
Python class to compute horizon angle. 
S. Filhol, July 2026

These classes build on the code examples from HORAYZON by C. Steger. 

TODO:
- [x] Implement horizon for Curved DEM
- [x] Implement horizon for Planar DEM
- [ ] Write a little example code 

"""
#%%
#========= Librairies imports ========
import pdb
import numpy as np
import xarray as xr
import horayzon as hray
from pathlib import Path


#%%
#========= Class definition ========

class PlanarDEMHorizon:
    def __init__(self, dem_file, domain, n_azimuth=360, dist_search=500):
        """
        Initialize class.

        Args:
            dem_file (str or Path): Path to DEM projected in WGS84 (ESPG:4326)
            domain (dict): X/Y (East/North) bounding box of area fo interest
            dist_search (float): buffer distance around the area of interest to use for shadow computation. Unit is the same as the DEM spatial unit

        TODO:
        - [ ] test and debug
        - [ ] add possibility to derive domain automatically from dem itself and dist_search

        """
        self.file_dem = Path(dem_file)
        self.path_out = self.file_dem.parent
        self.domain = domain            # same unit as DEM
        self.dist_search = dist_search  # same unit as DEM
        self.n_azimuth = n_azimuth

        self.slice_in = None
        self.hori = None
        self.n_azimuth = None
        self.vec_north = None
        self.vec_norm = None
        self.vert_grid = None

        self._load_dem()



    def _load_dem(self):

        # Check domain is within extent of the DEM with the buffer (dist_search)
        if ((self.domain["x_min"] >= self.domain["x_max"])
                or (self.domain["y_min"] >= self.domain["y_max"])):
            raise ValueError("Invalid self.domain specification")

        # Compute outer self.domain boundaries
        self.domain_outer = {"x_min": self.domain["x_min"] - self.dist_search,
                        "x_max": self.domain["x_max"] + self.dist_search,
                        "y_min": self.domain["y_min"] - self.dist_search,
                        "y_max": self.domain["y_max"] + self.dist_search}

        self.x, self.y, self.elevation = hray.load_dem.planar_dem(self.file_dem, self.domain_outer)

        # Compute indices of inner domain
        self.slice_in = (slice(np.where(self.y >= self.domain["y_max"])[0][-1],
                          np.where(self.y <= self.domain["y_min"])[0][0] + 1),
                    slice(np.where(self.x <= self.domain["x_min"])[0][-1],
                          np.where(self.x >= self.domain["x_max"])[0][0] + 1))
        print("Inner domain size: " + str(self.elevation[self.slice_in].shape))

        # Create directional unit vectors for inner domain
        offset_0 = self.slice_in[0].start
        offset_1 = self.slice_in[1].start
        dem_dim_0, dem_dim_1 = self.elevation.shape
        vec_norm = np.zeros((dem_dim_0 - (2 * offset_0),
                             dem_dim_1 - (2 * offset_1), 3), dtype=np.float32)
        vec_norm[:, :, 2] = 1.0
        vec_north = np.zeros((dem_dim_0 - (2 * offset_0),
                              dem_dim_1 - (2 * offset_1), 3), dtype=np.float32)
        vec_north[:, :, 1] = 1.0
        #pdb.set_trace()
        # Merge vertex coordinates and pad geometry buffer
        self.vert_grid = hray.auxiliary.rearrange_pad_buffer(*np.meshgrid(self.x, self.y), self.elevation)


    def compute_horizon(self, n_azimuth=None, return_ds=False):
        if n_azimuth is None:
            azim_num = self.n_azimuth
        else:
            azim_num = n_azimuth

        offset_0 = self.slice_in[0].start
        offset_1 = self.slice_in[1].start
        dem_dim_0, dem_dim_1 = self.elevation.shape

        # Compute horizon
        self.hori, self.azimuth = hray.horizon.horizon_gridded(
            self.vert_grid, dem_dim_0, dem_dim_1,
            self.vec_norm, self.vec_north,
            offset_0, offset_1,
            dist_search=self.dist_search,
            azim_num=azim_num)

        if return_ds:

            ds = xr.Dataset(
                coords=dict(
                    azim=(["azim"], self.azimuth, {"units": "radian"}),
                    y=(["y"], self.y[self.slice_in[0]], {"units": "meter"}),
                    x=(["x"], self.x[self.slice_in[1]], {"units": "meter"}),
                ),
                data_vars=dict(
                    horizon=(["azim", "y", "x"], np.moveaxis(self.hori, 2, 0),
                               {"units": "radian"})
                )
            )
            return ds


    def compute_slope_aspect_svf(self):
        if self.hori is None:
            self.compute_horizon()

        # Compute slope
        x_2d, y_2d = np.meshgrid(self.x, self.y)
        slice_in_a1 = (slice(self.slice_in[0].start - 1, self.slice_in[0].stop + 1),
                       slice(self.slice_in[1].start - 1, self.slice_in[1].stop + 1))
        vec_tilt = hray.topo_param.slope_plane_meth(x_2d[slice_in_a1],
                                                    y_2d[slice_in_a1],
                                                    self.elevation[slice_in_a1])[1:-1, 1:-1]
        del x_2d, y_2d

        # Compute Sky View Factor
        svf = hray.topo_param.sky_view_factor(self.azimuth, self.hori, vec_tilt)

        # Compute slope angle and aspect
        slope = np.arccos(vec_tilt[:, :, 2].clip(max=1.0))
        aspect = np.pi / 2.0 - np.arctan2(vec_tilt[:, :, 1],
                                          vec_tilt[:, :, 0])
        aspect[aspect < 0.0] += np.pi * 2.0  # [0.0, 2.0 * np.pi]

        ds = xr.Dataset(
                coords=dict(
                    y=(["y"], self.y[self.slice_in[0]], {"units": "meter"}),
                    x=(["x"], self.x[self.slice_in[1]], {"units": "meter"}),
                ),
                data_vars=dict(
                    elevation=(["y", "x"], self.elevation[self.slice_in],
                               {"xg_name": "ellipsoidal height", "units": "m"}),
                    slope=(["y", "x"], slope,
                           {"long_name": "slope angle", "units": "radian"}),
                    aspect=(["y", "x"], aspect,
                            {"long_name": "slope aspect (clockwise from North)",
                            "units": "radian"}),
                    svf=(["y", "x"], svf,
                         {"long_name": "sky view factor", "units": "-"}),
                )
            )

        return ds





class CurvedDEMHorizon:
    def __init__(self, dem_file, domain, n_azimuth=360, dist_search=0.5):
        """
        Initialize class.

        Args:
            dem_file (str or Path): Path to DEM projected in WGS84 (ESPG:4326)
            domain (dict): lat/lon bounding box of area fo interest
            dist_search (float): buffer distance around the area of interest to use for shadow computation

        TODO:
        - [ ] test and debug

        """
        self.file_dem = Path(dem_file)
        self.path_out = self.file_dem.parent
        self.domain = domain
        self.dist_search = dist_search
        self.ellps = "WGS84"
        self.n_azimuth = n_azimuth

        self.trans_ecef2enu = None
        self.slice_in = None
        self.vec_tilt_enu = None
        self.hori = None
        self.n_azimuth = None
        self.vec_north_enu = None
        self.vec_norm_enu = None

        self._load_dem()


    def _load_dem(self):

        self.domain_outer = hray.domain.curved_grid(self.domain, self.dist_search, self.ellps)
        self.lon, self.lat, self.elevation = hray.load_dem.wgs84_dem(self.file_dem, self.domain_outer)

        # Compute indices of inner domain
        self.slice_in = (slice(np.where(self.lat >= self.domain["lat_max"])[0][-1],
                          np.where(self.lat <= self.domain["lat_min"])[0][0] + 1),
                    slice(np.where(self.lon <= self.domain["lon_min"])[0][-1],
                          np.where(self.lon >= self.domain["lon_max"])[0][0] + 1))
        print("Inner domain size: " + str(self.elevation[self.slice_in].shape))

        # Compute ellipsoidal heights
        #elevation += hray.geoid.undulation(lon, lat, geoid="EGM96")  # [m]

        # Compute ECEF coordinates
        x_ecef, y_ecef, z_ecef = hray.transform.lonlat2ecef(*np.meshgrid(self.lon, self.lat), self.elevation, ellps=self.ellps)
        dem_dim_0, dem_dim_1 = self.elevation.shape

        # Compute ENU coordinates
        self.trans_ecef2enu = hray.transform.TransformerEcef2enu(lon_or=self.lon[int(len(self.lon) / 2)], lat_or=self.lat[int(len(self.lat) / 2)], ellps=self.ellps)
        x_enu, y_enu, z_enu = hray.transform.ecef2enu(x_ecef, y_ecef, z_ecef,
                                                      self.trans_ecef2enu)

        # Compute unit vectors (up and north) in ENU coordinates for inner domain
        vec_norm_ecef = hray.direction.surf_norm(*np.meshgrid(self.lon[self.slice_in[1]],
                                                              self.lat[self.slice_in[0]]))
        vec_north_ecef = hray.direction.north_dir(x_ecef[self.slice_in], y_ecef[self.slice_in],
                                                  z_ecef[self.slice_in], vec_norm_ecef,
                                                  ellps=self.ellps)
        del x_ecef, y_ecef, z_ecef
        self.vec_norm_enu = hray.transform.ecef2enu_vector(vec_norm_ecef, self.trans_ecef2enu)
        self.vec_north_enu = hray.transform.ecef2enu_vector(vec_north_ecef, self.trans_ecef2enu)
        del vec_norm_ecef, vec_north_ecef

        # Merge vertex coordinates and pad geometry buffer
        self.vert_grid = hray.auxiliary.rearrange_pad_buffer(x_enu, y_enu, z_enu)

         # Compute rotation matrix (global ENU -> local ENU)
        rot_mat_glob2loc = hray.transform.rotation_matrix_glob2loc(self.vec_north_enu,
                                                                   self.vec_norm_enu)

        # Compute slope (in global ENU coordinates!)
        slice_in_a1 = (slice(self.slice_in[0].start - 1, self.slice_in[0].stop + 1),
                       slice(self.slice_in[1].start - 1, self.slice_in[1].stop + 1))
        self.vec_tilt_enu = np.ascontiguousarray(hray.topo_param.slope_plane_meth(
                x_enu[slice_in_a1], y_enu[slice_in_a1], z_enu[slice_in_a1],
                rot_mat=rot_mat_glob2loc, output_rot=False)[1:-1, 1:-1])

        del rot_mat_glob2loc
        del x_enu, y_enu, z_enu


    def compute_horizon(self, n_azimuth=None, return_ds=False):
        if n_azimuth is None:
            azim_num = self.n_azimuth
        else:
            azim_num = n_azimuth

        offset_0 = self.slice_in[0].start
        offset_1 = self.slice_in[1].start
        dem_dim_0, dem_dim_1 = self.elevation.shape

        # Compute horizon
        self.hori, self.azimuth = hray.horizon.horizon_gridded(
            self.vert_grid, dem_dim_0, dem_dim_1,
            self.vec_norm_enu, self.vec_north_enu,
            offset_0, offset_1,
            dist_search=self.dist_search,
            azim_num=azim_num)

        if return_ds:

            ds = xr.Dataset(
                coords=dict(
                    azim=(["azim"], self.azimuth, {"units": "radian"}),
                    lat=(["lat"], self.lat[self.slice_in[0]], {"units": "degree"}),
                    lon=(["lon"], self.lon[self.slice_in[1]], {"units": "degree"}),
                ),
                data_vars=dict(
                    horizon=(["azim", "lat", "lon"], np.moveaxis(self.hori, 2, 0),
                               {"units": "radian"})
                )
            )
            return ds


    def compute_slope_aspect_svf(self):
        if self.hori is None:
            self.compute_horizon()


        # Compute Sky View Factor
        svf = hray.topo_param.sky_view_factor(self.azimuth, self.hori, self.vec_tilt)

        # Compute slope angle and aspect
        slope = np.arccos(self.vec_tilt[:, :, 2].clip(max=1.0))
        aspect = np.pi / 2.0 - np.arctan2(self.vec_tilt[:, :, 1],
                                          self.vec_tilt[:, :, 0])
        aspect[aspect < 0.0] += np.pi * 2.0  # [0.0, 2.0 * np.pi]

        ds = xr.Dataset(
                coords=dict(
                    lat=(["lat"], self.lat[self.slice_in[0]], {"units": "degree"}),
                    lon=(["lon"], self.lon[self.slice_in[1]], {"units": "degree"}),
                ),
                data_vars=dict(
                    elevation=(["lat", "lon"], self.elevation[self.slice_in],
                               {"long_name": "ellipsoidal height", "units": "m"}),
                    slope=(["lat", "lon"], slope,
                           {"long_name": "slope angle", "units": "radian"}),
                    aspect=(["lat", "lon"], aspect,
                            {"long_name": "slope aspect (clockwise from North)",
                            "units": "radian"}),
                    svf=(["lat", "lon"], svf,
                         {"long_name": "sky view factor", "units": "-"}),
                )
            )

        return ds


#%%
#========= Example usage ========

if __name__ == "__main__":

    print("WARNING: Example TBI")




