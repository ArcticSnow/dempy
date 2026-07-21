"""
Python class to compute shadow mask from a DEM and timestamps. 
S. Filhol, July 2026


TODO:
1. [x] Try online example for MB
2. [x] Try with EPSG:4326 projected DEM
3. [-] implement python class
4. [ ] apply to all images

"""
#%%
#========= Librairies imports ========
import pdb
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from skyfield.api import load, wgs84
import horayzon as hray
from pathlib import Path
import pandas as pd
from matplotlib.widgets import MultiCursor


#%%
#========= Class definition ========

class DEMShadowMapper:
    def __init__(self, dem_file, domain, dist_search=0.5):
        """
        Initialize class.

        Args:
            dem_file (str or Path): Path to DEM projected in WGS84 (ESPG:4326)
            domain (dict): lat/lon bounding box of area fo interest
            dist_search (float): buffer distance around the area of interest to use for shadow computation

        TODO:
        - [ ] add function to compute for multiple timestamp and store output in dataset
        - [ ] finish reproject_ds() function
        - [ ] use this class as a base to do compute horizons and svf
        - [ ] add check if DEM is in WGS84, otherwise raise error.

        """
        self.file_dem = Path(dem_file)
        self.path_out = self.file_dem.parent
        self.domain = domain
        self.dist_search = dist_search
        self.ellps = "WGS84"

        self.terrain = hray.shadow.Terrain()
        self.trans_ecef2enu = None
        self.slice_in = None
        self.vec_tilt_enu = None

        self._load_dem()


    def _load_dem(self):

        self.domain_outer = hray.domain.curved_grid(self.domain, self.dist_search, self.ellps)
        self.lon, self.lat, self.elevation = hray.load_dem.wgs84_dem(self.file_dem, self.domain_outer)

        # Compute indices of inner domain
        self.slice_in = (slice(np.where(self.lat >= self.domain["lat_max"])[0][-1],
                          np.where(self.lat <= self.domain["lat_min"])[0][0] + 1),
                    slice(np.where(self.lon <= self.domain["lon_min"])[0][-1],
                          np.where(self.lon >= self.domain["lon_max"])[0][0] + 1))
        offset_0 = self.slice_in[0].start
        offset_1 = self.slice_in[1].start
        print("Inner domain size: " + str(self.elevation[self.slice_in].shape))
        elevation_ortho = np.ascontiguousarray(self.elevation[self.slice_in])


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
                                                  ellps=ellps)
        del x_ecef, y_ecef, z_ecef
        vec_norm_enu = hray.transform.ecef2enu_vector(vec_norm_ecef, self.trans_ecef2enu)
        vec_north_enu = hray.transform.ecef2enu_vector(vec_north_ecef, self.trans_ecef2enu)
        del vec_norm_ecef, vec_north_ecef

        # Merge vertex coordinates and pad geometry buffer
        vert_grid = hray.auxiliary.rearrange_pad_buffer(x_enu, y_enu, z_enu)

        # Compute rotation matrix (global ENU -> local ENU)
        rot_mat_glob2loc = hray.transform.rotation_matrix_glob2loc(vec_north_enu,
                                                                   vec_norm_enu)
        del vec_north_enu

        # Compute slope (in global ENU coordinates!)
        slice_in_a1 = (slice(self.slice_in[0].start - 1, self.slice_in[0].stop + 1),
                       slice(self.slice_in[1].start - 1, self.slice_in[1].stop + 1))
        self.vec_tilt_enu = np.ascontiguousarray(hray.topo_param.slope_plane_meth(
                x_enu[slice_in_a1], y_enu[slice_in_a1], z_enu[slice_in_a1],
                rot_mat=rot_mat_glob2loc, output_rot=False)[1:-1, 1:-1])

        # Compute surface enlargement factor
        surf_enl_fac = 1.0 / (vec_norm_enu * self.vec_tilt_enu).sum(axis=2)
        # surf_enl_fac[:] = 1.0
        print("Surface enlargement factor (min/max): %.3f" % surf_enl_fac.min()
              + ", %.3f" % surf_enl_fac.max())

        # Initialise terrain
        mask = np.ones(self.vec_tilt_enu.shape[:2], dtype=np.uint8)
        dim_in_0, dim_in_1 = self.vec_tilt_enu.shape[0], self.vec_tilt_enu.shape[1]
        self.terrain.initialise(vert_grid, dem_dim_0, dem_dim_1,
                                offset_0, offset_1, self.vec_tilt_enu, vec_norm_enu,
                                surf_enl_fac, mask=mask, elevation=elevation_ortho,
                                refrac_cor=True)

    def get_solarposition(self, tstep):
        """
        Function to copmute sun position at a given timestep for the DEM

        Args:
            tstep (DateTime): timestamp in Python Datetime object (timezone aware)

        Return:
            sun position
            sun azimuth
            sun elevation angle 
        """

        # Load Skyfield data
        load.directory = self.path_out
        planets = load("de421.bsp")
        sun = planets["sun"]
        earth = planets["earth"]
        loc_or = earth + wgs84.latlon(self.trans_ecef2enu.lat_or, self.trans_ecef2enu.lon_or)    # -> position lies on the surface of the ellipsoid by default (elevation = 0)

        ts = load.timescale()
        t = ts.from_datetime(tstep)
        astrometric = loc_or.at(t).observe(sun)
        alt, az, d = astrometric.apparent().altaz()
        x = d.m * np.cos(alt.radians) * np.sin(az.radians)
        y = d.m * np.cos(alt.radians) * np.cos(az.radians)
        z = d.m * np.sin(alt.radians)
        sun_position = np.array([x, y, z], dtype=np.float32)

        return sun_position, alt, az


    def compute_sw_dir_cor_map(self, tstep):
        """
        Function to compute direct shortwave correction factor (cf. Horayzon) for a single timestep

        Args:
            tstep (pydatetime): timestamp in python datetime format

        Return:
            xarray dataset of the georeferenced map (WGS84, EPSG:4326). value bounded to 0 (shadow)
        """

        # convert pandas datetime to python datetime object
        #tvec = [x.to_pydatetime() for x in time_vec]

        # create array
        sw_dir_cor_buffer = np.zeros( self.vec_tilt_enu.shape[:2], dtype=np.float32)
        # compute solar position for the given timestamps 
        sun_position, sun_elevation, sun_azimuth = self.get_solarposition(tstep)
        self.terrain.sw_dir_cor(sun_position, sw_dir_cor_buffer)
        out = xr.Dataset(
            {
                "sw_dir_cor": (("lat", "lon", "time"), np.expand_dims(sw_dir_cor_buffer, axis=2)),
                "sun_azimuth": (("time"), [float(sun_azimuth.degrees)]),
                "sun_elevation": (("time"), [float(sun_elevation.degrees)])
            },
            coords={
                "lat": self.lat[self.slice_in[0]],
                "lon": self.lon[self.slice_in[1]],
                "time": [tstep]
            },
        )
        out["sw_dir_cor"].attrs["long_name"] = "Correction factor for direct downward shortwave radiation"
        out["sw_dir_cor"].attrs["more info"] = "sw_dir_cor() according to Müller and Scherer (2005"
        out["sun_azimuth"].attrs["unit"] = "degrees"
        out["sun_elevation"].attrs["unit"] = "degrees"
        out.attrs["reference"] = "Steger, C. R., Steger, B. and Schär, C. (2022): HORAYZON v1.2: an efficient and flexible ray-tracing algorithm to compute horizon and sky view factor, Geosci. Model Dev., 15, 6817–6840, https://doi.org/10.5194/gmd-15-6817-2022"
        out.attrs["horayzon_url"] = "https://github.com/ChristianSteger/HORAYZON"

        out = out.rio.write_crs(4326)
        del sw_dir_cor_buffer
        
        return out
        

    def compute_shadow_map(self, tstep):
        """
        Function to compute direct shadow map (cf. Horayzon) for a single timestep

        Args:
            tstep (pydatetime): timestamp in python datetime format

        Return:
            xarray dataset of the georeferenced map (WGS84, EPSG:4326). value bounded to 0 (shadow)
        """

        # create array
        shadow_buffer = np.zeros( self.vec_tilt_enu.shape[:2], dtype=np.float32)
        # compute solar position for the given timestamps 
        sun_position, sun_elevation, sun_azimuth = self.get_solarposition(tstep)
        self.terrain.shadow(sun_position, shadow_buffer)
        out = xr.Dataset(
            {
                "shadow": (("lat", "lon", "time"), np.expand_dims(shadow_buffer, axis=2)),
                "sun_azimuth": (("time"), [float(sun_azimuth.degrees)]),
                "sun_elevation": (("time"), [float(sun_elevation.degrees)])
            },
            coords={
                "lat": self.lat[self.slice_in[0]],
                "lon": self.lon[self.slice_in[1]],
                "time": [tstep]
            },
        )
        out["shadow"].attrs["long_name"] = "shadow mask (0: illuminated, 1: self-shaded, 2: terrain-shaded, 3: masked)"
        out["shadow"].attrs["more info"] = "shadow() according to Müller and Scherer (2005"
        out["sun_azimuth"].attrs["unit"] = "degrees"
        out["sun_elevation"].attrs["unit"] = "degrees"
        out.attrs["reference"] = "Steger, C. R., Steger, B. and Schär, C. (2022): HORAYZON v1.2: an efficient and flexible ray-tracing algorithm to compute horizon and sky view factor, Geosci. Model Dev., 15, 6817–6840, https://doi.org/10.5194/gmd-15-6817-2022"
        out.attrs["horayzon_url"] = "https://github.com/ChristianSteger/HORAYZON"

        out = out.rio.write_crs(4326)
        del shadow_buffer
        
        return out


    def reproject_ds(ds, epsg_dst, transform, shape, resampling='bilinear'):
        """
        Function to reproject from WGS84 to any projection supported by rasterio, and crop to extent (using transform and shape)
        
        Args:
            ds (dataset or dataarray): xarray dataset in WGS84 with lat/lon dimension
            epsg_dst (int): epsg code of the destination projection 
            transform (object): rasterio transform. 
            shape (array or tuple): shape of the array to reproject
        
        TODO:
        - [ ] add code to check if spatial dimention name are lat/lon or already x/y
        - [ ] integrate better with class usign self. Maybe write independent function that is reussed in class. 
        - [ ] identify why raster tif after this are not readable by QGIS.
        """
        # reproject output to desired EPSG and domain
        dd = ds.rename({'lon':'x','lat':'y'}).rio.reproject(epsg_dst, transform=transform, resampling=resampling, shape=shape)
        return dd

#%%
#========= Example usage ========

if __name__ == "__main__":

    domain = {"lon_min": 6.935, "lon_max": 7.019,
                  "lat_min": 45.91, "lat_max": 45.94}
    dem_file = "/home/filhols/PAPROG/SteepSnow/data/GIS/dem_4326.tif"
    dist_search = 0.25  # search distance for terrain shading [kilometre]
    ellps = "WGS84"  # Earth's surface approximation (sphere, GRS80 or WGS84)


    # Initialize shadow Mapper
    mapper = DEMShadowMapper(
        dem_file=dem_file,
        domain=domain,
        dist_search=0.5, # [km]
    )

    # Loop through all timestamp
    print(f"---> processing timestamp {row.tst.strftime('%Y-%m-%d')}")
    ds_series = mapper.compute_sw_dir_cor_map(row.tst.to_pydatetime())
    
    ds = xr.open_dataset(row.file)
    ds = ds.isel(band=0).rename({'band_data':'im'})
        
    # reproject output to desired EPSG and domain
    dd = ds_series.sw_dir_cor.isel(time=0).rename({'lon':'x','lat':'y'}).rio.reproject(2154, transform=ds.rio.transform(), resampling='bilinear', shape=ds.im.shape)
    dd.attrs['_FillValue']=-9999
    #pdb.set_trace()
    ((dd!=0)*1).rio.to_raster(f'data/zone_{massif}/shadow_mask_{row.tst.strftime("%Y%m%d")}.tif')

    # code to plot example output
    plot=True
    if plot:
        # plot example map
        fig, ax = plt.subplots(2,1,sharex=True, sharey=True)
        dd.plot(ax=ax[0])
        ds.im.plot(cmap=plt.cm.gray, ax=ax[1])
        cursor = MultiCursor(fig.canvas, (ax[0], ax[1]), color='r',lw=0.5, horizOn=True, vertOn=True)
        plt.show()
