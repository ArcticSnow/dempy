import xarray as xr
import numpy as np
from dempy.shadow import DEMShadowMapper
import matplotlib.pyplot as plt


dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"

dem = xr.open_dataset(dem_path).isel(band=0).band_data





domain = {"lon_min": 6.935, "lon_max": 7.019, "lat_min": 45.91, "lat_max": 45.94}
dist_search = 0.25  # search distance for terrain shading [kilometre]
ellps = "WGS84"  # Earth's surface approximation (sphere, GRS80 or WGS84)


# Initialize shadow Mapper
mapper = DEMShadowMapper(
    dem_file=dem_path,
    domain=domain,
    dist_search=0.5, # [km]
)

dem.rio.reproject(4326).plot()





dds = []
for i, row in df_meta.iterrows():
    print(f"---> processing timestamp {row.tst.strftime('%Y-%m-%d')}")
    ds_series = mapper.compute_sw_dir_cor_map(row.tst.to_pydatetime())

    ds = xr.open_dataset(row.file)
    ds = ds.isel(band=0).rename({'band_data':'im'})
    
    # reproject output to desired EPSG and domain
    dd = ds_series.sw_dir_cor.isel(time=0).rename({'lon':'x','lat':'y'}).rio.reproject(2154, transform=ds.rio.transform(), resampling='bilinear', shape=ds.im.shape)
    dd.attrs['_FillValue']=-9999
    #pdb.set_trace()
    
    # export to geotiff (not sure why not working with Qgis)
    #((dd!=0)*1).rio.to_raster(f'data/zone_{massif}/shadow_mask_{row.tst.strftime("%Y%m%d")}.tif')
    # TODO: store data as dataset instead of tiff
    dds.append(((dd!=0)*1))

fs = xr.concat(dds, dim='time')
fs.to_netcdf(f'data/zone_{massif}/ds_shadows.nc')


