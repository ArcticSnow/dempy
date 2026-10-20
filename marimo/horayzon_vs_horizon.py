import marimo

__generated_with = "0.25.1"
app = marimo.App(width="full")

with app.setup:
    import xarray as xr
    import pandas as pd
    import numpy as np
    from dempy.shadow import DEMShadowMapper
    import matplotlib.pyplot as plt
    from datetime import datetime
    from skyfield.api import utc
    from skyfield.api import load, wgs84
    from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
    from dempy.fast_shadows import HorizonComputer, compute_lon_lat_from_lambert_xy


    dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"
    dem = xr.open_dataset(dem_path).band_data
    # dem_4326 = dem.rio.reproject(4326)
    # dem_4326.to_netcdf("dem_4326.nc")
    dem_4326 = xr.open_dataset("dem_4326.nc")
    # assez capricieux le ShadowMapper pour le choix du domaine ...
    domain = {"lon_min": 6.935, "lon_max": 7.02, "lat_min": 45.915, "lat_max": 45.935}
    dist_search = 0.25  # search distance for terrain shading [kilometre]
    ellps = "WGS84"  # Earth's surface approximation (sphere, GRS80 or WGS84)
    # Initialize shadow Mapper
    mapper = DEMShadowMapper(
        dem_file="dem_4326.nc",
        domain=domain,
        dist_search=0.5, # [km]
    )
    date = datetime(2019,3,22,10,32,45, tzinfo=utc)
    ds_series = mapper.compute_sw_dir_cor_map(date)
    sw_dir_2154 = ds_series.rename(dict(lon="x", lat="y")).sw_dir_cor.transpose('time', 'y', 'x').rio.reproject(2154).isel(time=0)

    # calcul de lat et lon sur le DEM + lat/lon du point central
    _lon, _lat = compute_lon_lat_from_lambert_xy(dem.isel(band=0))
    dem = dem.assign_coords(lon=_lon, lat=_lat)
    j0 = dem.sizes["y"]//2
    i0 = dem.sizes["x"]//2
    lat_center = dem.lat.isel(y=j0,x=i0).item()
    lon_center = dem.lon.isel(y=j0,x=i0).item()


@app.cell
def _():
    sw_reproj_match = ds_series.rename(dict(lon="x", lat="y")).sw_dir_cor.transpose('time', 'y', 'x').rio.reproject_match(dem).isel(time=0)
    return (sw_reproj_match,)


@app.cell
def _(sw_reproj_match):
    (sw_reproj_match > 0).plot()
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    sw_dir_2154.x.values, dem.x.values
    return


@app.function
def get_sun_azimuth_elevation(lat0, lon0, date_):
    ts = load.timescale()
    t = ts.utc(date_.year, date_.month, date_.day, date_.hour, date_.minute, date_.second)
    planets = load("de421.bsp")
    sun = planets["sun"]
    earth = planets["earth"]
    # Position de l'observateur
    observer = earth + wgs84.latlon(
        latitude_degrees=lat0,
        longitude_degrees=lon0,
        elevation_m=300,
    )
    # Position du Soleil vue depuis l'observateur
    astrometric = observer.at(t).observe(sun)
    apparent = astrometric.apparent()
    alt, az, distance = apparent.altaz()
    return az.degrees, alt.degrees


@app.cell
def _(azi):
    np.cos(np.deg2rad(180-azi))
    return


@app.cell
def _(azi):
    sun_azimuth = np.deg2rad(azi)
    x_center = dem.x.mean().item()
    y_center = dem.y.mean().item()
    dem_centered = dem.copy().isel(band=0)
    dem_centered.coords['x'] = dem.x - x_center
    dem_centered.coords['y'] = dem.y - y_center
    target_grid = get_target_grid(dem_centered)
    target_grid    
    dem_rotation_angle = sun_azimuth - np.pi
    demrot = get_rotated_dem(dem_centered, target_grid, dem_rotation_angle)
    horizon = HorizonComputer(demrot)
    for _step_size in 20*[1] + 20*[4]:
        horizon.step(_step_size)
    horizon_for_date = get_rotated_dem(horizon.horizons[-1], dem_centered, -1*dem_rotation_angle)
    horizon_for_date.coords['x'] = horizon_for_date.x + x_center
    horizon_for_date.coords['y'] = horizon_for_date.y + y_center
    horizon_for_date = np.atan(horizon_for_date*np.cos(np.deg2rad(180-azi)))
    horizon_for_date.plot()
    return (horizon_for_date,)


@app.cell
def _():

    azi, elev = get_sun_azimuth_elevation(lat_center, lon_center, date)

    azi, elev, mapper.get_solarposition(date)
    return azi, elev


@app.cell
def _():
    42+23/60+19.2/3600, 157+16/60
    return


@app.cell
def _():
    dem_4326.band_data.plot()
    return


@app.cell
def _():
    ds_series.sw_dir_cor.plot()
    return


@app.cell
def _(sw_reproj_match):
    shadow1 = (sw_reproj_match > 0)

    shadow1.plot()
    return (shadow1,)


@app.cell
def _(elev, horizon_for_date):
    shadow2 = (horizon_for_date <= np.deg2rad(elev) +0.3)
    shadow2.plot()
    return (shadow2,)


@app.cell
def _(shadow1, shadow2):
    diff_shadows = shadow2.astype('int') - shadow1.astype('int')
    return (diff_shadows,)


@app.cell(expand_output=True)
def _(diff_shadows, mo):
    _fig, _ax = plt.subplots(figsize=(15,10))

    diff_shadows.plot(ax=_ax)

    mo.mpl.interactive(_ax)
    return


@app.cell
def _(diff_shadows):
    dom = dict(y = slice(6.5448e6, 6.544e6), x = slice(1.008e6, 1.010e6))
    (diff_shadows.sel(**dom) == 1).sum().item(), (diff_shadows.sel(**dom) == -1).sum().item(), (diff_shadows.sel(**dom) == 0).sum().item()
    return


@app.cell
def _(diff_shadows):
    diff_shadows
    return


if __name__ == "__main__":
    app.run()
