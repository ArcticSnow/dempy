import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")

with app.setup:
    import matplotlib.pyplot as plt
    import xarray as xr
    import numpy as np
    from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
    from dempy.fast_shadows import HorizonComputer, compute_lon_lat_from_lambert_xy
    from skyfield.api import load, wgs84
    from datetime import datetime
    import pandas as pd

    dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"

    dem = xr.open_dataset(dem_path).isel(band=0).band_data
    _lon, _lat = compute_lon_lat_from_lambert_xy(dem)
    dem = dem.assign_coords(lon=_lon, lat=_lat)

    date_ = datetime(2025,2,3,5,12,13)
    tst_str = [20170610103302, 20180602103656, 20190303102906, 20190322103245, 20190513103258,
               20190710103709, 20200708103707, 20210119103654, 20210126103305, 20210214103639,
               20210523103332, 20211203104105, 20211216104112, 20220111104107, 20220124104049,
               20220713103309]
    tst = pd.to_datetime(tst_str, format='%Y%m%d%H%M%S')
    # dem = xr.open_dataset(dem_path).isel(band=0).band_data.isel(x=slice(None, None, 10), y=slice(None, None, 10))


    j0 = dem.sizes["y"]//2
    i0 = dem.sizes["x"]//2
    # R = 6371 # rayon de la Terre en km

    lat0 = dem.lat.isel(y=j0,x=i0).item()
    lon0 = dem.lon.isel(y=j0,x=i0).item()

    # dem = dem[::50, ::50]


@app.cell
def _(mo):
    sun_azimuth, sun_elevation = get_sun_azimuth_elevation(lat0, lon0, tst[3])
    sun_azimuth = np.deg2rad(sun_azimuth)
    sun_elevation = np.deg2rad(sun_elevation)
    mo.md(f"""
          Picked date is **{tst[3]}**\n
          Sun Azimuth in radian **{sun_azimuth}**\n
          Sun Elevation in radian **{sun_elevation}**
          """)
    return sun_azimuth, sun_elevation


@app.cell
def _(sun_azimuth):
    x_center = dem.x.mean().item()
    y_center = dem.y.mean().item()
    dem_centered = dem.copy()
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
    horizon_for_date.plot()
    return (horizon_for_date,)


@app.cell
def _():
    dem
    return


@app.cell
def _(horizon_for_date):
    horizon_rad = np.atan(horizon_for_date)
    horizon_rad.plot()
    return (horizon_rad,)


@app.cell
def _(horizon_rad, sun_elevation):
    shadows = (horizon_rad < sun_elevation)
    return (shadows,)


@app.cell
def _(shadows):
    shadows
    return


@app.cell
def _(shadows):
    shadows.plot()
    return


@app.cell
def _():
    import rioxarray

    def save_as_tiff(mask, path):
        mask.astype("uint8").rio.to_raster(
            path,
            dtype="uint8",
            compress="DEFLATE",
        )

    return (save_as_tiff,)


@app.cell
def _(save_as_tiff, shadows):
    from pathlib import Path

    stor_path = Path("/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/shadows_simon")
    save_as_tiff(shadows,stor_path/"dem_corp_verte_2019_03_22T10_32_45.tif")
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


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Calcul des ombrages sans prise en compte de la rotation de la Terre
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - [ ] calculer lat et lon sur le dem
    - [ ]
    """)
    return


@app.cell
def _():
    import marimo as mo

    return (mo,)


if __name__ == "__main__":
    app.run()
