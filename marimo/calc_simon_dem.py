import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"
    return (dem_path,)


@app.cell
def _(dem_path):
    import matplotlib.pyplot as plt
    import xarray as xr
    import numpy as np
    from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
    from dempy.fast_shadows import compute_shadows, compute_lon_lat_from_lambert_xy


    dem = xr.open_dataset(dem_path).isel(band=0).band_data
    # dem = xr.open_dataset(dem_path).isel(band=0).band_data.isel(x=slice(None, None, 10), y=slice(None, None, 10))

    lon, lat = compute_lon_lat_from_lambert_xy(dem)
    # Attach as coordinates to the DEM for easy indexing
    dem = dem.assign_coords(lon=lon, lat=lat)
    return compute_shadows, dem, lat, lon, np, plt, xr


@app.cell
def _(dem):
    dem
    return


@app.cell
def _(dem, lat, lon, np):
    j0 = dem.sizes["y"]//2
    i0 = dem.sizes["x"]//2
    R = 6371 # rayon de la Terre en km

    lat0 = dem.lat.isel(y=j0,x=i0).item()
    lon0 = dem.lon.isel(y=j0,x=i0).item()

    h_tangent = 1000*R*(1-np.cos(np.deg2rad(lat0))*np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(lon-lon0))-np.sin(np.deg2rad(lat0))*np.sin(np.deg2rad(lat)))
    return h_tangent, lat0, lon0


@app.cell
def _(h_tangent):
    h_tangent.plot() # cétacé négligeable
    return


@app.cell
def _(dem, h_tangent):
    (dem - h_tangent).plot()
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
    from skyfield.api import load, wgs84

    def get_sun_azimuth_elevation(lat0, lon0, date_):
        ts = load.timescale()
        t = ts.utc(2026, 9, 28, 12, 0, 0)
        t = ts.utc(date_.year, date_.month, date_.day, date_.hour, date_.minute, date_.second)
        planets = load("de421.bsp")
        sun = planets["sun"]
        earth = planets["earth"]
        # Position de l'observateur
        observer = earth + wgs84.latlon(        latitude_degrees=lat0,
            longitude_degrees=lon0,
            elevation_m=300,
        )
        # Position du Soleil vue depuis l'observateur
        astrometric = observer.at(t).observe(sun)
        apparent = astrometric.apparent()
        alt, az, distance = apparent.altaz()
        return az.degrees, alt.degrees


    return (get_sun_azimuth_elevation,)


@app.cell
def _():
    from datetime import datetime
    import pandas as pd

    date_ = datetime(2025,2,3,5,12,13)
    tst_str = [20170610103302,
                    20180602103656,
                    20190303102906,   # image has cloud like feature leading to non homogeneous values
                    20190322103245,
                    20190513103258,
                    20190710103709,
                    20200708103707,
                    20210119103654,
                    20210126103305,
                    20210214103639,
                    20210523103332,
                    20211203104105,
                    20211216104112,
                    20220111104107,
                    20220124104049,
                    20220713103309]

    tst = pd.to_datetime(tst_str, format='%Y%m%d%H%M%S')
    return (tst,)


@app.cell
def _(shadow0):
    shadow0.attrs['elevation'] = 10
    shadow0.attrs
    return


@app.cell
def _(save_as_tiff, stor_path, tst, xr):
    for _time in tst:
        _path_shadow = stor_path / f"{_time.strftime('%Y%m%d%H%M%S')}.nc"
        if _path_shadow.exists():
            _shadow = xr.open_dataset(_path_shadow)
            _shadow_asbool = xr.Dataset(dict(mask=_shadow.__xarray_dataarray_variable__.astype(bool)))
            save_as_tiff(_shadow_asbool,stor_path / f"{_time.strftime('%Y%m%d%H%M%S')}.tif")
    return


@app.cell
def _(
    compute_shadows,
    dem,
    get_sun_azimuth_elevation,
    h_tangent,
    lat0,
    lon0,
    np,
    save_as_tiff,
    stor_path,
    tst,
    xr,
):
    for _time in tst:
        path_shadow = stor_path / f"{_time.strftime('%Y%m%d%H%M%S')}.tif"
        if path_shadow.exists():
            print(f"{path_shadow} already exist, ignoring")
            continue
        print(f"Computing mask for {path_shadow}")
        azi, elevation = get_sun_azimuth_elevation(lat0, lon0, _time)
        shadow = xr.Dataset(dict(mask=compute_shadows(dem - h_tangent, np.deg2rad(azi), np.deg2rad(elevation))))
        save_as_tiff(shadow, path_shadow)
    return


@app.cell
def _():
    # shadows_plots = {}
    # 
    # for _time in shadows_dict.keys():
    #     _fig, _ax = plt.subplots()
    #     shadows_plots[_time] = _ax
    #     (1-shadows_dict[_time]).plot(ax=_ax)
    #     
    #  
    # mo.tabs(shadows_plots)
    return


@app.cell
def _(stor_path):
    stor_path.exists()
    return


@app.cell
def _(Path):
    stor_path = Path("/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/shadows_simon")
    return (stor_path,)


@app.cell
def _():
    from pathlib import Path

    return (Path,)


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    - [ ] optimiser le netcdf? (300 Mo pour un tableau de booleens...)
    """)
    return


@app.cell
def _(stor_path, tst, xr):
    shadows_dict = {}


    _time = tst[7]
    path_shad = stor_path / f"{_time.strftime('%Y%m%d%H%M%S')}.tif"
    shadow0 = xr.open_dataset(path_shad)
    shadow0
    return (shadow0,)


@app.cell
def _(plt, shadow0):
    _fig, _ax = plt.subplots(figsize=(10,5))

    (1-shadow0.band_data).plot(ax=_ax)
    return


@app.cell
def _(shadow0):
    shadow0.sum()
    return


@app.cell
def _(shadow0, xr):
    # shadow0_asbool.astype("uint8").to_netcdf(
    #     stor_path / "mask.nc",
    #     engine="netcdf4",
    #     encoding={"mask": {"dtype": "uint8", "zlib": True, "complevel": 9, "shuffle": True}},
    # )
    shadow0_asbool = shadow0.__xarray_dataarray_variable__.astype(bool)
    shadow0_asbool = xr.Dataset(dict(mask=shadow0_asbool))
    # wow en tif on gagne un facteur 1000!!
    return (shadow0_asbool,)


@app.cell
def _():
    import rioxarray

    def save_as_tiff(mask, path):
        mask.astype("uint8").rio.to_raster(
            path,
            dtype="uint8",
            compress="DEFLATE",
        )

    return rioxarray, save_as_tiff


@app.cell
def _(shadow0_asbool):
    shadow0_asbool
    return


@app.cell
def _(rioxarray, stor_path):
    mask = rioxarray.open_rasterio(stor_path/"mask.tif")
    return (mask,)


@app.cell
def _(mask):
    (1-mask).plot()
    return


if __name__ == "__main__":
    app.run()
