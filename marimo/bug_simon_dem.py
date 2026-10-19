import marimo

__generated_with = "0.25.1"
app = marimo.App(width="columns")


@app.cell(column=0)
def _(dem_test):
    demo = dem_test[0]
    return (demo,)


@app.cell
def _(demo):
    demo
    return


@app.cell
def _(demo):
    demo.plot()
    return


@app.cell
def _(demo, np):
    (demo - demo.shift(y=-10) - 10*np.tan(0.15)).plot()
    return


@app.cell
def _(demo, xr):
    shadowso = xr.zeros_like(demo).astype(bool)
    return (shadowso,)


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
    (dem.x[1] - dem.x[0]).item()
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


@app.cell
def _(dem, h_tangent):
    ix0 = 1500
    ix1 = 3000
    iy0 = 0
    iy1 = 1500

    dem_test = [
        (dem - h_tangent).isel(x=slice(ix0,ix1), y = slice(iy0, iy1)),
        (dem - h_tangent).isel(x=slice(ix0,ix1,5), y = slice(iy0, iy1,5)),
        (dem - h_tangent).isel(x=slice(ix0,ix1,10), y = slice(iy0, iy1,10)),
    ]
    return (dem_test,)


@app.cell
def _(compute_shadows, dem_test, np, sun_positions, xr):
    _azi = sun_positions[0]["azi"]
    _elev = sun_positions[0]["elevation"]
    _azi = np.deg2rad(180)
    _elev = 0.15
    shadows = [xr.Dataset(dict(mask=compute_shadows(dem, _azi, _elev))) for dem in dem_test]
    return (shadows,)


@app.cell
def _(shadows):
    (1-shadows[0].mask).plot()
    return


@app.cell
def _(plt, shadows):
    fig_old, _axs = plt.subplots(1,3, figsize=(20,5))

    for _sh, _ax in zip(shadows, _axs):
        (1-_sh.mask).plot(ax=_ax)

    fig_old
    return (fig_old,)


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


    return (get_sun_azimuth_elevation,)


@app.cell
def _():
    from pathlib import Path

    return


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
def _(tst):
    tst
    return


@app.cell
def _(get_sun_azimuth_elevation, lat0, lon0, np, tst):
    sun_positions = []

    for _time in tst:
        _azi, _elev = get_sun_azimuth_elevation(lat0, lon0, _time)
        _azi = np.deg2rad(_azi)
        _elev = np.deg2rad(_elev)
        sun_positions.append(dict(azi=_azi, elevation=_elev))
    return (sun_positions,)


@app.cell
def _(sun_positions):
    sun_positions
    return


@app.cell(column=1)
def _(shadowso):
    shadowso
    return


@app.cell
def _(demo, np, shadowso):
    shadowsl = shadowso

    for shift in range(500):
        shadowsl = shadowsl | ((demo - demo.shift(y=-shift) - shift*np.tan(0.15)) < 0)
    return (shadowsl,)


@app.cell
def _(demo):
    demo.plot()
    return


@app.cell
def _(shadowsl):
    (1-shadowsl).plot()
    return


@app.cell
def _(dem, np, xr):
    horizons = xr.zeros_like(dem)

    for _shift in range(300):
        # diffdemo = np.atan((dem.shift(y=-_shift) - dem)/(1.5*_shift))
        diffdemo = (dem.shift(y=-_shift) - dem)/(1.5*_shift)
        horizons = diffdemo.where(diffdemo >=horizons, horizons)

    horizons = np.atan(horizons)
    horizons.plot()
    return (horizons,)


@app.cell
def _(horizons):
    (horizons < 0.05).plot()
    return


@app.cell
def _(horizons):
    horizons.plot.hist(bins=50)
    return


@app.cell
def _(shadows):
    (1-shadows[0]).mask.plot()
    return


@app.cell
def _(fig_old):
    fig_old
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
