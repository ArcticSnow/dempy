import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    import matplotlib.pyplot as plt
    import xarray as xr
    import numpy as np
    from dempy.data import open_dataset
    from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
    from dempy.fast_shadows import VectorisedShadowIterator


    dem = (open_dataset("dem_alp")
               .isel(band=0).band_data
               .isel(x=slice(3000, 3500), y=slice(3500,4000))
               )


@app.function
def compute_shadows(dem, sun_azimuth, sun_elevation):
    x_center = dem.x.mean().item()
    y_center = dem.y.mean().item()

    dem_centered = dem.copy()
    dem_centered.coords['x'] = dem.x - x_center
    dem_centered.coords['y'] = dem.y - y_center


    target_grid = get_target_grid(dem_centered)
    demrot = get_rotated_dem(dem_centered, target_grid, sun_azimuth)
    SI = VectorisedShadowIterator(demrot, sun_elevation)
    SI.compute()
    shadowsrot = SI.shadows
    shadows_centered = get_rotated_dem(SI.shadows, dem_centered, -1*sun_azimuth)

    shadows = shadows_centered.copy()
    shadows.coords['x'] = shadows.x + x_center
    shadows.coords['y'] = shadows.y + y_center
    return shadows


@app.cell
def _():
    theta = 0.1
    Z_sun = 0.3
    _fig0, _ax0 = plt.subplots()
    _fig1, _ax1 = plt.subplots()

    shadows = compute_shadows(dem, theta, Z_sun)

    (1-shadows).plot(ax=_ax1)
    dem.plot(ax=_ax0)

    mo.hstack([_ax0, _ax1])
    return


if __name__ == "__main__":
    app.run()
