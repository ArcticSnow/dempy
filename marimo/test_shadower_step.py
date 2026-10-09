import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    import matplotlib.pyplot as plt
    import xarray as xr
    import numpy as np
    from dempy.fast_shadows import VectorisedShadowIterator

    dem_path = "/Users/merzisenh/Library/CloudStorage/ProtonDrive-h.merzisen@proton.me-folder/DATA/dem_crop_verte_shadow.tif"

    dem = xr.open_dataset(dem_path).isel(band=0).band_data

    demtest = dem.isel(y=slice(None, None, 500), x=slice(None,None,500))
    shadow = VectorisedShadowIterator(demtest, 0.05)


@app.cell
def _():
    _fig, (_ax1, _ax2, _ax3) = plt.subplots(1,3, figsize=(16, 3.5))
    demtest.plot(ax=_ax1, cbar_kwargs={"orientation": "horizontal"})
    shadow.diff_dem.plot(ax=_ax2, cbar_kwargs={"orientation": "horizontal"})
    shadow.peaks_to_explore.plot(ax=_ax3, cbar_kwargs={"orientation": "horizontal"})
    for _ax in (_ax2, _ax3):
        _ax.set_ylabel("")
        _ax.set_yticklabels([])
    _fig
    return


@app.cell
def _():

    _fig, _ax1 = plt.subplots(figsize=(12, 5))
    demtest.isel(x=0).to_series().reset_index().plot(kind="scatter", x="y", y="band_data", ax=_ax1)
    demtest.isel(x=0).to_series().plot(ax=_ax1)
    shadow.top_solar_ray.isel(x=0).plot()
    (2000*shadow.peaks_to_explore).isel(x=0).to_series().reset_index().plot(kind="scatter", x="y", y="band_data", color="red", ax=_ax1)
    return


@app.cell
def _():
    shadow.compute(-1)
    shadow.top_solar_ray.plot()
    return


@app.cell
def _():
    shadow.new_shadows.plot()
    return


@app.cell
def _():
    _fig, (_ax1, _ax2) = plt.subplots(1,2, figsize=(10, 4))
    # demtest.plot(ax=_ax1)
    # shadow.diff_dem.plot(ax=_ax2)

    demtest.isel(x=3).plot(ax=_ax1)
    last_z = demtest.isel(x=3, y=-1).item()
    (-np.tan(0.15)*(demtest.y - demtest.y[-1])+last_z).plot(ax=_ax1)
    shadow.diff_dem.isel(x=3).plot(ax=_ax2)
    (0*demtest.y - np.tan(0.15)*shadow.resolution).plot(ax=_ax2)
    _fig
    # donc le critère c'est: diff_dem < -1*np.tan(0.15)*shadow.resolution
    return


@app.cell
def _():
    shadow.resolution
    shadow.diff_dem
    return


if __name__ == "__main__":
    app.run()
