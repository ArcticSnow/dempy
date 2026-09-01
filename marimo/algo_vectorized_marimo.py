import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Démonstration de l'algorithme vectorisé sur un DEM de 2000 par 2000

    Le calcul prend dans les 5s ☺️
    """)
    return


@app.cell
def _():
    import marimo as mo
    from dempy.data import open_dataset
    import numpy as np
    import matplotlib.pyplot as plt
    import xarray as xr
    from dempy.fast_shadows import VectorisedShadowIterator

    return VectorisedShadowIterator, mo, np, open_dataset, plt


@app.cell
def _(open_dataset):
    dem = open_dataset("dem_alp").isel(band=0).band_data
    return (dem,)


@app.cell
def _(dem):
    dem[::100, ::100].plot()
    return


@app.cell
def _(VectorisedShadowIterator, dem):
    SI = VectorisedShadowIterator(dem.isel(x=slice(2000, 4000), y=slice(2000,4000)), 0.1)
    SI.compute()
    return (SI,)


@app.cell
def _(SI):
    SI.shadows.plot()
    return


@app.cell
def _(SI, np, plt):
    remaining_steps = [step.sum().item() for step in SI.steps]

    _fig, _ax = plt.subplots(figsize=(12,8))
    _ax.plot(SI.dem.y.values[::-1], np.array(remaining_steps[:-1]))
    _in_steps = SI.steps[0].cumsum()
    # _ax.plot(SI.steps[0].cumsum().y.values, _in_.values)
    _in_steps.plot(ax=_ax)
    _ax.set_ylabel("nombre de points")
    _ax.set_title("nombre d'ombrages à calculer")
    _ax
    return


if __name__ == "__main__":
    app.run()
