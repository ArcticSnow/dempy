import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return


@app.cell
def _():
    from dempy.fast_shadows import compute_shadows
    from dempy.data import open_dataset
    import numpy as np

    dem = (open_dataset("dem_alp")
               .isel(band=0).band_data
                #.isel(x=slice(3000, 3500), y=slice(3500,4000))
               .isel(x=slice(None, None,10), y=slice(None, None, 10))
               )
    return compute_shadows, dem, np


@app.cell
def _(compute_shadows, dem, np):
    shadows = compute_shadows(dem, np.pi/2, np.pi/30)
    return (shadows,)


@app.cell
def _(shadows):
    (1-shadows).plot()
    return


if __name__ == "__main__":
    app.run()
