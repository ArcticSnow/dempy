import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")

with app.setup:
    import marimo as mo
    import matplotlib.pyplot as plt
    import xarray as xr
    import numpy as np
    from dempy.data import open_dataset
    from dempy.rotate import Transformer, get_rotated_dem, get_target_grid
    from dempy.fast_shadows import VectorisedShadowIterator


    dem = open_dataset("dem_alp").isel(band=0).band_data[::10,::10]
    # -------------------------------------------------
    # Compute latitude and longitude arrays for the DEM (Lambert‑93 → WGS‑84)
    # -------------------------------------------------
    from pyproj import Transformer as ProjTransformer

    # Get 1‑D coordinate vectors
    x_vals = dem.x.values          # shape (x)
    y_vals = dem.y.values          # shape (y)

    # 2️⃣ Create 2‑D grid matching the DEM layout
    X, Y = np.meshgrid(x_vals, y_vals)   # both shape (y, x)

    # 3️⃣ Set up transformer: Lambert‑93 (EPSG:2154) → WGS‑84 (EPSG:4326)
    lambert_to_wgs84 = ProjTransformer.from_crs(
        crs_from="EPSG:2154",
        crs_to="EPSG:4326",
        always_xy=True,
    )

    # 4️⃣ Transform the whole grid at once
    lon_vals, lat_vals = lambert_to_wgs84.transform(X, Y)   # both (y, x)

    # 5️⃣ Wrap as DataArrays, preserving DEM dimensions & coordinates
    lon = xr.DataArray(
        lon_vals,
        dims=dem.dims,
        coords=dem.coords,
        name="longitude",
    )
    lat = xr.DataArray(
        lat_vals,
        dims=dem.dims,
        coords=dem.coords,
        name="latitude",
    )

    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    # 6️⃣ Attach as coordinates to the DEM for easy indexing
    dem = dem.assign_coords(lon=lon, lat=lat)


@app.cell
def _():
    dem.plot()
    return


@app.cell
def _():
    print(dem)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    $$ \boxed{ |d|= R\left[ 1-\cos LAT_0\cos LAT\cos(LON-LON_0) -\sin LAT_0\sin LAT \right] } $$
    """)
    return


@app.cell
def _():
    j0, i0 = 500, 350

    lat0 = lat.isel(y=j0,x=i0).item()
    lon0 = lon.isel(y=j0,x=i0).item()

    R = 13000

    # distance de la sphère (terrestre) plan tangent au point choisi
    h_tangent = 1000*R*(1-np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)-np.sin(lat0)*np.sin(lat))

    h_tangent.plot()
    return (h_tangent,)


@app.cell
def _(h_tangent):
    dem_curv = (dem - h_tangent)
    dem_curv.plot()
    return (dem_curv,)


@app.cell
def _():
    dem.plot()
    return


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
def _(dem_curv):
    shadows = compute_shadows(dem_curv, 1.5, 0.04)
    return (shadows,)


@app.cell
def _(shadows):
    (1-shadows).plot()
    return


@app.cell
def _():
    # Compute hillshade from the DEM
    dx = dem.x.diff('x').mean().item()  # pixel spacing in x (meters)
    dy = dem.y.diff('y').mean().item()  # pixel spacing in y (meters)

    z = dem.values

    # Gradients (note y axis ordering)
    dz_dy, dz_dx = np.gradient(z, dy, dx)

    # Surface normal
    normal = np.dstack((-dz_dx, -dz_dy, np.ones_like(z)))

    # Sun direction
    sun_azimuth_deg = 315.0
    sun_elevation_deg = 45.0
    az = np.deg2rad(sun_azimuth_deg)
    el = np.deg2rad(sun_elevation_deg)

    sun_dir = np.array([
        np.cos(el) * np.sin(az),
        np.cos(el) * np.cos(az),
        np.sin(el),
    ])

    # Dot product for hillshade
    dot = np.sum(normal * sun_dir, axis=2)
    norm_mag = np.sqrt(np.sum(normal**2, axis=2))
    hillshade = np.clip(dot / (norm_mag + 1e-12), 0, 1)

    hillshade_da = xr.DataArray(
        hillshade,    dims=dem.dims,
        coords=dem.coords,
        name="hillshade",
    )

    hillshade_da.plot(cmap="gray", vmin=0, vmax=1)
    return


if __name__ == "__main__":
    app.run()
