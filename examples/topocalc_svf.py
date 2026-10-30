"""
# Demo for topocalc usage in our repository

Navigate using go to definition feature from your editor to dig into  the original code.

## Environment

This code will work inside the `topocalc` environment
- `pixi shell -e topocalc` then `python`
- or directly `pixi run -e topocalc python`



## Horizon calculation

- the dem is first skewed (skew_transpose function)
- horizon calc by a C function on skewed DEM (no mean to know what's inside)


"""



import pooch
import xarray as xr
import numpy as np
from spatialnc import ipw
import matplotlib.pyplot as plt
from topocalc.viewf import viewf
from topocalc.horizon import horizon
from topocalc.gradient import gradient_d8

path = pooch.retrieve(
    url="https://github.com/roughnemezis/dempy-data/releases/download/v0.0.1/DEM_alp_L93_bilinear.tif",
    known_hash="2559f325b65d079dab33ad5e6238352cbaf467c596178041298c2a5e7ead5e4a"
    )



dem = xr.open_dataset(path)

demtest = dem.isel(band=0).isel(x=slice(None, None, 10), y=slice(None, None,10))


# computing svf with topocalc
svf_topocalc = viewf(demtest.band_data.values.astype('float64'), spacing=30) # douloureusement long...

demtest['svf'] = xr.DataArray(svf_topocalc[0], dims= ('y','x'))

demtest.svf.plot(); plt.show()


gradient_topocalc = gradient_d8(demtest.band_data.values, 30, 30)

demtest['gradient_d8'] = xr.DataArray(gradient_topocalc[0], dims= ('y','x'))

demtest.gradient_d8.plot(); plt.show()


hcos = horizon(70, demtest.band_data.values, 30)

