import pooch
import xarray as xr



def open_dataset(name):
    if name == "dem_alp":
        path = pooch.retrieve(
            url="https://github.com/roughnemezis/dempy-data/releases/download/v0.0.1/DEM_alp_L93_bilinear.tif",
            known_hash="2559f325b65d079dab33ad5e6238352cbaf467c596178041298c2a5e7ead5e4a"
            )
        dem = xr.open_dataset(path)
        return dem

