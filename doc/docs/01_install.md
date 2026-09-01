# Installing

This project uses [pixi](https://pixi.prefix.dev/latest/) as a package manager. You don't have to use `pixi`(which itself builds a standard `pyproject.toml`) but we will only document installation using `pixi`.

- follow [installation instructions in pixi documentations](https://pixi.prefix.dev/latest/installation/) it should be a matter of seconds
- clone this project from github
```bash
git clone git@github.com:ArcticSnow/dempy.git
```
- `cd dempy`
- `pixi install` will create a `.pixi` folder containing the virtual env for the project, install everything in this virtualenv
- `pixi shell` will put you in the project environment (`python` executable is that of the venv)
- `pixi run cmd`  will run `cmd` in the project environment (run for example `pixi run marimo edit marimo/algo_vectorized_marimo.py`to tun your first notebook


# Working

## Getting data to work with

DEM data is shipped with the project, to run the examples or develop new features. It's a 30m DEM on the french Alps by IGN, the data module allows you to access it.

```python
from dempy.data import open_dataset
dem = open_dataset("dem_alp").isel(band=0).band_data
```

Et voilà. Some data to play with.

We use the `pooch` library to get the data. Data is stored online, `pooch` will download it and cache it on your computer so you only have to download it once.

