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

