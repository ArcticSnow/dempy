"""
Collection of classes to work with DEM using xarray
S. Filhol, March 2026

- [ ]
"""

import xarray as xr
import rioxarray
import zarr
import numpy as np
from morphometrics import (
    calculate_slope, calculate_aspect, calculate_tpi, calculate_tri,
    calculate_flow_accumulation, calculate_spi, calculate_twi
)

class DEMProcessor:
    def __init__(self, dem_path: str):
        self.dem_path = dem_path
        self.dem = None
        self.results = None

    def open_dem(self) -> xr.DataArray:
        """Open the DEM as an xarray DataArray."""
        self.dem = rioxarray.open_rasterio(self.dem_path, masked=True).squeeze()
        return self.dem

    def compute_morphometrics(self, tpi_radius: int = 3, tri_radius: int = 1) -> xr.Dataset:
        """Compute all morphometric indices."""
        if self.dem is None:
            self.open_dem()

        dem_values = self.dem.data
        slope = calculate_slope(dem_values)
        aspect = calculate_aspect(dem_values)
        tpi = calculate_tpi(dem_values, radius=tpi_radius)
        tri = calculate_tri(dem_values, radius=tri_radius)
        flow_acc = calculate_flow_accumulation(dem_values)
        spi = calculate_spi(flow_acc, slope)
        twi = calculate_twi(flow_acc, slope)

        self.results = xr.Dataset({
            'slope': (('y', 'x'), slope),
            'aspect': (('y', 'x'), aspect),
            'tpi': (('y', 'x'), tpi),
            'tri': (('y', 'x'), tri),
            'flow_accumulation': (('y', 'x'), flow_acc),
            'spi': (('y', 'x'), spi),
            'twi': (('y', 'x'), twi),
        })

        # Copy coordinates from the original DEM
        self.results = self.results.assign_coords({
            'x': self.dem.x,
            'y': self.dem.y,
        })

        return self.results

    def save_to_zarr(self, zarr_path: str):
        """Save the results as a zarr dataset."""
        if self.results is None:
            self.compute_morphometrics()

        self.results.to_zarr(zarr_path, mode='w')
        print(f"Results saved to {zarr_path}")

# Example usage
if __name__ == "__main__":
    dem_path = "path/to/your/dem.tif"
    zarr_path = "path/to/output/morphometrics.zarr"

    processor = DEMProcessor(dem_path)
    processor.open_dem()
    processor.compute_morphometrics()
    processor.save_to_zarr(zarr_path)
