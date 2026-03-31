"""
Collection of classes to work with DEM using xarray
S. Filhol, March 2026

- [ ]
"""

import xarray as xr
import rioxarray
import zarr
import numpy as np
from morphometrics import calculate_slope, calculate_aspect, calculate_tpi, calculate_tri, calculate_flow_accumulation, calculate_spi, calculate_twi
from typing import Tuple

class DEMProcessor:
    def __init__(self, dem_path: str):
        self.dem_path = dem_path
        self.dem = None
        self.results = None
        self.morphmetrics = ['slope', 'aspect', 'tpi', 'tri', ]

    def open_dem(self) -> xr.DataArray:
        """Open the DEM as an xarray DataArray."""
        self.dem = xr.open_dataset(self.dem_path, engine='rasterio')

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

    def _get_optimal_chunks(self, shape: Tuple[int, int], target_size: int = 1024) -> Tuple[int, int]:
        """
        Calculate optimal chunk sizes for zarr storage.

        Args:
            shape: Shape of the array (rows, cols).
            target_size: Target chunk size in KB (default: 1024 KB = 1 MB per chunk).

        Returns:
            Tuple of chunk sizes (chunk_rows, chunk_cols).
        """
        rows, cols = shape
        # Estimate the size of one element (float32 = 4 bytes)
        element_size = 4
        # Target elements per chunk
        target_elements = (target_size * 1024) // element_size

        # Calculate chunk sizes
        chunk_rows = min(rows, int(np.sqrt(target_elements * rows / cols)))
        chunk_cols = min(cols, int(np.sqrt(target_elements * cols / rows)))

        # Ensure chunk sizes are at least 1 and divide the array dimensions
        chunk_rows = max(1, chunk_rows)
        chunk_cols = max(1, chunk_cols)
        chunk_rows = min(rows, chunk_rows)
        chunk_cols = min(cols, chunk_cols)

        return (chunk_rows, chunk_cols)

    def save_to_zarr(self, zarr_path: str, target_chunk_size: int = 1024):
        """
        Save the results as a zarr dataset with smart chunking.

        Args:
            zarr_path: Path to save the zarr dataset.
            target_chunk_size: Target chunk size in KB (default: 1024 KB = 1 MB per chunk).
        """
        if self.results is None:
            self.compute_morphometrics()

        # Get the shape of the first variable to determine chunk sizes
        example_var = list(self.results.data_vars.values())[0]
        chunk_rows, chunk_cols = self._get_optimal_chunks(example_var.shape, target_chunk_size)

        # Set encoding with the calculated chunk sizes
        encoding = {var: {'chunks': (chunk_rows, chunk_cols), 'compressor': zarr.Blosc(cname='zstd', clevel=5)}
                    for var in self.results.data_vars}

        # Save to zarr
        self.results.to_zarr(
            zarr_path,
            mode='w',
            encoding=encoding,
            consolidated=True  # Improves metadata access performance
        )
        print(f"Results saved to {zarr_path} with chunks {chunk_rows}x{chunk_cols}")
# Example usage
if __name__ == "__main__":
    dem_path = "path/to/your/dem.tif"
    zarr_path = "path/to/output/morphometrics.zarr"

    processor = DEMProcessor(dem_path)
    processor.open_dem()
    processor.compute_morphometrics()
    processor.save_to_zarr(zarr_path)
