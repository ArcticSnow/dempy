import numpy as np

def winstral_fast(dem, cellsize, dmax, in_wind):
    '''
    Function to derive Winstral Sx surface from a dem. 
    :param dem: numpy array of a DEM
    :param cellsize: resolution of the DEM
    :param dmax: maximum distance for winstral parameter estimate
    :param in_wind: wind direction of incoming wind
    :return: Winstral Sx index 
    '''
    grid = np.zeros(dem.shape) * np.nan
    # pad the input edges with np.nan
    view_range = np.ceil(dmax / cellsize).astype(np.int32)
    pad_shape = ((view_range, view_range),) * 2
    dem_padded = np.pad(dem, pad_shape,'constant', constant_values=np.nan)
    # define wind sectors in degrees
    wind_inc = 5
    wind_width = 30
    wind1 = in_wind - wind_width / 2
    wind2 = in_wind + wind_width / 2
    winds = np.arange(in_wind - wind_width / 2, in_wind + wind_width / 2, wind_inc)
    # The angles we check. Add last dimension so we can broadcast the direction
    # samples.
    alpha_rad = np.expand_dims(winds * np.pi / 180, -1)
    # pre-compute the cell indices that are sampled for each direction
    y_offsets = -np.round(np.arange(1, view_range) * np.cos(alpha_rad)).astype(np.int32)
    x_offsets = np.round(np.arange(1, view_range) * np.sin(alpha_rad)).astype(np.int32)
    # pre-compute the distances for each sampled cell
    distances = np.sqrt(cellsize**2 * (x_offsets** 2 + y_offsets** 2))
    # set distances that are too large to np.nan so they're not considered
    distances[(distances == 0.) | (distances > dmax)] = np.nan
    for y in range(view_range, view_range + dem.shape[0]):
        for x in range(view_range, view_range + dem.shape[1]):
            # compute the difference in altitude for all cells along all angles
            altitude_diff = dem_padded[y + y_offsets, x + x_offsets] - dem_padded[y, x]
            # directions are in the first dimension, cells in the last
            slope = altitude_diff / distances
            amax = np.nanmax(slope, -1)
            amin = np.nanmin(slope, -1)
            result = np.where(-amin > amax, amin, amax)
            # maybe nanmean would be more correct, but we reproduce the
            # exisiting implementation for now
            result = np.nanmean(np.arctan(result))
            #result = np.nansum(np.arctan(result)) / len(winds)
            grid[y - view_range, x - view_range] = result
    return grid
