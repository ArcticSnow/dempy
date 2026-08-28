'''
Functions to derive Tabler's profile or Winstral Sx index over a DEM
S. Filhol, March 2026


Tabler's logic is as follow:
1. import DEM
2. rotate dem along wind direction
3. calculate for each line the tabler's profile (use Cython)
4. rotate back new map
5. export to raster (if necessary)


TODO:
- [ ] improve efficiency using cython or rust or numba
- [ ] mak sure tabler's function is complete 
- [ ] add Liu et al 2025 wind factor
- [ ] wrap each own methods into a class that can compute them based on a single DEM open as a dataset

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def resample_line_1m(line, dx):
    '''
    Function to resample a line with a resolution different than 1m.
    :param line: profile elevation, line to be resampled at 1m resolution
    :param dx: resolution of line
    :return: a 1m interpolated version of line
    '''
    xline = np.arange(0, line.__len__() * dx, dx)
    x_interp = np.arange(0, line.__len__() * dx, 1)
    line_int = np.interp(x_interp, xline, line)

    return line_int


def tabler_profile(line, dx=1, slope_adjust=1):

    if dx != 1:
        line = resample_line_1m(line, dx)

    dx = 1
    xline = np.arange(0, line.__len__() * dx, dx)

    # start first point at 45m from dem border
    xstart = 0
    xend = xline.max() - 45 - 44
    z_tab = np.copy(line)

    for xinc in np.arange(xstart, xend, 1):

        t1 = z_tab[xline == (xinc )]
        t2 = z_tab[xline == (xinc + 44)]
        t3 = z_tab[xline == (xinc + 44 + 15)]
        t4 = z_tab[xline == (xinc + 44 + 30)]
        t5 = z_tab[xline == (xinc + 44 + 45)]

        x1 = (t2 - t1) / 45
        x2 = max((t3 - t2) / 15, -0.2)
        x3 = max((t4 - t3) / 15, -0.2)
        x4 = max((t5 - t4) / 15, -0.2)

        # Same but using linear interpolation to estimate slopes:
        #
        # x = xline[(xline>=(xinc-45)) & (xline<xinc)]
        # xs = np.vstack([x, np.ones(len(x))]).T
        # x1 = np.linalg.lstsq(xs, z_tab[(xline>=(xinc-45)) & (xline<xinc)])[0][0]
        # x1 = max(-x1, -0.2)
        #
        # x = xline[(xline > xinc) & (xline <= (xinc+15))]
        # xs = np.vstack([x, np.ones(len(x))]).T
        # x2 = np.linalg.lstsq(xs, z_tab[(xline >= xinc) & (xline < (xinc + 15))])[0][0]
        # x2 = max(-x2, -0.2)
        #
        # x = xline[(xline >= xinc+15) & (xline < (xinc + 30))]
        # xs = np.vstack([x, np.ones(len(x))]).T
        # x3 = np.linalg.lstsq(xs, z_tab[(xline >= xinc+15) & (xline < (xinc + 30))])[0][0]
        # x3 = max(-x3, -0.2)
        #
        # x = xline[(xline >= xinc +30) & (xline < (xinc + 45))]
        # xs = np.vstack([x, np.ones(len(x))]).T
        # x4 = np.linalg.lstsq(xs, z_tab[(xline >= xinc + 30) & (xline < (xinc + 45))])[0][0]
        # x4 = max(-x4, -0.2)

        Nslope = 0.25 * x1 + 0.55 * x2 + 0.15 * x3 + 0.05 * x4
        z_tab[xline == (xinc+45)] = max(line[xline == (xinc+45)], z_tab[xline == (xinc + 44)] + Nslope * slope_adjust * dx)

    print('Tabler profile done')
    return z_tab, xline





if __name__ == "__main__":

    # Example on how to use the Tabler's model

    # Generate a random surface
    line = np.zeros(1000)
    for i in range(1, line.__len__()):
        line[i] = line[i - 1] + np.random.uniform(-1, 1)
    plt.plot(line)
    dx = 1

    # smooth the random surface
    lineM = pd.Series(line).rolling(window=10).mean()
    plt.plot(lineM)
    line = np.array(lineM)

    xline = np.arange(0, line.__len__() * dx, dx)
    # start first point at 45m from dem border

    xstart = 45
    xend = xline.max() - 45
    slope_adjust = 1

    z_tab = np.copy(line)
    for xinc in np.arange(xstart, xend, 1):
        print(xinc)

        t1 = z_tab[xline == (xinc - 44)]
        t2 = z_tab[xline == (xinc - 1)]
        t3 = z_tab[xline == (xinc + 14)]
        t4 = z_tab[xline == (xinc + 29)]
        t5 = z_tab[xline == (xinc + 44)]

        x1 = (t2 - t1) / 45
        x2 = max((t3 - t2) / 15, -0.2)
        x3 = max((t4 - t3) / 15, -0.2)
        x4 = max((t5 - t4) / 15, -0.2)

        Nslope = 0.25 * x1 + 0.55 * x2 + 0.15 * x3 + 0.05 * x4

        z_tab[xline == xinc] = max(line[xline == xinc], z_tab[xline == (xinc - 1)] + Nslope * slope_adjust * dx)

    plt.plot(z_tab, label='z')

        # return z_tab

    # z = tabler_eq(line, .5)
    plt.figure()
    plt.plot(line, label='line')
    plt.plot(z_tab, label='z')
    plt.legend()
