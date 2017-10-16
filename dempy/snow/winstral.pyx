from __future__ import division
cimport numpy as np
import numpy as np

from cython import wraparound, boundscheck

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


cdef float winstralSX_core2(np.ndarray[double, ndim=2] dem, float cellsize,float dmax, float wind, int i, int j):
    # convert wind angle to radians
    cdef float alpha_rad = wind * np.pi / 180

    # extract altitude of reference from the DEM
    cdef float ref_altitude = dem[i, j]

    # intialize variables for while-loop
    cdef int ii = i
    cdef int jj = j
    cdef int nrows, ncols 
    nrows = dem.shape[0]
    ncols = dem.shape[1]
    cdef float max_tan_sx = 0
    cdef int ll = ii
    cdef int mm = jj
    cdef int nb_cells = 0

    cdef float inv_distance
    cdef float delta_elev 
    cdef float altitude 
    cdef float tan_sx
    cdef float max_sx 

    # loop algorithm for every pixels
    while (ll > 0) and (ll < nrows - 1) and (mm > 0) and (mm < ncols - 1):
        altitude = dem[ll, mm]
        if np.isnan(altitude):
            break
        if (ll != ii) or (mm != jj):

            # Get elevation change between the two points
            delta_elev = altitude - ref_altitude

            # calculate distance
            inv_distance = 1 / (np.sqrt(cellsize ** 2 * ((ll - ii) ** 2 + (mm - jj) ** 2)))

            if inv_distance < (1 / dmax):
                break

            # estimate tangent
            tan_sx = delta_elev * inv_distance

            if np.abs(tan_sx) > np.abs(max_tan_sx):
                max_tan_sx = tan_sx

        nb_cells += 1
        ll = ii - np.round(nb_cells * np.cos(alpha_rad))
        mm = jj + np.round(nb_cells * np.sin(alpha_rad))
    max_sx = np.arctan(max_tan_sx)
    return max_sx

cdef float winstralSX_core1(np.ndarray[double, ndim=2]  dem, float cellsize, float dmax, float wind1, float wind2, float wind_inc, int ii, int jj):
    cdef float b

    # to make sure wind1 is smaller than wind2
    if wind1 > wind2:
        b = wind1
        wind1 = wind2
        wind2 = b

    # intialize counting variable
    cdef float mysum = 0
    cdef int mycount = 0

    cdef float max_sx
    cdef float sx

    # perform Sx-algorithm for every increment of wind direction
    for wind in np.arange(wind1, wind2, wind_inc):
        max_sx = winstralSX_core2(dem, cellsize, dmax, wind, ii, jj)
        mysum += max_sx
        mycount += 1
    if mycount == 0:
        sx = np.nan
    else:
        sx = mysum/mycount
    return sx



cpdef np.ndarray[double, ndim=2] winstralSX(np.ndarray[double, ndim=2] dem, float cellsize, float dmax, float in_wind):
    cdef np.ndarray[double, ndim=2] grid
    cdef int wind_inc = 5
    cdef int wind_width = 30
    cdef int nrows
    cdef int ncols
    cdef int ii
    cdef int jj

    nrows = dem.shape[0]
    ncols = dem.shape[1]


    # initialize a 2D matrice of the same size as the DEM
    grid = np.zeros((nrows,ncols)) - 9999

    # define wind sectors in degrees
    
    
    cdef float wind1 = in_wind - wind_width/2
    cdef float wind2 = in_wind + wind_width/2

    # loop algorithm for every pixel of the DEM
    for ii in np.arange(1,nrows,1):
        for jj in np.arange(1,ncols,1):
            grid[ii, jj] = winstralSX_core1(dem, cellsize, dmax, wind1, wind2, wind_inc, ii, jj)

    grid[grid == (-9999)] = np.nan
    return grid
