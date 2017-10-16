from __future__ import division
cimport numpy as np
import numpy as np

from cython import wraparound, boundscheck

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t


@boundscheck(False)
@wraparound(False)
cdef np.ndarray[double, ndim=1] resample_line_1m(np.ndarray[double, ndim=1]  line, float dx):
    '''
    Function to resample a line with a resolution different than 1m.
    :param line: profile elevation, line to be resampled at 1m resolution
    :param dx: resolution of line
    :return: a 1m interpolated versin of line
    '''

    cdef np.ndarray[double, ndim=1] xline = np.arange(0, line.__len__() * dx, dx)
    cdef np.ndarray[double, ndim=1] x_interp = np.arange(0, line.__len__() * dx, 1)
    cdef np.ndarray[double, ndim=1] line_int = np.interp(x_interp, xline, line)

    return line_int


@boundscheck(False)
@wraparound(False)
cpdef np.ndarray[double, ndim=2]  tabler_profile(np.ndarray[double, ndim=1] line, float dx=1, float slope_adjust=1):

    if dx != 1:
        line = resample_line_1m(line, dx)

    dx = 1
    cdef np.ndarray[double, ndim=1]  xline = np.arange(0, line.__len__() * dx, dx)

    # start first point at 45m from dem border
    cdef int xstart = 0
    cdef int xend = xline.max() - 45 - 44
    cdef np.ndarray[double, ndim=1]  z_tab = np.copy(line)
    cdef int xinc

    cdef float t1, t2, t3, t4, t5, x1, x2, x3, x4, Nslope

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

    
    cdef np.ndarray[double, ndim=2] ret = np.vstack((xline, z_tab)) 
    
    print('Tabler profile done')
    return ret
