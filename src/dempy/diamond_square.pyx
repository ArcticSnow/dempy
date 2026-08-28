"""
Cython code for generating random terrain using the Diamond Square algorithm
S. Filhol, March 2026

Installation:

create a file setup.py with the following
    from distutils.core import setup
    from Cython.Build import cythonize
    import numpy

    setup(
	    ext_modules = cythonize(['*.pyx']),
	    include_dirs=[numpy.get_include()]
	    )

Run python setup.py install 

Then you can test the compilation using:


import pyximport; pyximport.install()
import diamond_square as ds

from time import time

t0 = time()
a = ds.diamondSquare(8000, 8000, 1, 0.2)
t1 = time()
print('Processing time: ' + str(t1 - t0))

from matplotlib import pyplot as plt 
plt.imshow(a)
plt.show()




"""

cimport numpy as np
import numpy as np
from libc.math cimport isnan
from cython import wraparound, boundscheck
from libc.stdlib cimport rand, RAND_MAX

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@boundscheck(False)
@wraparound(False)
cdef np.ndarray[double, ndim=2] diamondStep(np.ndarray[double, ndim=2] TR, int tSize, float Zrange):
    cdef int sh = np.int(tSize / 2)
    cdef int maxIndex = TR[:, 0].__len__()-1 # size of TR
    cdef int row = sh
    cdef int col = sh # row, col are the indices of each square 's centerpoint
    cdef float value
    cdef float displacement

    while (row < maxIndex):
        while (col < maxIndex):
            # average heightvalue of 4 cornerpoints
            value = TR[row - sh, col - sh] + \
                    TR[row - sh, col + sh] + \
                    TR[row + sh, col - sh] + \
                    TR[row + sh, col + sh]
            value = value / 4

            # displacement
            displacement = rand()/RAND_MAX * Zrange - Zrange / 2
            value = value + displacement

            # set diamond - point( if not predefined)
            if isnan(TR[row, col]):
                TR[row, col] = value


            # next square in same row
            col = col + tSize
        # next row
        col = sh
        row = row + tSize

    return TR

@boundscheck(False)
@wraparound(False)
cdef np.ndarray[double, ndim=2] squareStep(np.ndarray[double, ndim=2] TR, int tSize, float Zrange):
    cdef int sh = np.int(tSize / 2)
    cdef int maxIndex = TR[:, 0].__len__()-1 # size of TR
    colStart = sh
    cdef int row = 0
    cdef int col = colStart # row, col are the indices of each diamond's centerpoint
    cdef float value
    cdef int nop
    cdef float displacement

    while (row <= maxIndex):
        while (col <= maxIndex):
            value = 0
            nop = 4 # number of points the following cases handle the boundary points,
                    #  i.e.the incomplete diamonds

            # north
            if row > 0:
                value = value + TR[row - sh, col]
            else:
                nop = nop - 1

            # east
            if col < maxIndex:
                value = value + TR[row, col + sh]
            else:
                nop = nop - 1

            # south
            if row < maxIndex:
                value = value + TR[row + sh, col]
            else:
                nop = nop - 1

            # west
            if col > 0:
                value = value + TR[row, col - sh]
            else:
                nop = nop - 1

            # displacement
            displacement = (rand()/RAND_MAX) * Zrange - Zrange / 2
            value = value / nop + displacement

            # set square point( if not predefined)
            if isnan(TR[row, col]):
                TR[row, col] = value

            # next diamond in same row
            col = col + sh

        # next row the starting column alternates between 1 and sh
        if colStart == 0:
            colStart = sh
        else:
            colStart = 0

        col = colStart
        row = row + sh

    return TR


@boundscheck(False)
@wraparound(False)
def diamondSquare(int nrows, int ncols, float Zrange, float Roughness):

    cdef int tSize 
    cdef np.ndarray[double, ndim=2] terrain
    cpdef np.ndarray[double, ndim=2] T

    tSize = np.int(1 + 2 ** (np.ceil(np.log(np.max([ncols, nrows])) / np.log(2))))

    if (Roughness<0) or (Roughness>1):
        raise ValueError('Roughness must be within [0:1]')

    # ==================================================================
    #  Functions
    # ==================================================================


    #==================================================================
    # 1. Initialize terrain
    #==================================================================

    terrain = np.zeros((tSize,tSize)) * np.nan
    terrain[0, 0] = 0
    terrain[0, tSize-1] = 0
    terrain[tSize-1, 0] = 0
    terrain[tSize-1, tSize-1] = 0

    tSize = tSize - 1

    # ==================================================================
    # 2. Loop
    # ==================================================================
    while tSize > 0:
        # perform diamond stepfor entire terrain
        terrain = diamondStep(terrain, tSize, Zrange)

        # perform square step for entire terrain
        terrain = squareStep(terrain, tSize, Zrange)

        # adjust parameters for next scale
        tSize = np.int(tSize / 2)
        Zrange = Zrange * (1 / (2 ** Roughness))

        if tSize == 1:
            break

    T = terrain[0:nrows, 0:ncols] # clip out a matrix of the requested size

    return T
