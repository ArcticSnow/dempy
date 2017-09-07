'''
Implementation of the Diamond square algorithm

Python transcription of the Matlab code developped by Perron for

Perron, J.T., J.W. Kirchner and W.E. Dietrich (2008), Spectral signatures of characteristic spatial scales
and non-fractal structure in landscapes, Journal of Geophysical Research - Earth Surface.

'''

import numpy as np

def diamondSquare(nrows, ncols, Zrange, Roughness, seed=False):


    tSize = np.int(1 + 2 ** (np.ceil(np.log(np.max([ncols, nrows])) / np.log(2))))

    if (Roughness<0) or (Roughness>1):
        raise ValueError('Roughness must be within [0:1]')

    if seed is True:
        print 'Seeding needs to be implemented!'


    # ==================================================================
    #  Functions
    # ==================================================================

    def diamondStep(TR, tSize, Zrange):
        sh = np.int(tSize / 2)
        maxIndex = TR[:, 0].__len__()-1 # size of TR
        row = sh
        col = sh # row, col are the indices of each square 's centerpoint

        while (row < maxIndex):
            while (col < maxIndex):
                # average heightvalue of 4 cornerpoints
                value = TR[row - sh, col - sh] + \
                        TR[row - sh, col + sh] + \
                        TR[row + sh, col - sh] + \
                        TR[row + sh, col + sh]
                value = value / 4

                # displacement
                displacement = np.random.rand(1) * Zrange - Zrange / 2
                value = value + displacement

                # set diamond - point( if not predefined)
                if np.isnan(TR[row, col]):
                    TR[row, col] = value


                # next square in same row
                col = col + tSize
            # next row
            col = sh
            row = row + tSize

        return TR

    def squareStep(TR, tSize, Zrange):
        sh = np.int(tSize / 2)
        maxIndex = TR[:, 0].__len__()-1 # size of TR
        colStart = sh
        row = 0
        col = colStart # row, col are the indices of each diamond's centerpoint

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
                displacement = np.random.rand(1) * Zrange - Zrange / 2
                value = value / nop + displacement

                # set square point( if not predefined)
                if np.isnan(TR[row, col]):
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

    T = terrain[0:nrows-1, 0:ncols-1] # clip out a matrix of the requested size

    return T

if __name__ == '__main__':
    diamondSquare(100,100,5,0.5)













































