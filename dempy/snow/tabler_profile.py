
'''
Functions to derive Tabler's profile over a DEM

1. import DEM
2. rotate dem along wind direction
3. calculate for each line the tabler's profile (use Cython)

4. rotate back new map
5. export to raster (if necessary)


LATER: implement tabler and sx coefficien using cython
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

#===============================================
# include functions here





#===============================================



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




