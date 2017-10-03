from __future__ import division
import pyximport; pyximport.install()
import ds
import numpy as np
import matplotlib.pyplot as plt
import cv2
import pandas as pd
from scipy.stats import chi2
import smooth2d as sm

class fftDecomposition(object):
    '''
    Class to perform analysis such as the one in Perron et al. 2008

    '''
    def __init__(self):
        self.dem = None
        self.dx = None
        self.dy = self.dx
        self.nx = None
        self.ny = None
        self.FreqMat = None
        self.DFTperiodogram = None
        self.Power_vec = None
        self.fshift = None
        self.Freq_vec = None
        self.DFTperiodogram_norm = None
        self.Power_vec_norm = None
        self.radial_power = None
        self.radial_freq = None
        self.dem_pad = None

    def smoothing(self, kernel_size=None, ret=False):
        '''
        Function to smooth
        :param Kernel_size:
        :param ret:
        :return:
        '''
        if kernel_size is None:
            kernel_size = np.int(np.min(self.dem.shape)/3)
        kernel = sm.kernel_square(kernel_size)
        self.dem_sm = sm.smooth(self.dem, kernel)
        self.dem = self.dem - self.dem_sm
        if ret is True:
            return self.dem_sm

    def fftmat(self, mat, dx=1, dy=1, pad=False, pad_window=True, openCV=True):
        '''
        Function to derive Fourier transform. Returns 2D and 1D periodogram
        :param mat:         input 2D matirx
        :param dx:          pixel resolution of the matrix in the x-dir
        :param dy:          pixel resolution of the matrix in the y-dir
        :param pad:         apply zeros padding around the matrix following. Default is False
        :param pad_window:  apply hanning padding
        :param openCV:      use openCV algorithm to calculate FFT. Faster than numpy
        :return:            Power_vec       -- 1D vector of the power from FFT
                            Freq_vec,       -- 1D vector of the Frequency
                            FreqMat         -- 2D matrix of frequency
                            DFTperiodogram  -- 2D matrix of periodogram
        '''
        nx, ny = mat.shape

        # apply padding if asked
        if pad_window is True:
            img_pad = self.hann2d(mat)

        if pad is True:
            print "Needs to be implemented"
            Lx = np.int(2 ** (np.ceil(np.log(np.max([nx, ny])) / np.log(2))))
            Ly = Lx
            img_pad = mat
        else:
            Lx = np.int(nx)
            Ly = Lx

        if (pad is False) and (pad_window is False):
            Lx = np.int(nx)
            Ly = Lx
            img_pad = mat

            print "image must be padded"

        # Frequency increments: from zero to Nyquist freq 1(2*dx)
        dfx = 1/(dx * Lx)
        dfy = 1/(dy * Ly)

        print 'Lx=' +str(Lx)
        print 'Ly=' + str(Ly)
        print img_pad.shape


        # calculate the 2D FFT
        if openCV:
            fft = cv2.dft(np.float32(img_pad), flags=cv2.DFT_COMPLEX_OUTPUT)
            fshift = np.fft.fftshift(fft)
            DFTperiodogram = np.copy(cv2.magnitude(fshift[:, :, 0], fshift[:, :, 1])**2)
            fft, fshift = None, None

        else:
            fft = np.fft.fft2(mat)
            fshift = np.fft.fftshift(fft)

            # Making sure the fft of the dem is detrented as the padding might add a bias from the previously detrended dem
            #fshift[np.int(Ly / 2), np.int(Lx / 2)] = 0

            # derive DFT periodogram
            #DFTperiodogram = np.copy(fshift * np.conj(fshift) / (Lx * Ly * Wss))
            DFTperiodogram = np.copy(np.abs(fshift) ** 2)
            fft, fshift = None, None

        # matrix of radial frequencies
        xc = np.int(Lx / 2)
        yc = np.int(Ly / 2)
        cols, rows = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
        FreqMat = np.sqrt((dfy * (rows - yc)) ** 2 + (dfx * (cols - xc)) ** 2)
        rows, cols = None, None

        Freq_vec = np.copy(FreqMat[xc:, yc:].flatten())
        Power_vec = np.copy(DFTperiodogram[xc:, yc:].flatten())

        # vector of sorted frequency and power
        # redesign this part!!!!!!!!!!!!!!!!


        # fft_part = np.copy(DFTperiodogram[:, 0:np.int(Lx / 2)])
        # fmat = np.copy(FreqMat[:, 0:np.int(Lx / 2)])
        # fmat[yc:Ly-1, xc-1] = -1
        #
        #
        # fvec = np.vstack((fmat.flatten(1), fft_part.flatten(1)))
        # print 'fvec shape: ' + str(fvec.shape)
        # fvec= np.copy(fvec[:, fvec[0, :].argsort()])
        # fvec = np.copy(fvec[:, (fvec[0, :] > 0)])
        #
        # # separate into power and frequency vectors
        # Power_vec = 2 * fvec[1, :]
        # Freq_vec = fvec[0, :]

        return Power_vec, Freq_vec, FreqMat, DFTperiodogram, img_pad

    def fftdem(self, pad=False, pad_window=True, openCV=True):
        '''
        Function to perforn fftmat() on the dem loaded into the class
        :param openCV: use
        :return: update class self variable
        '''
        self.nx, self.ny = self.dem.shape
        self.Power_vec, self.Freq_vec, self.FreqMat, self.DFTperiodogram, self.dem_pad = self.fftmat(self.dem, self.dx, self.dy, pad=pad, pad_window=pad_window, openCV=openCV)

    def hann2d(self, mat):
        '''
        Perform hanning filtering to a 2D matrix
        :param mat: 2D input matrix
        :return:    H -- 2D matrix including the filter
                    Wss --
        '''
        nx, ny = mat.shape

        # matrix coordinates of centroid
        a = (nx + 1) / 2
        b = (ny + 1) / 2

        X, Y = np.meshgrid(np.arange(0, ny), np.arange(0, nx))

        inter = a * (np.pi / 2) + (X != a)
        theta = inter * np.arctan2((Y - b), (X - a))  # angular polar coordinate
        inter = None

        r = np.sqrt((Y - b) ** 2 + (X - a) ** 2)  # radial polar coordinates
        X, Y = None, None


        # radius of ellipse for this theta
        rprime = np.sqrt((a ** 2) * (b ** 2) * (b ** 2 * (np.cos(theta)) ** 2 + a ** 2 * (np.sin(theta)) ** 2) ** (-1))
        theta = None

        ind = (r < rprime) * 0.5
        rrp = r / rprime
        r, rprime = None, None
        hanncoeff = ind * (1 + np.cos(np.pi * rrp))
        ind = None
        H = mat * hanncoeff
        hanncoeff, rpp = None, None

        return H

    def plot_2D_spec(self):
        return

    def bin_scatter(self,x,y, nbins=10):
        '''
        Function to bin Y data base on X
        :param x: 1D vector
        :param y:  1D vector
        :param nbins: number of evenly spaced bins
        :return: a dataframe with various stats of Y for each bin
        '''

        df = pd.DataFrame(np.transpose([x,y]))
        df.columns = ['freq', 'power']
        freq_cut, bins = pd.cut(df.freq, nbins, retbins=True)
        binned = df.groupby(freq_cut)

        spect1D_bin = pd.DataFrame()
        spect1D_bin['freq'] = [(a + b) / 2 for a, b in zip(bins[:-1], bins[1:])]
        spect1D_bin['power_min'] = binned.power.min().as_matrix()
        spect1D_bin['power_max'] = binned.power.max().as_matrix()
        spect1D_bin['power_mean'] = binned.power.mean().as_matrix()
        spect1D_bin['power_std'] = binned.power.std().as_matrix()

        return spect1D_bin

    def plot_1D_spec(self,x, y, nbins=10, errorbar=False, bin_only=True):

        # Create figure
        plt.figure()
        if bin_only is False:
            plt.scatter(x, y)

        scat_bin = self.bin_scatter(x,y,nbins=nbins)
        plt.plot(scat_bin.freq, scat_bin.power_mean,color='k')
        plt.scatter(scat_bin.freq, scat_bin.power_mean, s=50, c='k')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Frequency [1/m]')
        plt.ylabel('DFT mean squared amplitude')
        plt.title('Periodogram')
        if errorbar:
            plt.errorbar(scat_bin.freq, scat_bin.power_mean, yerr=scat_bin.power_std/2, color='k')

    def plot_1D_specDEM(self, nbins=10, errorbars=False):
        self.plot_1D_spec(self.Freq_vec, self.Power_vec, nbins=nbins, errorbar=errorbars)

    def plot_1D_specNORM(self, nbins=10, errorbars=False):
        self.plot_1D_spec(self.Freq_vec, self.Power_vec_norm, nbins=nbins, errorbar=errorbars)
        plt.title('Normalized periodogram')

    def normalized_spect(self, H, Zrange=1, nSynth=20, demVar=1):
        '''
        Function to perform the neormalization as in Perron et al. 2008 using the diamond square algorithm
        :param H: roughness parameter. 0<H<1
        :param Zrange: Elevation range of the synthetic terrain model
        :param nSynth: number of terrain model to simulate for deriving normalization
        :param demVar: Variance of the original DEM
        :return: return within the class variable self.DFTperiodogram_norm, and self.Power_vec_norm
        '''
        for i in range(1, nSynth+1):
            print 'Synthetic surface # ' + str(i) + ' of ' + str(nSynth)
            synthDEM = ds.diamondSquare(self.nx, self.ny, Zrange, H)
            synthDEM = synthDEM * np.sqrt(demVar)/np.std(synthDEM)

            Pvec, fvec, freqmat, Pmat, dem_pad = self.fftmat(synthDEM, dx=self.dx, dy=self.dy, pad_window=True)


            if i == 1:
                P = Pvec
                Pm = Pmat
            else:
                P = P + Pvec
                Pm = Pm + Pmat

        P = P / nSynth
        Pm = Pm / nSynth

        # scale the average spectra so they have total power = var. This step is
        # necessary because the average of N spectra, each of which has total power
        # X, will not necessarily have totalpower = X.
        P = P*demVar/np.nansum(P)
        Pm = Pm*demVar/np.nansum(Pm)

        if self.DFTperiodogram is None:
            raise ValueError('run fftdem() first to evaluate DEMs transform')

        self.DFTperiodogram_norm = self.DFTperiodogram / Pm
        self.Power_vec_norm = self.Power_vec / P

        # estimate misfit between synthetic dems and original one. This RMSE should be minize to by tuning H
        bin_spec = self.bin_scatter(self.Freq_vec, self.Power_vec, nbins=15)
        bin_spec_normed = self.bin_scatter(self.Freq_vec, P, nbins=15)

        self.synth_rmse = np.sqrt(np.mean((bin_spec.power_mean - bin_spec_normed.power_mean)**2))
        print 'rmse = ' + str(self.synth_rmse)

    def azimuthalAverage(self, center=None, ret=False):

        '''
        Calculate the azimuthally averaged radial profile.

        self.magnitude_spectrum - The 2D self.magnitude_spectrum
        center - The [x,y] pixel coordinates used as the center. The default is
                 None, which then uses the center of the self.magnitude_spectrum (including
                 fracitonal pixels).

        code from: http://www.astrobetter.com/wiki/tiki-index.php?page=python_radial_profiles
        '''
        # Calculate the indices from the self.magnitude_spectrum
        if self.DFTperiodogram is None:
            raise ValueError('You must run self.fftdem()')
        y, x = np.indices(self.DFTperiodogram.shape)

        if not center:
            center = np.array([(x.max() + 1 - x.min()) / 2.0, (x.max()+1 - x.min()) / 2.0])

        print center
        r = np.hypot(x - center[0], y - center[1])
        x, y =None, None

        # Get sorted radii
        ind = np.argsort(r.flat)
        r_sorted = r.flat[ind]
        i_sorted = self.DFTperiodogram.flat[ind]
        ind = None

        # Get the integer part of the radii (bin size = 1)
        r_int = r_sorted.astype(int)
        r_sorted = None

        # Find all pixels that fall within each radial bin.
        deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
        rind = np.where(deltar)[0]  # location of changed radius
        deltar = None

        nr = rind[1:] - rind[:-1]  # number of radius bin

        # Cumulative sum to figure out sums for each radius bin
        csim = np.cumsum(i_sorted, dtype=float)
        tbin = csim[rind[1:]] - csim[rind[:-1]]
        csim, rind, i_sorted = None, None, None

        self.radial_power = tbin / nr
        #self.radial_freq = (1/self.dx)*(np.linspace(0, self.radial_power.__len__()/2-1, self.radial_power.__len__()))
        tbin, nr = None, None

        Lx, Ly = self.dem_pad.shape
        xc = np.int(Lx / 2)
        yc = np.int(Ly / 2)
        dfx = 1 / (self.dx * Lx)
        dfy = 1 / (self.dy * Ly)
        self.radial_freq = np.linspace(0, np.sqrt((dfy * (Lx - yc)) ** 2 + (dfx * (Ly - xc)) ** 2), self.radial_power.__len__())

        if ret is True:
            return self.radial_power, self.radial_freq

if __name__ == '__main__':

    # script to demonstrate how to use fftDecomposition() based on a random surface
    t = ds.diamondSquare(100, 50, 20, .3)
    p = fftDecomposition()
    p.dem = t
    p.dx = 1
    p.dy = 1
    p.smoothing(kernel_size=20)
    p.fftdem()

    # To generate synthetic spectrum, use H that minimizes the RMSE (could be done incrementaly)
    p.normalized_spect(H=0.3, demVar=p.Power_vec.sum())

    # To look at confidence level of for the normalized periodogram, use a chi-square distribution
    plt.imshow(p.DFTperiodogram_norm>chi2.pdf(.90, 2))

    p.plot_1D_specDEM(30)
    p.plot_1D_specNORM(nbins=30)


    p.azimuthalAverage()
    plt.figure()
    plt.semilogy(p.radial_freq, p.radial_power)
    plt.xlabel('Frequency [1/m]')
    plt.ylabel('Power Spectrum')

    plt.figure()
    plt.imshow(p.dem)
