from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import cv2
import pandas as pd
import decomposition.diamondSquare as ds
from scipy.stats import chi2

class perron_fftDEM(object):
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

    def fftmat(self, mat,dx=1,dy=1, pad=False, pad_window=True, openCV=False, ret=False):
        nx, ny = mat.shape
        if pad_window is True:
            img_pad, Wss = self.hann2d(mat)
        else:
            Wss = np.ones((ny, nx)).sum()

        if pad is True:
            print "Needs to be implemented"
        else:
            Lx = np.int(nx)
            Ly = Lx

        if (pad is False) and (pad_window is False):
            Lx = np.int(2 ** (np.ceil(np.log(np.max([nx, ny])) / np.log(2))))
            Ly = Lx

            print "image must be padded"

        # Frequency increments: from zero to Nyquist freq 1(2*dx)
        dfx = 1/(dx * Lx)
        dfy = 1/(dy * Ly)


        # calculate the 2D FFT
        if openCV:
            fft = cv2.dft(np.float32(mat), flags=cv2.DFT_COMPLEX_OUTPUT)
            fshift = np.fft.fftshift(fft)
            DFTperiodogram = np.copy(cv2.magnitude(fshift[:,:,0],fshift[:,:,1])**2)

        else:
            fft = np.fft.fft2(mat)
            fshift = np.fft.fftshift(fft)

            # Making sure the fft of the dem is detrented as the padding might add a bias from the previously detrended dem
            #fshift[np.int(Ly / 2), np.int(Lx / 2)] = 0

            # derive DFT periodogram
            #DFTperiodogram = np.copy(fshift * np.conj(fshift) / (Lx * Ly * Wss))
            DFTperiodogram = np.copy(np.abs(fshift) ** 2)

        # matrix of radial frequencies
        xc = np.int(Lx / 2)
        yc = np.int(Ly / 2)
        cols, rows = np.meshgrid(np.arange(0, Lx), np.arange(0, Ly))
        FreqMat = np.sqrt((dfy * (rows - yc)) ** 2 + (dfx * (cols - xc)) ** 2)

        # vector of sorted frequency and power
        fft_part = np.copy(DFTperiodogram[:, 0:np.int(Lx / 2)])
        fvec = np.copy(FreqMat[:, 0:np.int(Lx / 2)])
        fvec[yc:Ly-1, xc-1] = -1

        fvec = np.vstack((fvec.flatten(), fft_part.flatten())).T
        fvec= fvec[fvec[:, 0].argsort(),]
        fvec = np.copy(fvec[fvec[:, 0] > 0, :])

        # separate into power and frequency vectors
        Power_vec = 2 * fvec[:, 1]
        Freq_vec = fvec[:, 0]

        return Power_vec, Freq_vec, FreqMat, DFTperiodogram

    def fftdem(self, openCV=True):
        self.nx, self.ny = self.dem.shape
        self.Power_vec, self.Freq_vec, self.FreqMat, self.DFTperiodogram = self.fftmat(self.dem, self.dx, self.dy, openCV=openCV)

    def hann2d(self, mat):
        nx, ny = mat.shape

        # matrix coordinates of centroid
        a = (nx + 1) / 2
        b = (ny + 1) / 2

        X, Y = np.meshgrid(np.arange(0, ny), np.arange(0, nx))

        theta = (X == a) * (np.pi / 2) + (X != a) * np.arctan2((Y - b), (X - a))  # angular polar coordinate

        r = np.sqrt((Y - b) ** 2 + (X - a) ** 2)  # radial polar coordinates

        # radius of ellipse for this theta
        rprime = np.sqrt((a ** 2) * (b ** 2) * (b ** 2 * (np.cos(theta)) ** 2 + a ** 2 * (np.sin(theta)) ** 2) ** (-1))

        hanncoeff = (r < rprime) * (0.5 * (1 + np.cos(np.pi * r / rprime)))
        H = mat * hanncoeff

        Wss = (hanncoeff ** 2).sum()

        return H, Wss

    def plot_2D_spec(self):
        return

    def bin_scatter(self,x,y, nbins=10):

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

    def plot_1D_spec(self,x, y, nbins=10, errorbar=False):

        # Create figure
        plt.figure()
        plt.scatter(x, y)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Frequency [1/m]')
        plt.ylabel('DFT mean squared amplitude')
        plt.title('Periodogram')

        scat_bin = self.bin_scatter(x,y,nbins=nbins)
        plt.plot(scat_bin.freq, scat_bin.power_mean,color='k')
        plt.scatter(scat_bin.freq, scat_bin.power_mean, s=50, c='k')
        if errorbar:
            plt.errorbar(scat_bin.freq, scat_bin.power_mean, yerr=scat_bin.power_std/2, color='k')

    def plot_1D_specDEM(self, nbins=10, errorbars=False):
        self.plot_1D_spec(self.Freq_vec, self.Power_vec, nbins=nbins, errorbar=errorbars)

    def plot_1D_specNORM(self, nbins=10, errorbars=False):
        self.plot_1D_spec(self.Freq_vec, self.Power_vec_norm, nbins=nbins, errorbar=errorbars)

    def normalized_spect(self, H, Zrange=1, nSynth=20, demVar=1):
        for i in range(1, nSynth+1):
            print 'Synthetic surface # ' + str(i) + ' of ' + str(nSynth)
            synthDEM = ds.diamondSquare(self.nx, self.ny, Zrange, H)
            synthDEM = synthDEM * np.sqrt(demVar)/np.std(synthDEM)

            Pvec, fvec, freqmat, Pmat = self.fftmat(synthDEM, dx=self.dx, dy=self.dy, pad_window=True)
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
        print self.synth_rmse

if __name__ == '__main__':

    # script to demonstrate how to use perron_fftDEM() based on a random surface
    t = ds.diamondSquare(100, 50, 20, .3)
    p = perron_fftDEM()
    p.dem = t
    p.dx = 1
    p.dy = 1
    p.fftdem()

    # To generate synthetic spectrum, use H that minimizes the RMSE (could be done incrementaly)
    p.normalized_spect(H=0.3, demVar=p.Power_vec.sum())

    # To look at confidence level of for the normalized periodogram, use a chi-square distribution
    plt.imshow(p.DFTperiodogram_norm>chi2.pdf(.90, 2))

    p.plot_1D_specDEM(30)
    p.plot_1D_spec(nbins=30)

    plt.figure()
    plt.plot(p.Power_vec)

