import numpy as np
import matplotlib.pyplot as plt
import cv2


class perron_fftDEM(object):
    def __init__(self):
        self.dem = None
        self.dx = None
        self.dy = self.dx
        self.nx = None
        self.ny = None


    def fftdem(self, pad=False, pad_window=True, openCV=True):
        self.nx, self.ny = self.dem.shape
        if pad_window is True:
            img_pad, Wss = self.hann2d(ret=True)
        else:
            Wss = np.sum(np.sum(np.ones((self.ny, self.nx))))

        if pad is True:
            print "Needs to be implemented"
        else:
            Lx = self.nx
            Ly = Lx

        if (pad is False) and (pad_window is False):
            Lx = 2 ** (np.ceil(np.log(np.max(self.nx, self.ny)) / np.log(2)))
            Ly = Lx
            print "image must be padded"

        # Frequency increments: from zero to Nyquist freq 1(2*dx)
        dfx = 1/(self.dx * Lx)
        dfy = 1/(self.dy * Ly)

        # calculate the 2D FFT
        if openCV:
            fft = cv2.dft(np.float32(self.dem), flags=cv2.DFT_COMPLEX_OUTPUT)
        else:
            fft = np.fft.fft2(self.dem)

        self.fshift = np.fft.fftshift(fft)

        # Making sure the fft of the dem is detrented as the padding might add a bias from the previously detrended dem
        self.fshift[Ly / 2 + 1, Lx / 2 + 1] = 0

        # derive DFT periodogram
        self.DFTperiodogram = self.fshift * np.conj(self.fshift) / (Lx * Ly * Wss)

        # matrix of radial frequencies
        xc = Lx / 2 + 1
        yc = Ly / 2 + 1
        cols, rows = np.meshgrid(yc, xc)
        self.FreqMat = np.sqrt((dfy * (rows - yc)) ** 2 + (dfx * (cols - xc)) ** 2)

        # vector of sorted frequency and power
        fft_part = self.DFTperiodogram[:, 1:(Lx / 2 + 1)]
        fvec = self.FreqMat[:, 1:(Lx / 2 + 1)]
        fvec[(yc + 1):Ly, xc] = -1
        fvec = np.sort(np.concatenate(fvec, fft_part, 0), 0)
        fvec = fvec[fvec[:, 1] > 0, :]

        # separate into power and frequency vectors
        self.Power_vec = 2 * fvec[:, 2]
        self.Freq_vec = fvec[:, 1]

        if ret is True:
            return self.Power_vec, self.Freq_vec, self.FreqMat, self.DFTperiodogram

    def hann2d(self, ret=False):

        # matrix coordinates of centroid
        a = (self.nx + 1) / 2
        b = (self.ny + 1) / 2

        X, Y = np.meshgrid(self.ny, self.nx)

        theta = (X == a) * (np.pi / 2) + (X != a) * np.arctan2((Y - b), (X - a))  # angular polar coordinate

        r = np.sqrt((Y - b) ** 2 + (X - a) ** 2)  # radial polar coordinates

        # radius of ellipse for this theta
        rprime = np.sqrt((a ** 2) * (b ** 2) * (b ** 2 * (np.cos(theta)) ** 2 + a ** 2 * (np.sin(theta)) ** 2) ** (-1))

        hanncoeff = (r < rprime) * (0.5 * (1 + np.cos(np.pi * r / rprime)))
        H = self.dem * hanncoeff

        Wss = np.sum(np.sum(hanncoeff ** 2))

        if ret is True:
            return H, Wss

    def plot_2D_spec(self):
        return

    def plot_1D_spec(self):

        # Add code to derive bining version of the periodogram
        print 'code for bining!'
        # Create figure
        plt.figure()
        plt.scatter(self.Freq_vec, self.Power_vec)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Frequency')
        plt.ylabel('DFT mean squared amplitude')
        plt.title('Periodogram')







p = perron_fftDEM()
p.dem = t
p.dx = 1
p.dy=1
p.fftdem()
