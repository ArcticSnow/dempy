from __future__ import division

import cv2
import numpy as np
import matplotlib.pyplot as plt
import smooth2d as sm

'''

#==========================================================================================
# 1. load DEM and pad with zeros
# 2. compute fft2d
# 3.

S. Filhol, Sepember 2016


# Ressources:
http://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_imgproc/py_transforms/py_fourier_transform/py_fourier_transform.html
https://www.cs.auckland.ac.nz/courses/compsci773s1c/lectures/ImageProcessing-html/topic1.htm
'''

class fftImage(object):
    def __init__(self):
        self.img = None
        self.rows, self.cols = img.shape
        self.dx = 1
        self.Fnyquist = 1/(2*self.dx)

    def smoothing(self, Kernel_size=None, ret=False):
        if Kernel_size is None:
            Kernel_size = np.int(np.min(self.img.shape)/3)
        kernel = sm.kernel_square(Kernel_size)
        self.img_sm = sm.smooth(self.img, kernel)
        self.img = self.img - self.img_sm
        if ret is True:
            return self.img_sm

    def window_pad(self, ret=False):
        w1 = np.cos(np.linspace(-np.pi / 2, np.pi / 2, cols))
        w2 = np.cos(np.linspace(-np.pi / 2, np.pi / 2, rows))
        W = np.outer(w2, w1)
        self.img_pad = self.img * W
        if ret is True:
            return self.img_pad

    def pad(self, ret=False):
        nrows = cv2.getOptimalDFTSize(rows)
        ncols = cv2.getOptimalDFTSize(cols)
        right = ncols - cols
        bottom = nrows - rows
        bordertype = cv2.BORDER_CONSTANT  # just to avoid line breakup in PDF file
        self.img_pad = cv2.copyMakeBorder(img, 0, bottom, 0, right, bordertype, value=0)

        if ret is True:
            return self.img_pad


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
        y, x = np.indices(self.power_spectrum.shape)

        if not center:
            center = np.array([(x.max() - x.min()) / 2.0, (x.max() - x.min()) / 2.0])

        r = np.hypot(x - center[0], y - center[1])

        # Get sorted radii
        ind = np.argsort(r.flat)
        r_sorted = r.flat[ind]
        i_sorted = self.power_spectrum.flat[ind]

        # Get the integer part of the radii (bin size = 1)
        r_int = r_sorted.astype(int)

        # Find all pixels that fall within each radial bin.
        deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
        rind = np.where(deltar)[0]  # location of changed radius
        nr = rind[1:] - rind[:-1]  # number of radius bin

        # Cumulative sum to figure out sums for each radius bin
        csim = np.cumsum(i_sorted, dtype=float)
        tbin = csim[rind[1:]] - csim[rind[:-1]]

        self.radial_prof = tbin / nr

        if ret is True:
            return self.radial_prof

    def power_spectrum(self, image, ret=False, openCV=True):

        if openCV:
            f = cv2.dft(np.float32(image),flags = cv2.DFT_COMPLEX_OUTPUT)
        else:
            f = np.fft.fft2(image)
        fshift = np.fft.fftshift(f)
        # magnitude_spectrum = 20*np.log(np.abs(fshift))
        self.power_spectrum = np.abs(fshift) ** 2

        if ret is True:
            return self.power_spectrum


    def plot_1D_spectrum(self):

        ps1d = self.azimuthalAverage(ret=True)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()

        ax1.semilogy((1/self.dx)*(np.linspace(0,ps1d.__len__()/2-1,ps1d.__len__()))/ps1d.__len__(), ps1d)

        ax1Ticks = ax1.get_xticks()
        ax2Ticks = ax1Ticks

        def tick_function(X):
            c = 3.e2
            V = 1 / X
            return ["%.3f" % z for z in V]

        ax2.set_xticks(ax2Ticks)
        ax2.set_xbound(ax1.get_xbound())
        ax2.set_xticklabels(tick_function(ax2Ticks))

        ax1.set_ylabel('Power Spectrum')
        ax1.set_xlabel('Frequency')
        ax2.set_xlabel('Wavelength (m)')

        ax1.grid(True)






