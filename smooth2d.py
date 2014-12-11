__author__ = 'svfilhol'

import matplotlib.pyplot as plt
import numpy as np
import cv2


def kernel_square(nPix):
    """
    Function to defin a square kernel of equal value for performing averaging
    :param nPix: size of the kernel in pixel
    :return: kernel matrix
    """
    print "Averaging kernel of " + str(nPix) + " by " + str(nPix)
    kernel = np.empty([nPix, nPix])
    kernel.fill(1)
    kernel /= kernel.sum()   # kernel should sum to 1!  :)
    return kernel

def smooth(mat, kernel):
    """
    Function that produce a smothed version of the 2D array
    :param mat: Array to smooth
    :param kernel: kernal array (output) from the function kernel_square()
    :return: smoothed array
    """
    r = cv2.filter2D(mat, -1, kernel)
    print "Done ..."
    return r

def crop_frame(mat, nPix):
    """
    Function to crop a 2D array by a frame all around the array by nPix pixels
    :param mat: 2D array to crop
    :param nPix: width of the crop in pixel
    :return: cropped array
    """
    matShape = mat.shape()
    matout = mat[(0+nPix):(matShape[0]-nPix), (0+nPix):(matShape[1]-nPix)]
    return matout

def compare_smooth(mat, nCompare = None, kPixMax = None):
    """
    Function to produce an array of smoothed DEM. Result can be displayed with plot_compare()
    :param mat:
    :param nCompare:
    :param kPixMax:
    :return: 3D array containing the smoothed array
    """
    if nCompare is None:
        nCompare = 4
    if kPixMax is None:
        kPixMax = 500
    kPix = np.round(np.linspace(10,kPixMax,nCompare))
    matSmoo = np.ones([mat.shape[0],mat.shape[1],nCompare])
    for i in range(0, nCompare):
        kernel = kernel_square(kPix[i])
        matSmoo[:,:,i] = smooth(mat, kernel)
    return matSmoo

def plot_compare(matSmoo, ncol, nrow, nCompare = None, kPixMax = None, residual = False, matRef = None):
    """
    Function to plot results from compare_smooth()
    :param matSmoo:
    :param ncol:
    :param nrow:
    :param nCompare:
    :param kPixMax:
    :param residual:
    :param matRef:
    :return: a plot
    """
    if nCompare is None:
        nCompare = 4

    if kPixMax is None:
        kPixMax = 500
    kPix = np.round(np.linspace(10, kPixMax, nCompare))

    if residual is True and matRef is not None:
        for i in range(0, nCompare):
            matSmoo[:,:,i] = matRef - matSmoo[:,:,i]
    plt.figure()

    for i in range(0, nCompare):
        plt.subplot(ncol,nrow,i)
        plt.imshow(matSmoo[:,:,i])
        plt.colorbar()
        plt.title("Kernel size of " + str(kPix[i]))
    plt.show()
