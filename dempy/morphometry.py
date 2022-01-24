# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:50:15 2014

@author: svfilhol

Geometry, morphology of surface


"""
import numpy as np
import cv2


def slope(Z, dx, neighbor=None):
    """
    **S = slope(Z)**\n
    Function to calculate the slope of a 2D matrix

    Parameters
    ==========
    **Z** - 2D matrix of elevation
    **neighbor** - number of neighboring pixel to consider (see numpy.gradient() function help)

    Returns
    =======
    **S** - 2D matrix containing slope values
    """
    if neighbor is None:
        neighbor = 1
    Zy, Zx = np.gradient(Z, neighbor)
    S = np.sqrt((Zx/(dx*neighbor))**2+(Zy/(dx*neighbor))**2) 
    S = np.arctan(S)*180/np.pi           
    return S


def mean_curvature(Z, neighbor=None):
    """
    **K = mean_curvature(Z, neighbor=None)**\n
    Function to calculate the mean curvature of a 2D matrix\n
    see  http://en.wikipedia.org/wiki/Mean_curvature

    Parameters
    ==========
    **Z** - 2D matrix of elevation
    **neighbor** - number of neighboring pixel to consider (see numpy.gradient() function help)

    Returns
    =======
    **K** - 2D matrix of curvature value
    """
    if neighbor is None:
        neighbor=1
    Zy, Zx = np.gradient(Z,neighbor)                                                     
    Zxy, Zxx = np.gradient(Zx,neighbor)                                                  
    Zyy, _ = np.gradient(Zy) 

    H = (Zx**2 + 1)*Zyy - 2*Zx*Zy*Zxy + (Zy**2 + 1)*Zxx
    H = -H/(2*(Zx**2 + Zy**2 + 1)**(1.5))

    return H


def gaussian_curvature(Z, neighbor=None):
    '''
    **K = gaussian_curvature(Z)**\n
    Function to calculate the gaussian curvature of a 2D matrix\n
    see  http://en.wikipedia.org/wiki/Gaussian_curvature 
        
    Parameters
    ==========
    **Z** - 2D matrix of elevation
    **neighbor** - number of neighboring pixel to consider (see numpy.gradient() function help)
    
    Returns
    =======
    **K** - 2D matrix of curvature value
    '''
    if neighbor is None:
        neighbor=1
    Zy, Zx = np.gradient(Z,neighbor)                                                     
    Zxy, Zxx = np.gradient(Zx,neighbor)                                                  
    Zyy, _ = np.gradient(Zy)                                                    
    K = (Zxx * Zyy - (Zxy ** 2)) / (1 + (Zx ** 2) + (Zy **2)) ** 2
    return K
    



def laplacian_of_gaussian_operator(n, sigma, print_sum=False):
    '''
    Function to compute operator Laplacian of Gaussian for a "smooth" curvature computation
    See https://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm

    :param n:       number of pixel. Ideally use an odd number
    :param sigma:   standard deviation to use for gaussian curve. This sets the scale at which 
                    curvature computatino occurs. Make sure n is large enough
    :return:        n*n array operator
    '''

    def lapogaus(x, y, sigma):
        '''
        Laplacian of Gaussian operator for curvature computation smoothing noise.
        https://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm
        
        :param x:       array of x-coordinate. Should be centered on 0.
        :param y:       array of y-coordinate. Should be centered on 0.
        :param sigma:   float, standard deviation setting the with of the bell curve
        :return: array containing the operator
        '''
        op = -1/(np.pi * sigma **4)*(1-((x**2 + y**2)/(2*sigma**2)))*np.exp(-((x**2 + y**2)/(2*sigma**2)))
        return op

    x = np.arange(-n,n)
    y = np.arange(-n,n)
    Xs, Ys = np.meshgrid(x,y)
    kernel = lapogaus(Xs,Ys,sigma)
    if print_sum:
        print('Kernel sum = ',np.sum(kernel))
    return kernel

def plot_kernel(kernel):
    plt.figure()
    plt.imshow(kernel)
    plt.colorbar()

def convolve_array(arr, kernel):
    '''
    function to convolve a 2D array with a kernel
    '''
    res = cv2.filter2D(arr, -1, kernel)
    return res





#def aspect(Z):
 #   Zy, Zx = np.gradient(Z)                                                                                                    
  #  S = np.sqrt(Zx**2+Zy**2)            
   # return S
    
    
    