# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:50:15 2014

@author: svfilhol

Geometry, morphology of surface


"""
import numpy as np


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
    

#def aspect(Z):
 #   Zy, Zx = np.gradient(Z)                                                                                                    
  #  S = np.sqrt(Zx**2+Zy**2)            
   # return S
    
    
    