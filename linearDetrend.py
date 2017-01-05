# -*- coding: utf-8 -*-
"""
Created on Wed May 14 18:50:20 2014

@author: simonfilhol
"""
import numpy as np

def detrend(M):
    '''
    **D, coeff = detrend(M)**
    
    fits a plane to the surface defined by the elements of matrix M and detrends the surface by subtracting the value of the planar fit at each element. Returns the detrended surface in matrix D.

    Dependencies: lsplane()

    by S. Filhol
    
    Parameters
    ----------
    **M** :      a 2D array to detrend
    
    Returns
    -------
    **D** :      detrended array (2D)
    
    **coeff** :  coefficients of the linear trend such as D = M - (a*X + b*Y + c). X and Y are 2D arrays each containing the X and Y coordinates for each grid cell.
    '''     
    nx = np.size(M,0)
    ny = np.size(M,1)
    [Y, X] = np.mgrid[0:nx,0:ny]

    if ny*nx>10000:
        inds=np.zeros((nx, ny),dtype='int8')
        inds=inds.flatten()
        inds[0:9000]=1
        np.random.shuffle(inds)
        A=X.flatten()        
        A=A[inds==1]        
        B=Y.flatten()        
        B=B[inds==1]
        C=M.flatten()        
        C=C[inds==1]
        ABC=np.array([A, B, C])
        points = np.reshape(ABC, (np.shape(ABC)[0], -1))
        assert points.shape[0] < points.shape[1]
        ctr = points.mean(axis=1)
        x = points - ctr[:,None]
        Mat = np.dot(x, x.T)
        norm = np.linalg.linalg.svd(Mat)[0][:,-1]
        
    else:
        XYM = np.array([X.flatten(), Y.flatten(), M.flatten()])
        points = np.reshape(XYM, (np.shape(XYM)[0], -1))
        assert points.shape[0] < points.shape[1]
        ctr = points.mean(axis=1)
        x = points - ctr[:,None]
        Mat = np.dot(x, x.T)
        norm=np.linalg.linalg.svd(Mat)[0][:,-1]
    # for a plane with equation z = ax + by + c
    # at each (x,y) point, subtract the value of the fitted plane from dem
    d=-norm[0]*ctr[0]-norm[1]*ctr[1]-norm[2]*ctr[2]
    PlaneFit= -(norm[0]*X+ norm[1]*Y + d)/norm[2]
    D=M-PlaneFit
    coeff=np.array([norm,d])
    print ('Matrix detrended!')
    return D, coeff

#======================= Test Zone ================

