# -*- coding: utf-8 -*-
"""
Created on Thu May 15 19:32:08 2014

@author: simonfilhol
"""
import numpy as np

def myformula(X,Y,M,formula):
    return eval(formula)

def surfaceMaker(MatSize,Coords,fun):
    '''
        **D, dx, dy =  surfaceMaker(MatSize, Coords, fun="f(X,Y,M)")**\n
    create a 2D array based on input dimension and arithmetic combinations\n
    Copyright S. Filhol
    
    Parameters
    ----------
    **Matsize** :      a tuple containing the size of the matrix (Xsize, Ysize)\n
    **Coords** : a tuple containing the range of coordinates along X and Y-axis (Xstart, Xend, Ystart, Yend)

    Returns
    -------
    **D** :      detrended array (2D) \n
    **transform** : transform paramters use for saving matrix as rater, tuple of size 6
    (topleft X, dx W-E, rotation 0 if Norh up, topleft Y, dy N-S, rotation 0 )
    
    '''
    if np.size(MatSize)!=2:
        print "MatSize argument requires 2 integers as input such as (Xsize, Ysize)"
        return
    if np.size(Coords)!=4:
        print "Coords argument requires 4 real number such as (Xstart, Xend, Ystart, Yend)"
        return
    Xsize=MatSize[0]
    Ysize=MatSize[1]
    Xstart=Coords[0]
    Xend= Coords[1]
    Ystart=Coords[2]
    Yend=Coords[3]
    M=np.ones((Ysize,Xsize),float)
    x=np.linspace(Xstart,Xend,num=Xsize)
    y=np.linspace(Ystart,Yend,num=Ysize)
    X=M*x
    Y=y*M.T
    Y=Y.T
    D=myformula(X,Y,M,fun)
    dx=x[1]-x[0]
    dy=y[1]-y[0]
    transform=(Xstart, dx, 0, Ystart, dy, 0) #(topleft X, dx W-E, rotation 0 if Norh up, topleft Y, )
    return D, transform
   



#===================== Test Zone ================
'''    
Xsize=6
Ysize=5
Xstart=0
Xend= 100
Ystart=0
Yend=100
myOnes=np.ones((Ysize,Xsize),float)
x=np.linspace(Xstart,Xend,num=Xsize)
y=np.linspace(Ystart,Yend,num=Ysize)
Xmat=myOnes*x
Ymat=y*myOnes.T
Ymat=Ymat.T
M= - 1.5*myOnes*Ymat+ myOnes*x

D = detrend(sinmat)
'''
