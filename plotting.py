# -*- coding: utf-8 -*-
"""
Created on Thu May 15 23:06:08 2014

@author: simonfilhol
"""
from matplotlib import pyplot as plt
def thinningCheckup(fvec,DispOption=None):
    '''
    Function part of plotPeriodogram()
    Applies default thinning if the input vectors are too large
    '''
    if DispOption!=None:
        print "There are " + str(fvec.size) + " data points \n"
        print "You should use a thinning of " + str(int(fvec.size/50000))
    thinning = int(fvec.size/10000)
    return thinning

def plotPeriodogram(fvec, pvec,  axes=None, thinning=None):
    '''
    **plotPeriodogram(fvec, pvec,  axes=None, thinning=None)**
    Function that plot a periodogram
    
    Parameters
    ----------
    **fvec** : vector containing  the frequencies resulting from fftdem() \n
    **pvec** : vector containing  the power values resulting from fftdem() \n
    **axes** : string indicating what type of axes to use. Possible options are: \n\t
        - "loglog" \n\t
        - "semilogx"\n\t
        - "semilogy" \n\t
        - None (default option)\n
    **thinning** : parameter to thin the data o plot as vectors can be ver large. It will plot only the number of dots indicated by thinning '''
    # Wvec=1/fvec
    if thinning == None:
        thinning = thinningCheckup(fvec)
    plt.figure()
    if axes == "loglog":
        plt.loglog(fvec[range(0,fvec.size,thinning)],pvec[range(0,pvec.size,thinning)])
    elif axes == "semilogx":
        plt.semilogx(fvec[range(0,fvec.size,thinning)],pvec[range(0,pvec.size,thinning)])
    elif axes == "semilogy":
        plt.semilogy(fvec[range(0,fvec.size,thinning)],pvec[range(0,pvec.size,thinning)])
    else:
        plt.plot(fvec[range(0,fvec.size,thinning)],pvec[range(0,pvec.size,thinning)])
    plt.title("Periodogram")
    plt.ylabel("DFT mean square amplitude")
    plt.xlabel("Frequency (1/m)")
    plt.show()
    
def plotMap(Map,Zmin,Zmax,Cmap='gray',Title=None):
    '''
    Function to plot map image
    **plotMap(Map,Zmin,Zmax,Cmap='gray',Title=None)**
    
    Parameters
    ----------
    **Map** : 
    '''
    plt.figure()
    plt.imshow(Map,vmin=Zmin, vmax=Zmax, cmap=Cmap,origin="lower")
    plt.colorbar()
    if Title != None:
        plt.title(Title)
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    