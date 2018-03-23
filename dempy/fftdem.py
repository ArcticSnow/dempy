# -*- coding: utf-8 -*-
import numpy as np
import pyfftw
from dempy import linearDetrend


#================================================== PROCESSING FUNCTIONS
def hann2d(M):
    '''
    Function to ... blabla \n
    **H, Wss = hann2d(M)** \n
    dependencies : numpy    \n
    Parameters
    ------------
    **M** : 2D array
    Returns
    -------
    **H** : 2D array \n
    **Ws** : coefficient
    '''
    nx = np.size(M, 1)
    ny = np.size(M, 0)
    a = (nx + 1) / 2
    b = (ny + 1) / 2
    [X, Y] = np.meshgrid(range(0,nx),range(0,ny))
    
    theta = ((X==a))*(np.pi/2) + ((X!=a))*np.arctan2((Y-b),(X-a))
    r = np.sqrt((Y-b)**2 + (X-a)**2)
    rprime = np.sqrt((a**2)*(b**2)*(b**2*(np.cos(theta))**2 + a**2*(np.sin(theta))**2)**(-1))
    hanncoeff = (r < rprime)*(0.5*(1 + np.cos(np.pi*r/rprime)))
    H = M*hanncoeff
    Wss = sum(sum(hanncoeff**2))
    return H, Wss


def fftdem(dem, dx, dy=None, pad=None, window=None):
    '''
    **M, Pmat, fmat, Pvec, fvec = fftdem(dem, dx, dy, pad, window)** \n
    Compute the Fourier transform of a Digital Elevation Model (DEM) (2D   array) \n
    dependency: hann2D()
    
    Parameters
    ----------
    **dem** : 2D array of a DEM \n
    **dx, dy** : spatial resolution of the DEM. if dy = None then dy = dx \n
    **pad** : boolean, pading the dem with zeros \n
    **window** : 
    
    Returns
    -------
    **M** : 2D array of the padded DEM \n
    **Pmat** : matrix of Fourier coefficients\n
    **fmat** : matrix of frequencies\n
    **Pvec** : vector contianing all Fourier coefficients sorted based on frequencies (fvec)\n
    **fvec** : vector containing all frequencies sorted in ascending order
    '''
    if pad==None:
        pad=0
    if window==None:
        window=0
    if dy==None:
        dy=dx
        
    nx = np.size(dem,1)
    ny = np.size(dem,0)
    if np.logical_not(pad) and (np.remainder(nx, 2) or np.remainder(ny, 2)):
        print('If either dimension of the input matrix is odd, it is recommended to pad with zeros.')
    if window:
        dem, Wss = hann2d(dem)
    else:
        Wss = sum(sum(np.ones([ny, nx])))
    
    if pad:
        #calculate the power of 2 to pad with zeros 
        Lx = int(2**(np.ceil(np.log(np.max([nx,ny]))/np.log(2))))      
        Ly = Lx
    else: # no zero padding
        Lx = nx
        Ly = ny
    
    # calculate the frequency increments: frequency goes from zero (DC) to
    # 1/(2*dx) (Nyquist in x direction) in Lx/2 increments analogous for y.
    dfx = 1/(dx*Lx)
    dfy = 1/(dy*Ly)
    # Do a 2D FFT, padding with zeros.
    # After fftshift, the power spectra have zero frequency (DC) at the grid
    # point (Ly/2 + 1, Lx/2 + 1), e.g., (257,257) for a 512x512 
    # result. Since the zero point is offset from the center of the grid, the 
    # min and max frequencies on the two axes are different (specifically, the 
    # point at (257,512) is one bin below the Nyquist frequency, whereas the 
    # point at (257,1) corresponds to the Nyquist frequency).
    xc = int(Lx/2)
    yc = int(Ly/2) # matrix indices of zero frequency
    
    dem2 = pyfftw.n_byte_align_empty((Lx,Ly), 16, 'complex128')
    dem2 = pyfftw.interfaces.numpy_fft.fftshift(pyfftw.interfaces.numpy_fft.fft2(dem,(Ly,Lx)))
    #dem = np.fft.fftshift(np.fft.fft2(dem,(Ly,Lx)))
    
    if window == 1 or pad == 1 :
        M = pyfftw.interfaces.numpy_fft.ifft2(pyfftw.interfaces.numpy_fft.ifftshift(dem2))
        M = np.real(M)
    else:
        M = dem
    #M = np.real(np.fft.ifft2(np.fft.ifftshift(dem)))
    dem2[yc, xc] = 0 # Although the mean of the detrended
                           # matrix is zero, the mean of the
                           # windowed matrix may not be,
                           # so here we zero out the DC
                           # component (the data mean)
                           
    # Calculate the DFT periodogram, which has units of amplitude^2, or m^2 for
    # topography. Dividing by Wss*twopower^2 corrects for the reduction of
    # amplitude by the windowing function (see Numerical Recipes 13.4)
    dem2 = np.real(dem2 * np.conj(dem2) / (Lx * Ly * Wss))
    
    # assign the power spectrum to the output argument
    Pmat = dem2
    
    # Create a matrix of radial frequencies
    
    [cols, rows] = np.meshgrid(range(0,Lx),range(0,Ly)) # matrices of column and row indices
    fmat = np.sqrt((dfy*(rows-yc))**2 + (dfx*(cols-xc))**2) # frequency matrix
    
    
    # Create sorted, non-redundant vectors of frequency and power 
    dem2 = dem2[:,range(0,xc+1)]
    # fvec = dfreq*sqrt((rows(:,1:xc)-yc)**2 + (cols(:,1:xc)-xc)**2)
    fvec = fmat[:,range(0,xc+1)]
    
    fvec[range((yc+1),Ly),xc] = -1 # This half-column is redundant. Set the 
                         # frequencies to negative values so they 
                         # will be clipped out below
    myfvec=np.ones((fvec.flatten().size,2),dtype=float)    
    myfvec[:,0]=fvec.flatten()
    myfvec[:,1]=dem2.flatten()
    fvec=myfvec
    myfvec=None
    fvec=fvec[fvec[:,0].argsort()]
    fvec = fvec[fvec[:,0]>0,:] # Take only positive frequencies. This gets rid 
                            # of DC (zero freq) as well as the redundant 
                            # frequencies we set to -1 above
                            
    # Separate into power and frequency vectors and assign to output arguments
    Pvec = 2*fvec[:,1] # the factor of 2 corrects for the fact that we have
                    # taken only half of the 2D spectrum. sum(Pvec) should
                    # now equal sum(Pmat(:)).
    fvec = fvec[:,0]
    
    return M, Pmat, fmat, Pvec, fvec


def gaussian(freqmat, mu, sigma):
    '''
    Function used in filtdem
    Translation in Python by S. Filhol
    dependecies : numpy
    '''
    G=np.exp(-(freqmat-mu)**2/(2*sigma**2))
    G=G/np.max(G.flatten())
    return G


def specfilt2d(fmat,f,filtype):
    '''
    Function used in filtdem
    
    dependecies : numpy, gaussian()
    '''
    if filtype=="lowpass":
        if np.size(f)!=2:
            print("For lowpass filter, -f- requires 2 frequency indications such as: (flow,fhigh)=(.5,1)")
            return
        flo = f[0]
        fhi = f[1]
        sigma=np.abs(fhi-flo)/3
        F=gaussian(fmat,flo,sigma)
        F[fmat<flo]=1
        
    elif filtype=="highpass":
        print('\n ... Need to check code for Highpass because does not produce good product ... \n')
        return
        if np.size(f)!=2:
            print("For highpass filter, -f- requires 2 frequency indications such as: (flow,fhigh)=(.5,1)")
            return
        flo = f[0]
        fhi = f[1]
        sigma=np.abs(fhi-flo)/3
        F=gaussian(fmat,fhi,sigma)
        F[fmat>=fhi]=1
        
    else:
        if filtype != "bandpass":
            print("filtype arguments must either be -lowpass-, -bandpass-, or -highpass- ")
            return
        if np.size(f)!=4:
            print("For bandpass filter, -f- requires 4 frequency indications such as: (flow1, flow2,fhigh1, fhigh2)=(.5,1,2,2.5)")
            return
        flo1 = f[0]
        flo2 = f[1]
        fhi1 = f[2]
        fhi2 = f[3]
        sigmalo = np.abs(flo2-flo1)/3
        sigmahi = np.abs(fhi2-fhi1)/3
        Flo=gaussian(fmat,flo2,sigmalo)
        Fhi=gaussian(fmat,fhi1,sigmahi)
        # Conditionale statement below will break... need to be adapted to python syntax!!!!!!!
        F = Flo*((fmat<=flo2)*fmat) + Fhi*((fmat>=fhi1)*fmat) + 1*((fmat>flo2 & fmat<fhi1)*fmat)
    return F
        



def filtdem(M,dx,dy,filtype, f):
    '''
    **M = filtdem(Pmat,Fmat,filtype, f)** \n
    filter a fourier transformed DEM into a topo 2D array \n
    dependencies: numpy, specfilt2d(), gaussian()
    
    Parameters
    ----------
    **Pmat** : 2D array, Fourier transform resutling from fftdem(). Contains the Fourier coefficients \n
    **Fmat** : 2D array, containing the corresponding frequencies \n
    **filtype** : string, 3 possibilities: "lowpass", "highpass", and "bandpass" \n
    **f** : tuple of frequencies for the filter. 
            - If filtype = "lowpass", f = (flow, fhigh)
            - If filtype = "highpass", f = (flow, fhigh)
            - If filtype = "bandpass", f = (flow1, flow2, fhigh1, fhigh2) with flow1 < flow2 < fhigh1 < fhigh2
    
    Returns
    -------
    **M** : 2D array, inversed transformed of Pmat
    '''
    nx = np.size(M,1)
    ny = np.size(M,0)
    pad=0
    window=0
    _,Pmat, fmat, _, _ = fftdem(M, dx, dy, pad, window)    
    M = linearDetrend.detrend(M)[0]    
    M_filt=pyfftw.n_byte_align_empty((nx,ny), 16, 'complex128')
    M_filt= pyfftw.interfaces.numpy_fft.fftshift(pyfftw.interfaces.numpy_fft.fft2(M))
    F = specfilt2d(fmat, f, filtype)
    # take the inverse FFT of the filtered spectrum
    M = np.real(pyfftw.interfaces.numpy_fft.ifft2(pyfftw.interfaces.numpy_fft.ifftshift(M_filt*F)))
    return M
    
    
    
    
#====================== SCRIPTING ======================
'''
myRaster=openRaster("/Users/simonfilhol/Desktop/Dune_project/s13w10_processing/height_grid_raster3_shade.tif")

mydata, dx, dy, transform = raster2array(myRaster)
dem, Pmat, freqmat, Pvec, fvec = fftdem(mydata,dx,dx,1,1)
plotPeriodogram(fvec,Pvec)

"""
Created on Mon May 12 15:16:27 2014

@author: simonfilhol, May 2014
FFT processing based on a matlab implementation written by Perron et al.
"""

# fftdem function. Function to estimate fourier transform of a DEM.
# work to accomplish:
#   - correct Hann2 funtion: problem of array size in the line H=m*hanncoef


sinmat=np.ones((100,100),float)
x=range(0,100)

sinmat=sinmat[:,range(0,100)]*x
sinmat=sin(sinmat)+sin(sinmat*3)+cos(sinmat*5)
# Pmat1, freqmat1, Pvec1, fvec1 = fftdem(sinmat,1,1,0,0)
# plotPeriodogram(fvec1,Pvec1)
# Pmat2, freqmat2, Pvec2, fvec2 = fftdem(sinmat,1,1,1,0)
# plotPeriodogram(fvec2,Pvec2)
Dem3, Pmat3, freqmat3, Pvec3, fvec3 = fftdem(sinmat,1,1,0,0)  # option that provides best results for this artificial surface
plotPeriodogram(fvec3,Pvec3)
# Pmat4, freqmat4, Pvec4, fvec4 = fftdem(sinmat,1,1,0,1)
# plotPeriodogram(fvec4,Pvec4)

sinmatFiltered=filtdem(Pmat3,freqmat3,(.38,.4),"lowpass")
_, _, _, Pvec4, fvec4 = fftdem(sinmatFiltered,1,1,0,0) 
plotPeriodogram(fvec4,Pvec4)
'''
