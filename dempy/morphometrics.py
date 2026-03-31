"""
Tools to copmute terrain morphometric variables
S. Filhol March 2026


"""
import numpy as np
import cv2
from numba import njit, prange


@njit(parallel=True)
def compute_aspect(dem, radius=1):
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    aspect = np.zeros_like(dem, dtype=np.float32)

    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            dz_dx = (padded_dem[i, j+1] - padded_dem[i, j-1]) / 2.0
            dz_dy = (padded_dem[i+1, j] - padded_dem[i-1, j]) / 2.0
            aspect[i-radius, j-radius] = np.degrees(np.arctan2(-dz_dy, dz_dx)) % 360

    return aspect


@njit(parallel=True)
def compute_twi(dem, radius=3):
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    tpi = np.zeros_like(dem, dtype=np.float32)

    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
            mean_neighbor = np.mean(window)
            tpi[i-radius, j-radius] = padded_dem[i, j] - mean_neighbor

    return tpi

@njit(parallel=True)
def compute_curvature_multiscale(dem, radius=1):
    rows, cols = dem.shape
    profile_curv = np.zeros_like(dem, dtype=np.float32)
    plan_curv = np.zeros_like(dem, dtype=np.float32)

    # Pad the DEM to handle edges
    padded_dem = np.pad(dem, radius, mode='edge')

    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            # Extract the neighborhood
            window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
            n = window.shape[0]

            # Fit a quadratic surface: z = ax^2 + bxy + cy^2 + dx + ey + f
            # Solve for coefficients using least squares
            x, y = np.meshgrid(np.arange(-radius, radius+1), np.arange(-radius, radius+1))
            x = x.flatten()
            y = y.flatten()
            z = window.flatten()

            # Design matrix for quadratic surface
            A = np.vstack([
                x**2, x*y, y**2, x, y, np.ones_like(x)
            ]).T

            # Solve for coefficients (a, b, c, d, e, f)
            coeffs, _, _, _ = np.linalg.lstsq(A, z, rcond=None)

            # Second derivatives
            a, b, c, d, e, f = coeffs
            d2z_dx2 = 2 * a
            d2z_dy2 = 2 * c
            d2z_dxdy = b

            # First derivatives (for slope and aspect)
            dz_dx = d
            dz_dy = e
            slope = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))
            aspect_rad = np.arctan2(-dz_dy, dz_dx)

            # Profile and plan curvature
            profile_curv[i-radius, j-radius] = (
                d2z_dx2 * np.cos(aspect_rad)**2 +
                2 * d2z_dxdy * np.sin(aspect_rad) * np.cos(aspect_rad) +
                d2z_dy2 * np.sin(aspect_rad)**2
            )
            plan_curv[i-radius, j-radius] = (
                d2z_dx2 * np.sin(aspect_rad)**2 -
                2 * d2z_dxdy * np.sin(aspect_rad) * np.cos(aspect_rad) +
                d2z_dy2 * np.cos(aspect_rad)**2
            )

    return profile_curv, plan_curv

@njit(parallel=True)
def compute_tri_classical(dem):
    rows, cols = dem.shape
    tri = np.zeros_like(dem, dtype=np.float32)
    for i in prange(1, rows - 1):
        for j in prange(1, cols - 1):
            center = dem[i, j]
            neighbors = [
                dem[i-1, j-1], dem[i-1, j], dem[i-1, j+1],
                dem[i, j-1],               dem[i, j+1],
                dem[i+1, j-1], dem[i+1, j], dem[i+1, j+1]
            ]
            tri[i, j] = np.sqrt(np.sum((neighbors - center)**2))
    return tri


@njit(parallel=True)
def calculate_tri_radius(dem, radius=1):
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    tri = np.zeros_like(dem, dtype=np.float32)

    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            center = padded_dem[i, j]
            window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
            tri[i-radius, j-radius] = np.sqrt(np.sum((window - center)**2) / (window.size - 1))

    return tri



@njit(parallel=True)
def compute_slope(dem, radius=1):
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    slope = np.zeros_like(dem, dtype=np.float32)

    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            dz_dx = (padded_dem[i, j+1] - padded_dem[i, j-1]) / 2.0
            dz_dy = (padded_dem[i+1, j] - padded_dem[i-1, j]) / 2.0
            slope[i-radius, j-radius] = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))

    return slope


@njit(parallel=True)
def compute_flow_accumulation(dem):
    rows, cols = dem.shape
    flow_acc = np.zeros_like(dem, dtype=np.float32)
    # Simple D8 flow routing (for demonstration; consider using a library for large DEMs)
    for i in prange(1, rows - 1):
        for j in prange(1, cols - 1):
            center = dem[i, j]
            neighbors = [
                (dem[i-1, j-1], (i-1, j-1)), (dem[i-1, j], (i-1, j)), (dem[i-1, j+1], (i-1, j+1)),
                (dem[i, j-1], (i, j-1)),                     (dem[i, j+1], (i, j+1)),
                (dem[i+1, j-1], (i+1, j-1)), (dem[i+1, j], (i+1, j)), (dem[i+1, j+1], (i+1, j+1))
            ]
            # Find the steepest downslope neighbor
            min_val = center
            min_idx = (i, j)
            for val, (ni, nj) in neighbors:
                if val < min_val:
                    min_val = val
                    min_idx = (ni, nj)
            # Accumulate flow
            if min_idx != (i, j):
                flow_acc[min_idx] += 1
    return flow_acc



@njit
def compute_spi(flow_acc, slope):
    return flow_acc * np.tan(slope)


@njit(parallel=True)
def calculate_spi(flow_acc: np.ndarray, slope: np.ndarray, radius: int = 1) -> np.ndarray:
    padded_flow_acc = np.pad(flow_acc, radius, mode='edge')
    padded_slope = np.pad(slope, radius, mode='edge')
    rows, cols = flow_acc.shape
    spi = np.zeros_like(flow_acc, dtype=np.float32)
    for i in prange(radius, rows + radius):
        for j in prange(radius, cols + radius):
            if padded_slope[i, j] > 0 and padded_flow_acc[i, j] > 0:
                spi[i-radius, j-radius] = padded_flow_acc[i, j] * np.tan(padded_slope[i, j])
    return spi

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


def convolve_array(arr, kernel):
    '''
    function to convolve a 2D array with a kernel
    '''
    res = cv2.filter2D(arr, -1, kernel)
    return res
