"""
Tools to copmute terrain morphometric variables
S. Filhol March 2026


"""
import numpy as np
import cv2
from numba import njit, prange, config, threading_layer

config.THREADING_LAYER = 'threadsafe'


def compute_aspect(dem, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    aspect = np.zeros_like(dem, dtype=np.float32)

    @njit(parallel=True)
    def _compute_aspect(arr, padded_dem, radius):
        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                dz_dx = (padded_dem[i, j+1] - padded_dem[i, j-1]) / 2.0
                dz_dy = (padded_dem[i-1, j] - padded_dem[i+1, j]) / 2.0  # Note the change in dz_dy
                arr[i-radius, j-radius] = np.arctan2(dz_dx, dz_dy)  # Swapped arguments
        return arr

    aspect = _compute_aspect(aspect, padded_dem, radius)

    # Convert the aspect to have 0 at North and rotate clockwise
    aspect = (np.pi/2 - aspect) % (2 * np.pi)

    return aspect



def compute_twi(dem, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    tpi = np.zeros_like(dem, dtype=np.float32)

    @njit(parallel=True)
    def _compute_tpi(arr, padded_dem, radius):

        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
                mean_neighbor = np.mean(window)
                arr[i-radius, j-radius] = padded_dem[i, j] - mean_neighbor
        return arr

    tpi = _compute_tpi(tpi, padded_dem, radius)
    return tpi


@njit()
def meshgrid_3D(x, y, z):
    xx = np.empty(shape=(x.size, y.size, z.size), dtype=x.dtype)
    yy = np.empty(shape=(x.size, y.size, z.size), dtype=y.dtype)
    zz = np.empty(shape=(x.size, y.size, z.size), dtype=z.dtype)
    for i in range(z.size):
        for j in range(y.size):
            for k in range(x.size):
                xx[i,j,k] = x[k]  # change to x[k] if indexing xy
                yy[i,j,k] = y[j]  # change to y[j] if indexing xy
                zz[i,j,k] = z[i]  # change to z[i] if indexing xy
    return zz, yy, xx

@njit()
def meshgrid_2D(x, y):
    xx = np.empty(shape=(x.size, y.size), dtype=x.dtype)
    yy = np.empty(shape=(x.size, y.size), dtype=y.dtype)
    for j in range(y.size):
        for k in range(x.size):
            xx[j,k] = x[k]  # change to x[k] if indexing xy
            yy[j,k] = y[j]  # change to y[j] if indexing xy
    return yy, xx



def compute_curvature_multiscale(dem, window_size=3, pixel_size=1.0):
    """
    Compute metric, downslope profile and plan curvature from a DEM
    using a local quadratic fit over a moving window.

    Parameters
    ----------
    dem : 2D ndarray
        Digital elevation model, elevation in meters.
    window_size : int, optional
        Odd size of moving window (3, 5, 7, ...). Default is 3.
    pixel_size : float, optional
        DEM cell size in meters (assumes square pixels). Default is 1.0.

    Returns
    -------
    profile_curv : 2D ndarray (float32)
        Downslope profile curvature (1/m).
    plan_curv : 2D ndarray (float32)
        Plan curvature (1/m).
    """
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")

    radius = window_size // 2
    radius_float = float(radius)

    dem64 = dem.astype(np.float64, copy=False)
    profile_curv = np.zeros_like(dem64, dtype=np.float64)
    plan_curv = np.zeros_like(dem64, dtype=np.float64)

    # Pad DEM
    padded_dem = np.pad(dem64, radius, mode="edge")

    # Precompute coordinates in *pixel units* (indices)
    # We will convert to metric inside the kernel using pixel_size.
    x_coords, y_coords = np.meshgrid(
        np.arange(-radius_float, radius_float + 1, dtype=np.float64),
        np.arange(-radius_float, radius_float + 1, dtype=np.float64),
        indexing="xy"
    )
    x_flat = x_coords.ravel()
    y_flat = y_coords.ravel()
    len_x = x_flat.shape[0]

    @njit(parallel=True)
    def _compute_curv(profile_arr, plan_arr, padded_dem,
                      radius, x_flat, y_flat, len_x, pixel_size):

        rows = profile_arr.shape[0]
        cols = profile_arr.shape[1]

        for i in prange(radius, rows + radius):
            for j in range(radius, cols + radius):
                # Neighborhood window
                window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
                z = window.ravel()

                # Design matrix in *pixel* coordinates:
                # z = a i^2 + b i j + c j^2 + d i + e j + f
                A = np.empty((len_x, 6), dtype=np.float64)
                A[:, 0] = x_flat * x_flat    # i^2
                A[:, 1] = x_flat * y_flat    # i j
                A[:, 2] = y_flat * y_flat    # j^2
                A[:, 3] = x_flat             # i
                A[:, 4] = y_flat             # j
                A[:, 5] = 1.0                # constant

                # Least squares fit
                coeffs = np.linalg.lstsq(A, z, rcond=-1.0)[0]
                a = coeffs[0]
                b = coeffs[1]
                c = coeffs[2]
                d = coeffs[3]
                e = coeffs[4]
                # f = coeffs[5]

                # ------------------------------------------------------
                # Convert derivatives from pixel to metric units
                # ------------------------------------------------------
                # Coordinates:
                #   X = i * pixel_size, Y = j * pixel_size
                # Chain rule:
                #   p_phys = ∂z/∂X = (1/pixel_size) * ∂z/∂i
                #   q_phys = ∂z/∂Y = (1/pixel_size) * ∂z/∂j
                #   r_phys = ∂²z/∂X² = (1/pixel_size²) * ∂²z/∂i², etc.

                h = pixel_size

                # First derivatives at center (i=0,j=0), in pixel space:
                # ∂z/∂i = d, ∂z/∂j = e
                # In metric:
                p = d / h          # z_X (dimensionless slope)
                q = e / h          # z_Y

                # Second derivatives in pixel space:
                # ∂²z/∂i² = 2a, ∂²z/∂i∂j = b, ∂²z/∂j² = 2c
                # In metric (1/m):
                r = 2.0 * a / (h * h)   # z_XX
                s = b      / (h * h)    # z_XY
                t = 2.0 * c / (h * h)   # z_YY

                # ------------------------------------------------------
                # Metric geomorphometric curvatures (Evans / Z&T)
                # ------------------------------------------------------
                # p = z_X, q = z_Y, r = z_XX, s = z_XY, t = z_YY
                p2 = p * p
                q2 = q * q
                G = p2 + q2

                if G == 0.0:
                    profile = 0.0
                    plan = 0.0
                else:
                    sqrtG  = np.sqrt(G)
                    sqrt1G = np.sqrt(1.0 + p2 + q2)

                    # Upslope profile curvature (1/m)
                    num_prof = p2 * r + 2.0 * p * q * s + q2 * t
                    k_prof_up = num_prof / (sqrtG * (sqrt1G**3))

                    # Plan curvature (1/m)
                    num_plan = q2 * r - 2.0 * p * q * s + p2 * t
                    k_plan = num_plan / (G * sqrt1G)

                    # Convert to *downslope* profile curvature
                    profile = -k_prof_up
                    plan = k_plan

                profile_arr[i-radius, j-radius] = profile
                plan_arr[i-radius, j-radius] = plan

        return plan_arr, profile_arr

    plan_curv, profile_curv = _compute_curv(
        profile_curv, plan_curv,
        padded_dem,
        radius,
        x_flat, y_flat, len_x,
        float(pixel_size),
    )

    return profile_curv, plan_curv





def compute_tpi(dem, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    tpi = np.zeros_like(dem, dtype=np.float32)
    
    @njit(parallel=True)
    def _compute_tpi(arr, padded_dem, radius):
        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
                mean_neighbor = np.mean(window)
                arr[i-radius, j-radius] = padded_dem[i, j] - mean_neighbor
        return arr
    
    tpi = _compute_tpi(tpi, padded_dem, radius)
    return tpi


def compute_tri(dem, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    tri = np.zeros_like(dem, dtype=np.float32)

    @njit(parallel=True)
    def _compute_tri(arr, padded_dem, radius):
        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                center = padded_dem[i, j]
                window = padded_dem[i-radius:i+radius+1, j-radius:j+radius+1]
                arr[i-radius, j-radius] = np.sqrt(np.sum((window - center)**2) / (window.size - 1))
        return arr

    return _compute_tri(tri, padded_dem, radius)


def compute_slope(dem, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2  # Convert window_size to radius
    padded_dem = np.pad(dem, radius, mode='edge')
    rows, cols = dem.shape
    slope = np.zeros_like(dem, dtype=np.float32)

    @njit(parallel=True)
    def _compute_slope(arr, padded_dem, radius):
        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                dz_dx = (padded_dem[i, j+1] - padded_dem[i, j-1]) / 2.0
                dz_dy = (padded_dem[i+1, j] - padded_dem[i-1, j]) / 2.0
                arr[i-radius, j-radius] = np.arctan(np.sqrt(dz_dx**2 + dz_dy**2))

        return arr
    
    slope = _compute_slope(slope, padded_dem, radius)
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


def compute_spi_classical(flow_acc, slope):
    return flow_acc * np.tan(slope)


def compute_spi(flow_acc, slope, window_size=3):
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd (e.g., 3, 5, 7).")
    radius = window_size // 2  # Convert window_size to radius
    padded_flow_acc = np.pad(flow_acc, radius, mode='edge')
    padded_slope = np.pad(slope, radius, mode='edge')
    rows, cols = flow_acc.shape
    spi = np.zeros_like(flow_acc, dtype=np.float32)


    @njit(parallel=True)
    def _compute_spi(arr, padded_slope, padded_flow_acc, radius):
        for i in prange(radius, rows + radius):
            for j in prange(radius, cols + radius):
                if padded_slope[i, j] > 0 and padded_flow_acc[i, j] > 0:
                    arr[i-radius, j-radius] = padded_flow_acc[i, j] * np.tan(padded_slope[i, j])
        return arr
    return _compute_spi(spi, padded_slope, padded_flow_acc, radius)


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
