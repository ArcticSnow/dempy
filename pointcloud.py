from __future__ import division
from laspy.file import File
import numpy as np
import pandas as pd
import time, math

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0))
        return ret
    return wrap

@timing
def loadLAS2XYZ(filepath):
    '''
    Function to load in console the pointcloud of a LAS file
    :param filepath: filepath of the LAS file
    :return: xyz array containing coordinate of the points
    '''
    print('Start loading...')
    inFile = File(filepath, mode='r')
    coords = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
    print('Data loaded')
    return coords


@timing
def loadLAS2XYZ(filepath):
    '''
    Function to load in console the pointcloud of a LAS file with points attributes
    :param filepath: filepath of the LAS file
    :return: xyz array containing coordinate of the points
    '''
    print('Start loading...')
    inFile = File(filepath, mode='r')
    coords = np.vstack((inFile.x, inFile.y, inFile.z, inFile.amplitude, inFile.Intensity, inFile.reflectance, inFile.num_returns)).transpose()
    print('Data loaded')
    return coords


def xyz2binarray(xyz, xstart, xend, ystart, yend, nx=1000, ny=1000, method='min'):
    '''
    Function to extract projected grid on the XY-plane of point cloud statistics

    :param xyz: a 3 column vector containing the point location in cartesian coordinate system
    :param xstart: x-minimum of the grid
    :param xend: x-maximum of the grid
    :param ystart: y-minimm of the grid
    :param yend: y-maximum of the grid
    :param nx: number of grid cell in the x directions
    :param ny: number of grid cell in the y directions
    :param method: statistics to extract from each gridcell
    :return: returns a 2D array, xmin, and ymax

    TO IMPLEMENT:
        - being able to choose to input dx dy instead of nx ny
    '''
    binned, bins_x, bins_y, bin_xmin, bin_ymin = binData2D(xyz, xstart, xend, ystart, yend, nx, ny)

    if method == 'min':
        ret = binned.Z.min().unstack().T  # .iloc[::-1]
    elif method == 'max':
        ret = binned.Z.max().unstack().T  # .iloc[::-1]
    elif method == 'mean':
        ret = binned.Z.mean().unstack().T  # .iloc[::-1]
    elif method == 'median':
        ret = binned.Z.median().unstack().T  # .iloc[::-1]
    elif method == 'count':
        ret = binned.Z.count().unstack().T  # .iloc[::-1]

    xmin = bins_x[ret.columns.min().astype(int)]
    ymax = bins_y[ret.index.get_values().max().astype(int)]

    newIndy = np.arange(ret.index.get_values().min(), ret.index.get_values().max() + 1)
    newIndx = np.arange(ret.columns.min(), ret.columns.max() + 1)
    a = ret.reindex(newIndy, newIndx)
    mat = np.zeros((ny, nx)) * np.nan
    mat[bin_ymin:bin_ymin + a.shape[0], bin_xmin:bin_xmin + a.shape[1]] = a

    return mat[::-1], xmin, ymax

def LAS2txt(filepath,newfile):
    '''
    Function to convert a pointcloud save in LAS format into a .txt format
    :param filepath: filepath of the LAS file
    :param newfile: name of the new file
    :return: save data into a text file
    '''
    inFile = File(filepath, mode='r')
    coords = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
    if newfile[-4] != '.txt':
        newfile = newfile + '.txt'
    np.savetxt(newfile,coords)
    print('File saved: ' + newfile)

def xyz_subsample(xyz, length_out):
    '''
    Function to subsample a 3 columm matrix.
    :param xyz: 3 column matrix
    :param length_out: number of sample to output
    :return: a 3 column matrix
    '''
    ind = np.random.randint(0,xyz.shape[0],length_out)
    xyz_new = xyz[ind,:]
    print('xyz subsampled!')
    return xyz_new

def xyz_stat(xyz):
    print('Shape of array: ' + str(xyz.shape))
    print('Min of xyz: ')
    print(np.min(xyz, axis=0))
    print('Max of xyz: ')
    print(np.max(xyz, axis=0))
    print('Mean of xyz: ')
    print(np.mean(xyz, axis=0))
    print('Extent')
    print(np.max(xyz, axis=0)-np.min(xyz, axis=0))

def trans(xyz,trans_vec):
    '''
    Function to translate an xyz 3 column matrix
    :param xyz: a 3 column matrix
    :param trans_vec: a translation vector of length 3
    :return: a 3 column matrix translated
    '''

    xyz[:,0] = xyz[:,0] - trans_vec[0]
    xyz[:,1] = xyz[:,1] - trans_vec[1]
    xyz[:,2] = xyz[:,2] - trans_vec[2]
    return xyz

def translate_coords(coords, xyz_trans = None ,ask = True):
    '''
    Function to translate a point cloud
    :param coords: an xyz array
    :param xyz_trans: vector of translation in [x,y,z]
    :param ask: if True (default) brings an interactive console for approving the translation
    :return: translated xyz array
    '''
    if xyz_trans is None:
        xyz_trans = [coords[:,0].min(), coords[:,1].min(), coords[:,2].min()]
    if ask is True:
        print('Default translation:')
        print(str(xyz_trans) + '\n')
        res = input('Do you want to translate? 0/1')
        if res is 0:
            print('No Translation applied')
            return None
        if res is 1:
            return trans(coords, xyz_trans)
    if ask is not True:
        return trans(coords, xyz_trans)

def truncate(xyz, Xextent, Yextent):
    '''
    Function to truncate a point cloud with a rectangular shape
    :param xyz: a 3 column matrix containing the points coordinate
    :param Xextent: a vector of Xmin and Xmax (e.g. [Xmin,Xmax])
    :param Yextent: a vector of Ymin and Ymax (e.g. [Ymin, Ymax])
    :return: a 3 colum matrix containing the points coordiante within the specified rectangle
    '''
    xcut = xyz[xyz[:,0]>=Xextent[0]]
    xcut1 = xcut[xcut[:,0]<Xextent[1]]
    ycut = xcut1[xcut1[:,1]>=Yextent[0]]
    ycut1 = ycut[ycut[:,1]<Yextent[1]]
    return ycut1

def cart2cyl(xyz, xy_axis=None):
    '''
        function to convert cartesian coordinates to cylindrical

    :param xyz: a 3-column matrix containing the points coordinates expressed in a cartesian system
    :param xy_axis: an array of x and y coordinate for the center of the new cylindrical coordinate
    :return: a 3 colum matrix with the point coordinates are expressed in a cylindrical coordinate system
    '''
    if xy_axis is not None:
        xyz[:,0] =  xyz[:,0] - xy_axis[0]
        xyz[:,1] =  xyz[:,1] - xy_axis[1]
    rho = np.sqrt(xyz[:,0]**2 + xyz[:,1]**2)
    phi = np.arctan2(xyz[:,1], xyz[:,0])
    rpz = np.vstack((rho,phi,xyz[:,2]))
    return rpz.transpose()

def cyl2cart(rpz):
    '''
    convert cylindrical coordinate to cartesian
    :param rpz: a 3-column matrix containing the points coordinates expressed in a cylindrical system
    :return: a 3-column matrix containing the points coordinates expressed in a cartesian system
    '''
    x = rpz[:,0] * np.cos(rpz[:,1])
    y = rpz[:,0] * np.sin(rpz[:,1])
    xyz = np.vstack((x,y,rpz[:,2]))
    return xyz.transpose()

def rotate_cloud(xyz, angle, center_coord=None):
    '''
    Function to rotate a point cloud
    :param xyz: n*3 array containing point cloud coordinates in a cartesian system
    :param angle: angle of rotation in degrees
    :param center_coord: tuple with xy coordiantes of the center of rotation. Default is None
    :return: the rotated xyz point cloud
    '''
    if center_coord is None:
        center_coord = [np.mean(xyz[:,0]),np.mean(xyz[:,1])]
    rpz = cart2cyl(xyz, xy_axis=center_coord)
    rpz[:,1] = rpz[:,1] + angle
    xyz = cyl2cart(rpz)
    return xyz

def get_slice(xyz, thick, dir=0, center_coord=None):
    '''
    Function to extract a slice of the point cloud xyz
    :param xyz: n*3 array containing point cloud coordinates in a cartesian system
    :param thick: thickness of the slice
    :param dir: direction of the slice in degrees (default is 0)
    :param center_coord: tuple with xy coordinates of the center of rotation. Default is None
    :return: return slice in xyz format.
    '''
    if center_coord is None:
        center_coord = [np.mean(xyz[:,0]),np.mean(xyz[:,1])]
        print(center_coord)
    if dir % 180 != 0:
        xyz = rotate_cloud(xyz, (dir*math.pi/180), center_coord= center_coord)
    myslice = xyz[xyz[:,0]>=-(thick/2)]
    myslice = myslice[myslice[:,0]<=(thick/2)]
    return myslice

def get_slice_df(df_xyz, thick, dir=0, center_coord=None):
    '''
    Function to extract a slice of points from a dataframe
    :param xyz: n*3 array containing point cloud coordinates in a cartesian system
    :param thick: thickness of the slice
    :param dir: direction of the slice in degrees (default is 0)
    :param center_coord: tuple with xy coordinates of the center of rotation. Default is None
    :return: return slice in xyz format.
    '''
    df = df_xyz.copy()
    df_xyz=None
    if center_coord is None:
        center_coord = [df['x'].mean(),df['y'].mean()]
        print(center_coord)
    if dir % 180 != 0:
        xyz = rotate_cloud(np.array(df[['x','y','z']]), (dir*math.pi/180), center_coord = center_coord)
        df[['x','y']] = xyz[:,[0,1]]
        myslice = df[df.x >= - (thick / 2)]
        myslice = myslice[df.x <= (thick/2)]
    else:
        myslice = df[df.x >= (center_coord[0] - thick / 2)]
        myslice = myslice[df.x <= (center_coord[0] + thick / 2)]
        myslice['x'] = myslice['x'] - center_coord[0]
        myslice['y'] = myslice['y'] - center_coord[1]
    print('Data Sliced')
    return myslice

def center_pc_coord_df(df_xyz, center_coord=None):
    if center_coord is None:
        center_coord = [(df_xyz['x'].max()-df_xyz['x'].min())/2 + df_xyz['x'].min(),
                        (df_xyz['y'].max()-df_xyz['y'].min())/2 +df_xyz['y'].min()]
        print(center_coord)
    df_xyz['x'] = df_xyz['x'] - center_coord[0]
    df_xyz['y'] = df_xyz['y'] - center_coord[1]
    return df_xyz

@timing
def binData2D(myXYZ, xstart, xend, ystart, yend, nx, ny):
    '''
    Fucntion to bin a scatter point cloud (xyz) into a 2d array
    :param myXYZ: xyz array containings the point cloud coordiantes
    :param xstart:
    :param xend:
    :param ystart:
    :param yend:
    :param nx: number of cells along the x-axis
    :param ny: number of cells along hte y-axis
    :return: a group object (pandas library) with all points classified into bins
    '''
    # note, the division requires:     from _future_ import division
    x = myXYZ[:,0].ravel()
    y = myXYZ[:,1].ravel()
    z = myXYZ[:,2].ravel()
    df = pd.DataFrame({'X' : x , 'Y' : y , 'Z' : z})
    bins_x = np.linspace(xstart, xend, nx+1)
    x_cuts = pd.cut(df.X,bins_x, labels=False)
    bins_y = np.linspace(ystart,yend, ny+1)
    y_cuts = pd.cut(df.Y,bins_y, labels=False)
    bin_xmin, bin_ymin = x_cuts.min(), y_cuts.min()
    print('Data cut in a ' + str(bins_x.__len__()) + ' by ' + str(bins_y.__len__()) + ' matrix')
    dx = (xend - xstart)/nx
    dy = (yend - ystart)/ny
    print 'dx = ' + str(dx) + ' ; dy = ' + str (dy)
    grouped = df.groupby([x_cuts,y_cuts])
    print('Data grouped, \nReady to go!!')
    return grouped, bins_x, bins_y, int(bin_xmin), int(bin_ymin)

#=====================================================================
#=====================================================================
#       Function in PROGRESS !!!  Use at your own risk
#===================================================================
@timing
def binData3D(xyz,xstart, xend, ystart, yend, zstart, zend,nx,ny,nz):
    # not ready !!!!
    x = xyz[:,0].ravel()
    y = xyz[:,1].ravel()
    z = xyz[:,2].ravel()
    df = pd.DataFrame({'X' : x , 'Y' : y , 'Z' : z})
    bins_x = np.linspace(xstart,xend,nx)
    x_cuts = pd.cut(df.X,bins_x, labels=False)
    bins_y = np.linspace(ystart,yend,ny)
    y_cuts = pd.cut(df.Y,bins_y, labels=False)
    bins_z = np.linspace(zstart, zend, nz)
    z_cuts = pd.cut(df.Z,bins_z, labels=False)
    print('Data cut in a ' + str(bins_x.__len__()) + ' by ' + str(bins_y.__len__()) + ' by ' + str)(bins_z.__len__()) + ' matrix'
    dx = (xend-xstart)/nx
    dy = (yend - ystart)/ny
    dz = (zend - zstart)/nz
    print('dx = ' + str(dx) + ' ; dy = ' + str (dy) + ' ; dz = ' + str (dz))

    # create a 3D array
    my3d = np.zeros((len(x_cuts),len(y_cuts),len(z_cuts))) * np.nan

    # for loop through the vertical cuts of the poitn clouf to extrac 2d array for each
    for i in np.arange(z_cuts.min(),z_cuts.max()):
        subdf = df[z_cuts==i]
        # group layer into false and true depnding if presence of points or not
        grouped = subdf.groupby([x_cuts,ycuts]).filter(lambda x: np.shape(x)[0]>=1)

        #unstack group into a 2d array
        z = grouped.Z
        
        # add 2d array to 3d array
        my3d[:,:,i] = my2d
    print('Data grouped, \nReady to go!!')
    return my3d #3D array

# def cart2sphere():
#     # write function to convert xyz point coordinates to spehrical coordiantes
#
# def sphere2cart():
#     # write the reverse operation from cart2sphere()