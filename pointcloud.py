from __future__ import division
from laspy.file import File
import numpy as np
import pandas as pd
import time

def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print '%s function took %0.3f ms' % (f.func_name, (time2-time1)*1000.0)
        return ret
    return wrap

@timing
def loadLAS2XYZ(filepath):
    print 'Start loading...'
    inFile = File(filepath, mode='r')
    coords = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
    print 'Data loaded'
    return coords

def LAS2txt(filepath,newfile):
    inFile = File(filepath, mode='r')
    coords = np.vstack((inFile.x, inFile.y, inFile.z)).transpose()
    if newfile[-4] != '.txt':
        newfile = newfile + '.txt'
    np.savetxt(newfile,coords)
    print 'File saved: ' + newfile

def xyz_subsample(xyz, length_out):
    ind = np.random.randint(0,xyz.shape[0],length_out)
    xyz_new = xyz[ind,:]
    print 'xyz subsampled!'
    return xyz_new

def xyz_stat(xyz):
    print 'Shape of array: ' + str(xyz.shape)
    print 'Min of xyz: '
    print np.min(xyz, axis=0)
    print 'Max of xyz: '
    print np.max(xyz, axis=0)
    print 'Mean of xyz: '
    print np.mean(xyz, axis=0)
    print 'Extent'
    print np.max(xyz, axis=0)-np.min(xyz, axis=0)

def trans(xyz,trans_vec):
    xyz[:,0] = xyz[:,0] - trans_vec[0]
    xyz[:,1] = xyz[:,1] - trans_vec[1]
    xyz[:,2] = xyz[:,2] - trans_vec[2]
    return xyz

def truncate(xyz, Xextent, Yextent):
    xcut = xyz[xyz[:,0]>=Xextent[0]]
    xcut1 = xcut[xcut[:,0]<Xextent[1]]
    ycut = xcut1[xcut1[:,1]>=Yextent[0]]
    ycut1 = ycut[ycut[:,1]<Yextent[1]]
    return ycut1

def get_slice(xyz,):


def translate_coords(coords, xyz_trans = None ,ask = True):
    if xyz_trans is None:
        xyz_trans = [coords[:,0].min(), coords[:,1].min(), coords[:,2].min()]
    if ask is True:
        print 'Default translation:'
        print str(xyz_trans) + '\n'
        res = input('Do you want to translate? 0/1')
        if res is 0:
            print 'No Translation applied'
            return None
        if res is 1:
            return trans(coords, xyz_trans)
    if ask is not True:
        return trans(coords, xyz_trans)

@timing
def binData2D(myXYZ, xstart, xend, ystart, yend, nx, ny):
    # note, the division requires:     from __future__ import division
    x = myXYZ[:,0].ravel()
    y = myXYZ[:,1].ravel()
    z = myXYZ[:,2].ravel()
    df = pd.DataFrame({'X' : x , 'Y' : y , 'Z' : z})
    bins_x = np.linspace(xstart,xend,nx)
    x_cuts = pd.cut(df.X,bins_x, labels=False)
    bins_y = np.linspace(ystart,yend,ny)
    y_cuts = pd.cut(df.Y,bins_y, labels=False)
    print 'Data cut in a ' + str(bins_x.__len__()) + ' by ' + str(bins_y.__len__()) + ' matrix'
    dx = (xend - xstart)/nx
    dy = (yend - ystart)/ny
    print 'dx = ' + str(dx) + ' ; dy = ' + str (dy)
    grouped = df.groupby([x_cuts,y_cuts])
    print 'Data grouped, \nReady to go!!'
    return grouped

@timing
def binData3D(xyz,xstart, xend, ystart, yend, zstart, zend,nx,ny,nz):
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
    print 'Data cut in a ' + str(bins_x.__len__()) + ' by ' + str(bins_y.__len__()) + ' by ' + str(bins_z.__len__()) + ' matrix'
    dx = (xend-xstart)/nx
    dy = (yend - ystart)/ny
    dz = (zend - zstart)/nz
    print 'dx = ' + str(dx) + ' ; dy = ' + str (dy) + ' ; dz = ' + str (dz)
    grouped = df.groupby([x_cuts,y_cuts,z_cuts])
    print 'Data grouped, \nReady to go!!'
    return grouped

def voxelPCL(grouped):
    grouped