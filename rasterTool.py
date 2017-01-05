# -*- coding: utf-8 -*-
"""
Created on Thu May 15 19:53:00 2014

@author: simonfilhol
"""

import gdal, sys, osr, os, time
from gdalconst import *
from tkinter import *
import tkinter.filedialog as tkFileDialog
import matplotlib.pyplot as plt
import numpy as np

# function to load a raster file
def openRaster(InPath=None):
    '''
    openRaster()   Load a raster file into Python console using gdal libraries.\n
    **myRaster, InPath = openRaster("filepath")**\n
    Copyright S Filhol \n
    Dependencies: Tkinter, gdal, sys, os  
    
    Parameters
    ----------    
    **InPath** :  string indicating the file path including the filename      

    Returns
    -------
    **myRaster** :      Raster class object (gdal object) \n
    **InPath** : string indicating the file path including the filename
    '''
    if InPath == None:
        root = Tk()
        root.withdraw()
        InPath = tkFileDialog.askopenfilename()
        time.sleep(.9)
        root.destroy()
        root.quit()
    if InPath is None:
        print('Could not open image')
        sys.exit(1)
    print(InPath)
    myRaster=gdal.Open(InPath)
    return myRaster
    
# function to convert raster data into an array ready for processing
def raster2array(myRaster):
    '''
    raster2array()   Convert Gdal raster object to  2 array. The raster object has to have only one band\n
    **data, dx, dy, transform = raster2array(myRaster)**\n
    Copyright S Filhol \n
    Dependencies: gdal  
    
    Parameters
    ----------    
    **myRaster** :  Raster class object (gdal object) \n    

    Returns
    -------
    **data** :      2D array containing the raster band data \n
    **dx, dy** : spatial resolution of th raster \n
    **transform** : raster transformation information. Very usefull information if you save the data back to raster form or plot data
    '''
    transform = myRaster.GetGeoTransform()
    dx=transform[1]
    dy=transform[5]
    Xsize=myRaster.RasterXSize
    Ysize=myRaster.RasterYSize
    data=myRaster.ReadAsArray(0, 0, Xsize, Ysize)
    myRaster=None
    return data, dx, dy, transform


def makeGeotransform(Xmin, dx, Ymax, dy):
    '''

    :param Xmin: x-coordinate of the upper left corner
    :param dx: pixel resolution x-direction
    :param Ymax: y-coordiante of the upper left corner
    :param dy: pixel resolution y-direction
    :return: a geotranform object ready to use in GDAL functions

    Note: more info http://www.perrygeo.com/python-affine-transforms.html
    '''
    return [Xmin, dx, 0, Ymax, 0, -dy]


# function to save results as a geotiff raster file
def saveArray2rasterTif(fname, array, rasterGeotransform, _FillValue=-9999, OutPath=None, epsg=32632):
    '''
    Save to a GeoTiff file the array\n
    **saveArray2rasterTif(filename, transform, myArray, OutPath)**
    Copyright S Filhol
    Dependencies: gdal, os, Tkinter, osr

    :param fname:  string of the new file name
    :param array: 2D matrix containing the data to save as a raster file. If orientation of array in final raster is lipped try using    np.flipud(array)  as input
    :param rasterGeotransform:
    :param _FillValue: default -9999
    :param OutPath: optional. String indicating the path where to save the file \n
    :param epsg: projecton epsg code. default is 32632 for finse UTM 32N. 32606 for 6N (northern Alaska)
    :return:
    '''
    
    if OutPath==None:
        root = Tk()
        root.withdraw()
        OutPath = tkFileDialog.askdirectory(parent=root)
        time.sleep(0.9)
        root.destroy()
        root.quit()
    os.chdir(OutPath)
    cols = array.shape[1]
    rows = array.shape[0]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(fname, cols, rows, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform(rasterGeotransform)
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(_FillValue)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def plot_raster(raster, band):
    mat = raster.GetRasterBand(band).ReadAsArray()
    plt.figure()
    plt.imshow(mat)
    plt.colorbar()
    plt.show()


def get_pt_value_rasterfile(rasterfile, Xs, Ys):
  gdata = gdal.Open(rasterfile)
  gt = gdata.GetGeoTransform()
  data = gdata.ReadAsArray().astype(np.float)
  gdata = None
  x = (Xs - gt[0])/gt[1]
  y = (Ys - gt[3])/gt[5]
  return data[y.astype('int'), x.astype('int')]

def get_pt_value_array(myarray, geotransform, Xs, Ys):
  x = (Xs - geotransform[0])/geotransform[1]
  y = (Ys - geotransform[3])/geotransform[5]
  return myarray[y.astype('int'), x.astype('int')]


def get_pt_value_raster(myraster, Xs, Ys):
  gt = myraster.GetGeoTransform()
  data = myraster.ReadAsArray().astype(np.float)
  gdata = None
  x = (Xs - gt[0])/gt[1]
  y = (Ys - gt[3])/gt[5]
  return data[y.astype('int'), x.astype('int')]


def plot_array(mat):
    plt.figure()
    plt.imshow(mat, origin='lower')
    plt.colorbar()
    plt.show()

def plot_hist(map,nbin):
    plt.figure()
    plt.hist(map.flatten(), bins=nbin)
    plt.show()


def detrend(M):
    '''
    **D, coeff = detrend(M)**

    fits a plane to the surface defined by the 2D array M and detrends the surface by subtracting the value of the planar fit at each element.
    Returns the detrended surface in 2D array D.

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

    nx = np.size(M, 0)
    ny = np.size(M, 1)
    [Y, X] = np.mgrid[0:nx, 0:ny]

    if ny * nx > 10000:
        inds = np.zeros((nx, ny), dtype='int8')
        inds = inds.flatten()
        inds[0:9000] = 1
        np.random.shuffle(inds)
        A = X.flatten()
        A = A[inds == 1]
        B = Y.flatten()
        B = B[inds == 1]
        C = M.flatten()
        C = C[inds == 1]
        ABC = np.array([A, B, C])
        points = np.reshape(ABC, (np.shape(ABC)[0], -1))
        assert points.shape[0] < points.shape[1]
        ctr = points.mean(axis=1)
        x = points - ctr[:, None]
        Mat = np.dot(x, x.T)
        norm = np.linalg.linalg.svd(Mat)[0][:, -1]

    else:
        XYM = np.array([X.flatten(), Y.flatten(), M.flatten()])
        points = np.reshape(XYM, (np.shape(XYM)[0], -1))
        assert points.shape[0] < points.shape[1]
        ctr = points.mean(axis=1)
        x = points - ctr[:, None]
        Mat = np.dot(x, x.T)
        norm = np.linalg.linalg.svd(Mat)[0][:, -1]
    # for a plane with equation z = ax + by + c
    # at each (x,y) point, subtract the value of the fitted plane from dem
    d = -norm[0] * ctr[0] - norm[1] * ctr[1] - norm[2] * ctr[2]
    PlaneFit = -(norm[0] * X + norm[1] * Y + d) / norm[2]
    D = M - PlaneFit
    coeff = np.array([norm, d])
    print ('Matrix detrended!')
    return D, coeff
#======================= Test Zone =======================================

'''
myArray,transform=surfaceMaker((1000,1000),(0,100,0,150),"M*sin(X/5)+M*Y")
rasterOrigin=(0,0)
saveArray2rasterTif("dddd2.tif",myArray,rasterOrigin,.1,.15,NAN)
'''





