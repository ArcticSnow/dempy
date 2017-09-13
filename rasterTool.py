# -*- coding: utf-8 -*-
"""
Created on Thu May 15 19:53:00 2014

@author: simonfilhol
"""

import gdal, sys, osr, os, time
from gdalconst import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cv2

# function to load a raster file
def openRaster(InPath):
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
def saveArray2rasterTif(fname, array, rasterGeotransform, OutPath, _FillValue=-9999, epsg=32632):
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

    cwd = os.getcwd()
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
    os.chdir(cwd)

def fill_nodata(inPath, fileIn, filledPath, fileOut=None):
    '''
    Function to run the Fill_nodata algorithm of GDAL.

    WARNING: This function uses os.system that many do not recomend because of security reason. Anything entered within os.system() will be run into the comand line.

    !!!!!!!!!!!!!!!!!!!!!!!!!USE AT YOUR OWN RISK!!!!!!!!!!!!!!!!!!!

    :param inPath:
    :param fileIn:
    :param filledPath:
    :param fileOut:
    :return:
    '''
    if fileOut is None:
        fileOut = fileIn[:-4] + '_f.tif'

    os.system('gdal_fillnodata.py -md 5 -nomask -of GTiff ' + inPath + fileIn + ' ' + filledPath + fileOut)



def pad_nan_raster(myraster, newXmin, newYmax, newNx, newNy, fname, OutPath, _fill_na=-9999, epsg=32632):
    '''
    Function to pad raster with NAN
    :param myraster:
    :param newXmin:
    :param newYmax:
    :param newNx:
    :param newNy:
    :param fname:
    :param OutPath:
    :param _fill_na:
    :param epsg:
    :return:
    '''
    geoT = raster2array(myraster)[3]
    myarray = raster2array(myraster)[0]
    print '==================='
    print geoT
    print myarray.shape

    geoT = np.round(np.array(geoT), 1)

    if geoT[0] > newXmin:
        myarray = np.concatenate((np.full([myarray.shape[0], np.int(np.round((np.abs(newXmin - geoT[0])) / geoT[1]))], np.nan), myarray), axis=1)
        geoT[0] = newXmin

    if geoT[0] == newXmin:
        if myarray.shape[1] < newNx:
            myarray = np.concatenate((myarray, np.full([myarray.shape[0], newNx - myarray.shape[1]], np.nan)), axis=1)

    if geoT[3] < newYmax:
        myarray = np.concatenate(
            (np.full([np.int(np.round((newYmax - geoT[3]) / geoT[1])), myarray.shape[1]], np.nan), myarray), axis=0)
        geoT[3] = newYmax

    if geoT[3] == newYmax:
        if myarray.shape[1] < newNy:
            myarray = np.concatenate((myarray, np.full([newNy - myarray.shape[0], myarray.shape[1]], np.nan)), axis=0)

    print myarray.shape

    saveArray2rasterTif(fname, myarray, makeGeotransform(geoT[0], geoT[1], geoT[3], -geoT[5]), OutPath, _FillValue=_fill_na, epsg=epsg)




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


def get_pt_value_from_df(rasterMat, geoTransform, df_XsYs, colName='sample'):
    '''

    :param rastermat:
    :param geoTransform:
    :param df_XsYs:
    :return:
    '''
    df_XsYs = df_XsYs.drop(df_XsYs[(((df_XsYs["X"]-geoTransform[0])/geoTransform[1])<=0)|(((df_XsYs["X"]-geoTransform[0])/geoTransform[1])>=len(rasterMat[0]))].index)
    df_XsYs = df_XsYs.drop(df_XsYs[(((df_XsYs["Y"]-geoTransform[2])/geoTransform[3])<=0)|(((df_XsYs["Y"]-geoTransform[2])/geoTransform[3])>=len(rasterMat))].index)
    x = (df_XsYs["X"]-geoTransform[0])/geoTransform[1]
    y = (df_XsYs["Y"]-geoTransform[2])/geoTransform[3]


    sample=pd.DataFrame(rasterMat[y.astype('int'), x.astype('int')], columns=["slope"])
    df_sampled = pd.concat([x,y,sample],axis=1)

    df[colName]
    df_sampled.columns = ['Xs', 'Ys', 'sample']

    return df_sampled



def get_pt_value_from_df_guillaume(rastermat, gt, df):
    '''

    :param rastermat:
    :param gt:
    :param df:
    :return:
    '''
    data = rastermat
    df = df.drop(df[(((df["X"]-gt[0])/gt[1])<=0)|(((df["X"]-gt[0])/gt[1])>=len(data[0]))].index)
    df = df.drop(df[(((df["Y"]-gt[2])/gt[3])<=0)|(((df["Y"]-gt[2])/gt[3])>=len(data))].index)
    x = (df["X"]-gt[0])/gt[1]
    y = (df["Y"]-gt[2])/gt[3]
    rast_min=pd.DataFrame(data[:,:,0][y.astype('int'), x.astype('int')], columns=["Z_min"])
    rast_max=pd.DataFrame(data[:,:,1][y.astype('int'), x.astype('int')], columns=["Z_max"])
    rast_mean=pd.DataFrame(data[:,:,2][y.astype('int'), x.astype('int')], columns=["Z_mean"])
    rast_med=pd.DataFrame(data[:,:,3][y.astype('int'), x.astype('int')], columns=["Z_med"])
    rast_std=pd.DataFrame(data[:,:,4][y.astype('int'), x.astype('int')], columns=["std"])
    rast_count=pd.DataFrame(data[:,:,5][y.astype('int'), x.astype('int')], columns=["count"])
    rast_slope=pd.DataFrame(data[:,:,6][y.astype('int'), x.astype('int')], columns=["slope"])
    return df,pd.concat([rast_min,rast_max,rast_mean,rast_med,rast_std,rast_count,rast_slope],axis=1)



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

def rotateArray(mat, angle):
    '''
    Function to rotate a 2D matrix of a given angle
    :param mat: 2d input matrix
    :param angle: angle in
    :return:
    '''
    mat_center = tuple(np.array(mat.shape)/2)
    rot_mat = cv2.getRotationMatrix2D(mat_center,angle,1.0)
    result = cv2.warpAffine(mat, rot_mat, mat.shape,flags=cv2.INTER_LINEAR)
    return result


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





