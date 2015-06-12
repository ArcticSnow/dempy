# -*- coding: utf-8 -*-
"""
Created on Thu May 15 19:53:00 2014

@author: simonfilhol
"""

import gdal, sys, osr, os, time
from gdalconst import *
from Tkinter import *
import tkFileDialog
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
        print 'Could not open image'
        sys.exit(1)
    print InPath
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

# function to save results as a geotiff raster file
def saveArray2rasterTif(fname, array, rasterGeotransform, _FillValue,OutPath=None):
    '''
    Save to a GeoTiff file the array\n
    **saveArray2rasterTif(filename, transform, myArray, OutPath)** \n
    Copyright S Filhol \n
    Dependencies: gdal, os, Tkinter, osr  
    
    Parameters
    ----------    
    **fname** :  string of the new file name \n
    **array** : 2D matrix containing the data to save as a raster file \n
    **rasterOrigin**
    **OutPath** : optional. String indicating the path where to save the file \n
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
    outRasterSRS.ImportFromEPSG(32606)  #EPGS code for WGS 84/ UTM zone 6N
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def plot_raster(raster, band):
    mat = raster.GetRasterBand(band).ReadAsArray()
    plt.figure()
    plt.imshow(mat)
    plt.colorbar()
    plt.show()


def plot_matrix(mat):
    plt.figure()
    plt.imshow(mat)
    plt.colorbar()
    plt.show()

def plot_hist(map,nbin):
    plt.figure()
    plt.hist(map.flatten(), bins=nbin)
    plt.show()


#======================= Test Zone =======================================

'''
myArray,transform=surfaceMaker((1000,1000),(0,100,0,150),"M*sin(X/5)+M*Y")
rasterOrigin=(0,0)
saveArray2rasterTif("dddd2.tif",myArray,rasterOrigin,.1,.15,NAN)
'''





