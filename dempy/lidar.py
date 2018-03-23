from __future__ import division
import liblas.file as lf
import numpy as np
import pandas as pd
import gdal, osr


def openlas(fname):
    '''
    Function to open laz/las file
    :param fname:
    :return:
    '''
    f = lf.File(fname, mode='r')

    # Get number of points from header
    num_points = int(f.__len__())
    # Create empty numpy array
    PointsXYZIC = np.empty(shape=(num_points, 5))

    for i, p in enumerate(f):
        newrow = [p.x, p.y, p.z, p.intensity, p.classification]
        PointsXYZIC[i] = newrow
    p = f.read(0)
    return xyzrc

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
    # note, the division requires:     from __future__ import division
    x = myXYZ[:, 0].ravel()
    y = myXYZ[:, 1].ravel()
    z = myXYZ[:, 2].ravel()
    df = pd.DataFrame({'X' : x , 'Y' : y , 'Z' : z})
    bins_x = np.linspace(xstart,xend,nx)
    x_cuts = pd.cut(df.X, bins_x, labels=False)
    bins_y = np.linspace(ystart,yend,ny)
    y_cuts = pd.cut(df.Y, bins_y, labels=False)
    print('Data cut in a ' + str(bins_x.__len__()) + ' by ' + str(bins_y.__len__()) + ' matrix')
    dx = (xend - xstart)/nx
    dy = (yend - ystart)/ny
    print('dx = ' + str(dx) + ' ; dy = ' + str (dy))
    grouped = df.groupby([x_cuts,y_cuts])
    print('Data grouped, \nReady to go!!')
    return grouped

def plas2raster(plas, xstart, xend, ystart, yend, nx=1000, ny=1000, rasterFname='myraster.tif', method='min', _FillValue=-9999, epsg=32632):
    xyz = np.vstack(plas.x, plas.y, plas.z)
    binned = binData2D(myXYZ, xstart, xend, ystart, yend, nx, ny)

    if method == 'min':
        ret = binned.Z.min().unstack()
    elif method == 'max':
        ret = binned.Z.max().unstack()
    elif method == 'mean':
        ret = binned.Z.mean().unstack()
    elif method == 'median':
        ret = binned.Z.median().unstack()
    elif method == 'count':
        ret = binned.Z.count().unstack()

    # put in gdal raster format
    rasterGeotransform = [xstart, np.round((xend - xstart)/nx), 0, yend, 0, -np.round((yend-ystart)/ny)]
    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(rasterFname, nx, ny, 1, gdal.GDT_Float64)
    outRaster.SetGeoTransform(rasterGeotransform)
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(_FillValue)
    outband.WriteArray(ret)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

    return outRaster

fname = '/home/arcticsnow/Github/Finse_lidar/data/ScanPos003-SINGLESCANS-170324_103651.las'
mylas = openlas(fname)

print(np.min(mylas.x))
print(np.max(mylas.x))
print(np.min(mylas.y))
print(np.max(mylas.y))


raster = plas2raster(mylas, )








