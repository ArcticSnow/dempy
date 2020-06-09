'''
Class to work with data stored as raster file. This keeps data and meta together


'''


import gdal
import morphometry

class raster_map(object):
    
    def __init__(self, fname=None):
        if (fname is not None) and (fname[-4:]='.tif'):
            self.fname = fname
            self.gdal_raster = gdal.Open(self.fname)
            self.nx = self.gdal_raster.RasterXSize
            self.ny = self.gdal_raster.RasterYSize
            self.geot = self.gdal_raster.GetGeoTransform()
            self.dx = self.geot[1]
            self.dy = self.geot[5]
            self.values = self.gdal_raster.ReadAsArray(0, 0, self.nx, self.ny)
            self.epsg = 
    
    #################################################################
    #    Morphometry
    
    def slope(self, neighbor=None):
        self.slope = slope(Z=self.values, dx=self.dx, neighbor=neighbor)
        
    def mean_curvature(self):
    
    def gaussian_curvature(self):
    
    ################################################################
    # crop: update geotransform

    # detrend
    
    # decomposition
    
    # Geotstatisitcs
    
    # plotting