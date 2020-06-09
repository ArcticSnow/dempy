# -*- coding: utf-8 -*-
"""
Created on Thu May 15 2020

@author: simonfilhol

Tool box to compute, analyze and work with spatial observation of snow
"""



import pandas as pd
import numpy as np
import gdal
import scipy.interpolate as interp


###############################################################################################################
#                            SNOW from a Drone
###############################################################################################################





###############################################################################################################
#                            SNOW from a dGPS track
###############################################################################################################

# Function to aggregate non uniform GPS points to a regular grid
def aggregate_SD_track(df, extent=[419347, 420736, 6716038, 6717426], spatialResolution=0.2):
    '''
    Function to resample/aggregate snow depth derived from dGPS track (inherently not homogenous in space) at a given spatial resolution.
    :param df: dataframe containing the at least the following columns ['East', 'Elev', 'North', 'SD']
    :param extent: extent to consider, drop all measurment outside this range. [xmin, xmax, ymin, ymax]
    :param spatialResolution: spatial resolution in meter
    :return: a resampled version of the original dataframe
    '''
    E_bin = np.arange(extent[0], extent[1], spatialResolution)
    N_bin = np.arange(extent[2], extent[3], spatialResolution)
    
    dates = np.unique(df.Date)
    
    resampled = pd.DataFrame()
    for i, date in enumerate(dates):
        tmp = df.loc[(df.Date == date) & ((df.East > extent[0])|(df.East<extent[1])|(df.North>extent[2])|(df.North < extent[3]))]
        
        if tmp.shape[0]>0:
            E_cuts = pd.cut(tmp.East, E_bin, labels=False)
            N_cuts = pd.cut(tmp.North, N_bin, labels=False)

            tmpp = tmp.groupby([E_cuts, N_cuts]).mean()
            tmpp.index.set_names(('index_cut1','index_cut2'), inplace=True)
            tmpp.reset_index(inplace=True)
            tmpp['Nb_samples'] = tmp.groupby([E_cuts, N_cuts])['SD'].count().values
            tmpp['STD_SD'] = tmp.groupby([E_cuts, N_cuts])['SD'].std().values
            
            tmpp.drop(columns=['index_cut1','index_cut2'], inplace=True)

            resampled = resampled.append(tmpp)
    return resampled

def extract_sd_from_point(df_point, ref_dem_file, resol_agg=None):
    '''
    Function to compute snow depth from a dGPS point measurement in reference to a given DEM
    
    :param df_point: Pandas dataframe with East, North, and Elev columns of the measurements
    :param ref_dem_file: reference dem filename
    :param resol_agg: resolution at which final snowtdepth are aggregated. By default resol_agg takes the resolution of the reference DEM
    :retrun:  pandas dataframe with the columns East, Norht, Elev, Elev_ref, SD, STD_SD, Nb_samples. Each row is the value aggregated to the resolution indicated as parameter (resol_agg)
    '''
    
    # Open reference DEM
    dem_mic = ref_dem_file
    myRaster=gdal.Open(dem_mic)
    geot = myRaster.GetGeoTransform()
    Xsize=myRaster.RasterXSize
    Ysize=myRaster.RasterYSize
    data=myRaster.ReadAsArray(0, 0, Xsize, Ysize)

    # define extent and resoltuion from geotiff metadata
    extent = [geot[0], geot[0] + np.round(geot[1],2)*Xsize, geot[3], geot[3] + np.round(geot[5],2)*Ysize]
    if resol_agg is None:
        resol_agg = np.round(geot[1],2)
    
    # Create the X,Y coordinate meshgrid
    Xs = np.arange(extent[0],extent[1], np.round(geot[1],2))
    Ys = np.arange(extent[2], extent[3], np.round(geot[5],2))
    XX, YY = np.meshgrid(Xs, Ys)
    
    # make mask for the ROI based on extent provided above
    mask = np.copy(XX) * 0
    mask[(XX> extent[0]) & (XX < extent[1]) & (YY > extent[2])& (YY < extent[3])] = 1
    mask = mask.astype(bool)
    XY = np.vstack((XX[mask],YY[mask])).T 
    Z = data[mask]
    
    # crop dGPS data to area of interest (extent)
    df_point.loc[(df_point.East < extent[0])|(df_point.East > extent[1])|(df_point.North < extent[2])|(df_point.North > extent[3])] = np.nan
    df_point.dropna(inplace=True)
        
    # Linear interpolation of the reference DEM at the location of the dGPS measurement
    df_point['Elev_ref'] = interp.griddata(XY, Z, (df_point.East, df_point.North), method='linear')

    # compute snowdepth
    df_point['SD'] = df_point.Elev - df_point.Elev_ref

    # aggregate values to wanted resolution
    df_gps_resample = aggregate_SD_track(df_point, extent=extent, spatialResolution=resol_agg)
    
    return df_gps_resample
