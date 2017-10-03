__author__ = 'svfilhol'

import pandas as pd
import scipy.spatial as sp
import pyproj as Proj
import geopandas
from shapely.geometry import Point

# Function open shapefile
def loadSHP(filename):
    return ogr.Open(filename)

# Function to extract coordinates of points in the shapefile (valid for point shapefile)
def extract_coord(shp):
    layer = shp.GetLayer
    points = []
    for index in xrange(layer.GetFeatureCount()):
        feature = layer.GetFeature(index)
        geometry = feature.GetGeometryRef()
        points.append((geometry.GetX(),geometry.GetY()))
        return pd.DataFrame(points)

# Function to extract fields attributes pf points from the shapefile (valid fro point shapefile)
def extract_fields(shp):
    layer = shp.GetLayer
    fields = []
    for index in xrange(layer.GetFeatureCount()):
        feature = layer.GetFeature(index)
        f = feature.items()
        fields.append(f)

# Function to extract and compile data from shapefile into a usable pandas dataframe
def as_df(shp):
    points = extract_coord(shp)
    fields = extract_fields(shp)
    df = [points, fields]
    df.colmns = ['x', 'y', 'Specie', 'id']
    return df

def latlong2TUM(lat,long,epsg=32632):
    '''
    Function to convert lat long point coordinate to UTM (epsg)

    :param lat: array containing the latitudes
    :param long: array containing the latitudes
    :param epsg: epsg code of the final projection
    :return: Xs, Ys in array format
    '''
    pUTM = Proj.Proj(init='epsg:'+str(epsg))
    Xs, Ys = pUTM(lat, long)
    return Xs, Ys


def UTM2latlong(Xs, Ys, epsg=32632):
    '''
    Function to convert UTM (epsg) point coordinate to lat long

    :param Xs: array containing the latitudes
    :param Ys: array containing the latitudes
    :param epsg: epsg code of the original projection
    :return: long, lat in array format
    '''
    pUTM = Proj.Proj(init='epsg:' + str(epsg))
    long, lat = pUTM(Xs, Ys, inverse=True)
    return lat, long

#============== Working Zone =========
def p2p_dist(df):
    '''
    Function to derive distances of all pairs of points from delaunay (to be done)

    :param df:
    :return:
    '''
    sp.spatial.delaunay()
    #1 find all pairs
    triang = sp.Delaunay(xy)
    pairs = []
    for ind in triang.simplices:
        pair1 = ind[0,1]
        pair2 = ind[1,2]
        pair3 = ind[0,2]

    sp.distance.pdist(pairs,'euclidean')