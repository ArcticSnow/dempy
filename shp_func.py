__author__ = 'svfilhol'

import pandas as pd
import scipy.spatial as sp

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