__author__ = 'svfilhol'

import pandas as pd

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
def df(shp):
    points = extract_coord(shp)
    fields = extract_fields(shp)
    df = [points, fields]
    df.colmns = ['x', 'y', 'Specie', 'id']
    return df
