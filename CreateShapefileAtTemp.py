#!/usr/bin/python

import osgeo.ogr as ogr
import osgeo.osr as osr
import numpy as np
from scipy import interpolate


def FindTempAtTime(t, T, x):
    f = interpolate.interp1d(t, T, bounds_error=False, fill_value=-9999)
    Tx = f(x)
    return Tx

ti = 400.0
inputFile = "KingsCanyon_out.txt"

f = open(inputFile,"r")
# Skip header
header = f.readline()

listofsamples = []
for line in f:
    line = line.strip()
    columns = line.split()
    source = {}
    source['Name'] = columns[0]
    source['Misfit'] = float(columns[1])
    source['Latitude'] = float(columns[2])
    source['Longitude'] = float(columns[3])
    source['Time'] = [float(columns[i]) for i in range(4,8)]
    source['Temperature'] = [float(columns[i]) for i in range(8,12)]
    listofsamples.append(source)


# set up the shapefile driver
driver = ogr.GetDriverByName("ESRI Shapefile")

# create the data source
data_source = driver.CreateDataSource("Temperatures.shp")

# create the spatial reference, WGS84
srs = osr.SpatialReference()
srs.ImportFromEPSG(4326)

# create the layer
layer = data_source.CreateLayer("Temperatures", srs, ogr.wkbPoint)

# Add the fields we're interested in
field_name = ogr.FieldDefn("Name", ogr.OFTString)
field_name.SetWidth(24)
layer.CreateField(field_name)
layer.CreateField(ogr.FieldDefn("Latitude", ogr.OFTReal))
layer.CreateField(ogr.FieldDefn("Longitude", ogr.OFTReal))
layer.CreateField(ogr.FieldDefn("Temp", ogr.OFTReal))
layer.CreateField(ogr.FieldDefn("CoolRate", ogr.OFTReal))

# Process the text file and add the attributes and features to the shapefile
for sample in listofsamples:
    # create the feature
    feature = ogr.Feature(layer.GetLayerDefn())
    # Set the attributes using the values from the delimited text file
    feature.SetField("Name", sample['Name'])
    feature.SetField("Latitude", sample['Latitude'])
    feature.SetField("Longitude", sample['Longitude'])
    # Get Temperature at time ti
    Temperature = FindTempAtTime(sample['Time'], sample['Temperature'],ti) 
    feature.SetField("Temp", float(Temperature))
    # Calculate Cooling Rate
    T1 = FindTempAtTime(sample['Time'], sample['Temperature'],ti-1.0) 
    T2 = FindTempAtTime(sample['Time'], sample['Temperature'],ti+1.0) 
    dT = T2 - T1
    dt = 2.0
    feature.SetField("CoolRate", float(dT/dt))

    # create the WKT for the feature using Python string formatting
    wkt = "POINT(%f %f)" %  (float(sample['Longitude']) , float(sample['Latitude']))

    # Create the point from the Well Known Txt
    point = ogr.CreateGeometryFromWkt(wkt)

    # Set the feature geometry using the point
    feature.SetGeometry(point)
    # Create the feature in the layer (shapefile)
    layer.CreateFeature(feature)
    # Destroy the feature to free resources
    feature.Destroy()

# Destroy the data source to free resources
data_source.Destroy()

