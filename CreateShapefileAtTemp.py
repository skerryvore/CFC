
import osgeo.ogr as ogr
import osgeo.osr as osr

# use a dictionary reader so we can access by field name
reader = csv.DictReader(open("volcano_data.txt","rb"),
        delimiter='\t',
        quoting=csv.QUOTE_NONE)

# set up the shapefile driver
driver = ogr.GetDriverByName("ESRI Shapefile")

# create the data source
data_source = driver.CreateDataSource("Temp.shp")

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
layer.CreateField(ogr.FieldDefn("Temperature", ogr.OFTInteger))

# Process the text file and add the attributes and features to the shapefile
for row in reader:
    # create the feature
    feature = ogr.Feature(layer.GetLayerDefn())
    # Set the attributes using the values from the delimited text file
    feature.SetField("Name", row['Name'])
    feature.SetField("Region", row['Region'])
    feature.SetField("Latitude", row['Latitude'])
    feature.SetField("Longitude", row['Longitude'])
    feature.SetField("Elevation", row['Elev'])

    # create the WKT for the feature using Python string formatting
    wkt = "POINT(%f %f)" %  (float(row['Longitude']) , float(row['Latitude']))

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

