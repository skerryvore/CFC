#!/usr/bin/python

import gdal
import numpy as np
import pylab as plt
import argparse
from matplotlib.collections import LineCollection

parser = argparse.ArgumentParser(description='Script by Romain Beucher')
parser.add_argument('-f','--filename', help = 'NA Result file',  required=True)
parser.add_argument('-o','--output', help = 'NA Result file',  required=True)
args = parser.parse_args()
filename = args.filename
output = args.output

ax = plt.subplot(111)
ds = gdal.Open('./DEM/KingsCanyon.tif')
dem = ds.ReadAsArray()
gt = ds.GetGeoTransform()
ds = None

xres = gt[1]
yres = gt[5]

X = np.linspace(gt[0], gt[0] + dem.shape[1] * xres, dem.shape[1])
Y = np.linspace(gt[3], gt[3] + dem.shape[0] * yres, dem.shape[0])
ax.set_xlim((min(X), max(X)))
ax.set_ylim((min(Y), max(Y)))

X, Y = np.meshgrid(X, Y)

plt.pcolormesh(X,Y,dem, cmap="gist_earth")

ax.set_ylabel('Latitude')
ax.set_xlabel('Longitude')
ax.ticklabel_format(useOffset=False)


# Get number of samples:
#filename ="./results/KINGSCANYON/NA_results_sample1.txt" 
f = open(filename, 'r')
string = f.readline()
nsamples = int(string.split(sep=":")[-1])
num_lines = sum(1 for line in f)
samples = np.genfromtxt(filename, skip_header=4, skip_footer=num_lines-(nsamples+3))
f.close()

xs = samples[:,2]
ys = samples[:,1]

ax.plot(xs,ys,'ro', markersize=10)
ax.plot(xs[0], ys[0], 'bo', markersize=20)

xpairs = []
ypairs = []
data = np.genfromtxt("list_segments.txt")
for I in range(len(data[:,1])):
        xends = [data[I,0], data[I,2]]
        yends = [data[I,1], data[I,3]]
        xpairs.append(xends)
        ypairs.append(yends)

segs = list(zip(zip(data[:,1],data[:,0]), zip(data[:,3],data[:,2])))
line_segments = LineCollection(segs,linewidths=1)

ax.add_collection(line_segments)



if(output == "file"):
    string = filename
    string2 = string.split(sep=".")[0]
    plt.savefig(string2+".png")

if(output == "screen"):
    plt.show()


