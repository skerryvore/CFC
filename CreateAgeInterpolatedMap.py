#!/usr/bin/pytho    return final_age.value, fmean.value, fdist[0:200]n
import gdal
import numpy as np
from scipy import interpolate
from ctypes import *

def ReadDEMasArray(dem_file):
    ds = gdal.Open(dem_file)
    dem = ds.ReadAsArray()
    ds = None
    return dem

def FindTempAtTime(t, T, x):
    f = interpolate.interp1d(t, T, bounds_error=False, fill_value=-9999)
    Tx = f(x)
    return Tx

def GetDEMGrid(dem):
    ds = gdal.Open(dem)
    gt = ds.GetGeoTransform()
    ds = None
    xres = gt[1]
    yres = gt[5]
    X = np.linspace(gt[0], gt[0] + dem.shape[1] * xres, dem.shape[1])
    Y = np.linspace(gt[3], gt[3] + dem.shape[0] * yres, dem.shape[0])
    X, Y = np.meshgrid(X, Y)
    return X, Y

def GetKetcham(time,temp):
    dll = cdll.LoadLibrary("./libketcham.so")
    ntime = c_int(len(time))
    time = (c_float * ntime)(*time)
    temp = (c_float * ntime)(*temp)
    alo = c_double(16.3)
    # outputs
    final_age = c_double()
    oldest_age = c_double()
    fmean = c_double()
    fdist = [0]*200
    fdist = (c_double * 200)(*fdist)
    redDensity = c_double()
    dll.ketch_main_(byref(ntime), ketchtime, ketchtemp, byref(alo), byref(final_age),
                    byref(oldest_age), byref(fmean), fdist, byref(redDensity))
    return final_age.value, fmean.value, fdist[0:200]

time = [100,80,20,0]
temp = [150,80,40,0]


geogradient = 14.0
t0 = 100.0
dt = 1.0

t = np.arange(100.0,0.0,-1.0)

original_dem = ReadDEMasArray("./DEM/KingsCanyon.tif")
filtered_dem = ReadDEMasArray("./DEM/KingsCanyon_filtered2.tif")
offset_dem = original_dem - filtered_dem

inversions = np.genfromtxt("KingsCanyon_out.txt", skip_header=1, dtype=None)


for i in t:
    newtime = FindTempAtTime(time, temp, i)
    print(newtime)


# Call Ketcham routine for all pixels...
