#!/usr/bin/python

import numpy as np
import pylab as plt
from matplotlib.collections import LineCollection
import argparse
from itertools import islice
import matplotlib.colors as mc
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
import gdal

__author__ = "Romain Beucher"


def helfrag_plot(filename="./RUN00/NA_HELFRAG_results_step1.txt",range_time=[0,500], range_temp=[0,150],skip=0):

    data=np.genfromtxt(filename, skip_header=skip+1)
    reject = 1e21
    data = data[data[:,0] < reject]

    # Sort the data
    data = data[data[:,0].argsort()]
    # Reverse sort
    data = data[::-1]
    
    for i in [1,2,3]:
        data[:,i] = abs(data[:,i] - max(range_time))
    
    v = []
    colors = []
    for i in np.arange(data.shape[0]):
        x = data[i,1:4]
        x = np.append(x,0)
        y = data[i,4:8]
        d = tuple(zip(x,y))
        v.append(d)
    
    line_segments=LineCollection(v)
    
    x = data[-1,1:4]
    x = np.append(x,0)
    y = data[-1,4:8]

    misfits = data[:,0]

    return line_segments, x, y, misfits


def NA_convergence(filename="./RUN00/NA_HELFRAG_results_step1.txt",skip=0):

    data=np.genfromtxt(filename,skip_header=skip+1)
    threshold = 1e21
    num = np.zeros((len(data[:,0]),1))
    num[:,0] = range(1,len(data[:,0])+1)
    data = np.hstack((data,num))
    data = data[data[:,0] < threshold]
    return data[:,0], data[:,-1]

def find_line(filename, string):

    with open(filename, 'r') as myfile:
        for index, line in enumerate(myfile):
            if string in line:
                skip = index
                break
    return skip

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Script by Romain Beucher')
    parser.add_argument('-f','--filename', help = 'NA Result file',  required=True)
    parser.add_argument('-o','--output', help = 'NA Result file',  required=True)
    parser.add_argument('-rt1','--range_time', help = 'time range',  required=False)
    parser.add_argument('-rt2','--range_temp', help = 'temp range',  required=False)
    args = parser.parse_args()

    filename = args.filename
    output = args.output
    range_time = args.range_time
    range_temp = args.range_temp


    range_time = [0,500]
    range_temp = [0,150]

    skip = find_line(filename, "List of") + 3

    line_segments, xbest, ybest, mvals = helfrag_plot(filename, range_time, range_temp,skip)

    line_segments.set_array(mvals)
    line_segments.set_cmap(cm.jet)
    line_segments.set_alpha(0.3)
    line_segments.set_norm(LogNorm())
    ax = plt.subplot()
    ax.set_xlim(0,500)
    ax.set_ylim(0,150)
    ax.add_collection(line_segments)
    fig =plt.gcf()
    axcb = fig.colorbar(line_segments)
    plt.plot(xbest,ybest,color="red", linewidth=2)
    plt.title(filename)

    plt.xlim(range_time)
    plt.ylim(range_temp)
    plt.xlabel(r'Time (Ma)')
    plt.ylabel(r'Temperature $^{\circ}C$')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
   
    a = plt.axes([0.15, 0.67,0.2,0.2], axisbg="w")
    misfits, indexes = NA_convergence(filename,skip)
    plt.plot(indexes, misfits)
    plt.setp(a, xticks=[], yticks=[])

    # Get number of samples:
    f = open(filename, 'r')
    string = f.readline()
    nsamples = int(string.split(sep=":")[-1])
    num_lines = sum(1 for line in f)
    dtypes=[('Name', 'S10'), ('Lat', 'f'), ('Lon', 'f'), ('Elevmeas', 'i'), ('Elevdem', 'i'), ('Elevfilt', 'i'),('FTage', 'f'), ('MTL', 'f')]
    samples = np.genfromtxt(filename, skip_header=4, skip_footer=num_lines-(nsamples+3),dtype=None)

    if(output == "screen"):
        plt.show()

    if(output == "file"):
        string=filename
        string=string.split(sep=".")[0]
        plt.savefig(string+".png")



