#!/usr/bin/python

import numpy as np
import pylab as plt
from matplotlib.collections import LineCollection
import argparse
from itertools import islice

__author__ = "Romain Beucher"


def helfrag_plot(filename="./RUN00/NA_HELFRAG_results_step1.txt", threshold=150,range_time=[0,500], range_temp=[0,150],skip=0):

    data=np.genfromtxt(filename, skip_header=skip+1)
   
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
        if(data[i,0] < threshold): color="blue"
        if(data[i,0] > threshold): color="gray"
        v.append(d)
        colors.append(color)
    
    line_segments=LineCollection(v,colors=colors)
    
    x = data[-1,1:4]
    x = np.append(x,0)
    y = data[-1,4:8]

    return line_segments, x, y


def NA_convergence(filename="./RUN00/NA_HELFRAG_results_step1.txt",skip=0):

    data=np.genfromtxt(filename,skip_header=skip+1)
    return data[:,1]

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
    parser.add_argument('-t','--threshold', help = 'threshold',  required=True)
    parser.add_argument('-rt1','--range_time', help = 'time range',  required=False)
    parser.add_argument('-rt2','--range_temp', help = 'temp range',  required=False)
    args = parser.parse_args()

    filename = args.filename
    threshold = int(args.threshold)
    range_time = args.range_time
    range_temp = args.range_temp


    range_time = [0,500]
    range_temp = [0,150]

    skip = find_line(filename, "List of")

    line_segments, xbest, ybest = helfrag_plot(filename, threshold, range_time, range_temp,skip)

   
    ax = plt.subplot2grid((2,4),(0,0),colspan=2)
    ax.set_xlim(0,500)
    ax.set_ylim(0,150)
    ax.add_collection(line_segments)
    plt.plot(xbest,ybest,color="red", linewidth=2)
    
    plt.xlim(range_time)
    plt.ylim(range_temp)
    plt.xlabel(r'Time (Ma)')
    plt.ylabel(r'Temperature $^{\circ}C$')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
   
    misfits = NA_convergence(filename,skip)
    plt.subplot2grid((2,4),(1,0),colspan=2)
    plt.xlim([0,len(misfits)])
    plt.ylim([0,max(misfits)])
    plt.plot(misfits)



    # Get number of samples:
    f = open(filename, 'r')
    string = f.readline()
    nsamples = int(string.split(sep=":")[-1])
    num_lines = sum(1 for line in f)
    samples = np.genfromtxt(filename, skip_header=2, skip_footer=num_lines-nsamples-1,dtype=None)
    
    ages = []
    elevation = []
    for sample in samples:
        ages.append(sample[4])
        elevation.append(sample[3])

    plt.subplot2grid((2,4),(0,2), colspan=2, rowspan=2)
    plt.plot(ages,'o')

    plt.show()

