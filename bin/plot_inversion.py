#!/usr/bin/python

import numpy as np
import pylab as plt
from matplotlib.collections import LineCollection
import argparse

__author__ = "Romain Beucher"


def helfrag_plot(filename="./RUN00/NA_HELFRAG_results_step1.txt", threshold=150,range_time=[0,500], range_temp=[0,150]):

    data=np.genfromtxt(filename)
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
    
    ax = plt.axes()
    ax.set_xlim(0,500)
    ax.set_ylim(0,150)
    line_segments=LineCollection(v,colors=colors)
    ax.add_collection(line_segments)
    
    
    x = data[-1,1:4]
    x = np.append(x,0)
    y = data[-1,4:8]
    print(abs(x- max(range_time)))
    print(y)
    plt.plot(x,y,color="red", linewidth=2)
    
    plt.xlim(range_time)
    plt.ylim(range_temp)
    plt.xlabel(r'Time (Ma)')
    plt.ylabel(r'Temperature $^{\circ}C$')
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    

def NA_convergence(filename="./RUN00/NA_HELFRAG_results_step1.txt"):

    data=np.genfromtxt(filename)
    plt.plot(data[:,1])

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

    helfrag_plot(filename, threshold, range_time, range_temp)
    plt.show()

