import uwpyKepler as kep
import numpy as num
import pylab
import scipy
import sys
import os

def geteDataFromFile(kid):
    fileDir = '/astro/store/student-scratch1/johnm26/dFiles/'
    name = 'eDataDiscoveries.txt'
    dFile = open(fileDir + name, 'r')
    # first line is useless
    dFile.readline()
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            line   = line.split()
            period = float(line[1])
            t0     = float(line[2])
            q      = float(line[3])
            break
    return period, t0, q