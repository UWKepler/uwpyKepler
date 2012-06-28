import uwpyKepler as kep
import sys
import numpy as num
import pylab
import os

def createLC(kid):
    lc = kep.keplc.keplc(kid)
    kw = kep.keplc.kw(\
    ctype='LC',\
    gsize=2,\
    owin=15,\
    othresh=5,\
    dwin=50,\
    polyorder=6,\
    agap=1,\
    durfac=2)
    lc.runPipeline(kw)
    return lc

def quickKW(**kwargs):
    ctype='LC'
    gapSize=2
    oWin=15
    oThreshold=5
    dWin=50
    dPolyorder=6
    agap=1
    durfac=2
    for key in kwargs:
        if key == 'ctype':
            ctype = kwargs[key]
        elif key == 'gsize':
            gapSize = kwargs[key]
        elif key == 'owin':
            oWin = kwargs[key]
        elif key.lower().startswith('othresh'):
            oThreshold = kwargs[key]
        elif key == 'dwin':
            dWin = kwargs[key]
        elif key == 'polyorder':
            dPolyorder = kwargs[key]
        elif key == 'agap':
            agap = kwargs[key]
        elif key == 'durationfactor':
            durfac = kwargs[key]
    kw = kep.keplc.kw(\
    ctype=ctype,\
    gsize=gapSize,\
    owin=oWin,\
    othresh=oThreshold,\
    dwin=dWin,\
    polyorder=dPolyorder,\
    agap=agap,\
    durfac=durfac)
    return kw

def getCarterQatsData(kid):
    lc = createLC(kid)
    return lc.lcFinal['ydt'], lc.lcFinal['cadence']