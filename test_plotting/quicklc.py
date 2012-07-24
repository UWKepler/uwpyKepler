import uwpyKepler as kep
import sys
import numpy as num
import pylab
import os

def createLC(kid, **kwargs):
    lc = kep.keplc.keplc(kid)
    kw = quickKW(**kwargs)
    lc.runPipeline(kw)
    return lc

def quickKW(**kwargs):
    ctype='LC'
    gapSize=2
    oWin=15
    oThreshold=5
    dWin=150
    dPolyorder=6
    agap=2
    durfac=2
    detChoice = 'polynomial'
    
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
        elif key == 'detChoice':
            detChoice = kwargs[key]
    kw = kep.keplc.kw(\
    ctype=ctype,\
    gsize=gapSize,\
    owin=oWin,\
    othresh=oThreshold,\
    dwin=dWin,\
    polyorder=dPolyorder,\
    agap=agap,\
    durfac=durfac,\
    detChoice=detChoice)
    return kw

def getCarterQatsData(kid):
    lc = createLC(kid)
    return lc.lcFinal['ydt'], lc.lcFinal['cadence']