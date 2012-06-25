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

def getCarterQatsData(kid):
    lc = createLC(kid)
    return lc.lcFinal['ydt'], lc.lcFinal['cadence']