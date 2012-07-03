import uwpyKepler as kep
import numpy as num
import pylab
from scipy import optimize
import sys
import os

fileDir = '/astro/store/student-scratch1/johnm26/dFiles/'
name = 'eDataDiscoveries.txt'

def geteDataFromFile(kid):
    dFile = open(fileDir + name, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            line   = line.split()
            period = float(line[1])
            t0     = float(line[2])
            q      = float(line[3])
            return period, t0, q
    return -1, -1, -1
    
def eDataInFile(kid):
    dFile = open(fileDir + name, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            return True
    return False
    
class modelLC:
    def __init__(self, kid, lc, guess):
        self.kid    = kid
        self.xdata  = lc.lcFinal['x'] + lc.BJDREFI - 2454900e0
        self.ydata  = lc.lcFinal['ydt']
        key = lc.eData['KOI'].keys()[0]
        self.period = lc.eData['KOI'][key]['Period']
        self.t0     = lc.eData['KOI'][key]['T0']
        self.q      = lc.eData['KOI'][key]['Duration']
        self.inc    = guess[0]
        self.aRs    = guess[1]
        self.RpRs   = guess[2]
        self.u1     = guess[3]
        self.u2     = guess[4]
        self.setPhase()
        self.updateModel()
            
    def setPhase(self, **kwargs):
        for key in kwargs:
            if key == 'period':
                self.period = kwargs[key]
            if key == 't0':
                self.t0 = kwargs[key]
        self.phase = kep.func.foldPhase(\
            self.xdata, \
            self.t0 + 0.5*self.period - 2454900e0, \
            self.period)

    def updateModel(self):
        self.model = kep.tquick.TransitLC(\
            self.xdata,\
            1.,\
            self.inc,\
            self.aRs,\
            self.period,\
            self.RpRs,\
            self.u1,\
            self.u2,\
            self.t0)
            
    def getChiSqr(self):
        return num.sum((self.ydata - self.model)**2)


    def fitIncArsRprs(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.inc, self.aRs, self.RpRs])
            modelChiSqr = lambda args: \
                num.sum((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    args[0],\
                    args[1],\
                    self.period,\
                    args[2],\
                    self.u1,\
                    self.u2,\
                    self.t0))**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, full_output=True)
            self.inc  = output[0][0]
            self.aRs  = output[0][1]
            self.RpRs = output[0][2]
            warnflag  = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs = guess
            print 'fit one discarded'


    def fitPeriodT0(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.period, self.t0])
            modelChiSqr = lambda args: \
                num.sum((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    args[0],\
                    self.RpRs,\
                    self.u1,\
                    self.u2,\
                    args[1]))**2)
            resetPhase = lambda args: self.setPhase(\
                period=args[0], t0=args[1])
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, \
                    callback=resetPhase, full_output=True)
            self.period = output[0][0]
            self.t0     = output[0][1]
            warnflag    = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs = guess
            print 'fit two discarded'


    def fitu1(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.u1])
            modelChiSqr = lambda args: \
                num.sum((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    self.period,\
                    self.RpRs,\
                    args[0],\
                    self.u2,\
                    self.t0))**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, full_output=True)
            self.u1  = output[0][0]
            warnflag = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs = guess
            print 'u1 fit discarded'


    def fitu2(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.u2])
            modelChiSqr = lambda args: \
                num.sum((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    self.period,\
                    self.RpRs,\
                    self.u1,\
                    args[0],\
                    self.t0))**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, full_output=True)
            self.u2  = output[0][0]
            warnflag = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs = guess
            print 'u2 fit discarded'
            
            
    def finalFit(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.inc, self.aRs, self.RpRs, \
                 self.period, self.t0, \
                 self.u1, self.u2])
            modelChiSqr = lambda args: \
                num.sum((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    args[0],\
                    args[1],\
                    args[3],\
                    args[2],\
                    args[5],\
                    args[6],\
                    args[4]))**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, \
                    full_output=True)
            self.inc    = output[0][0]
            self.aRs    = output[0][1]
            self.RpRs   = output[0][2]
            self.period = output[0][3]
            self.t0     = output[0][4]
            self.u1     = output[0][5]
            self.u2     = output[0][6]
            warnflag    = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs = guess
            print 'final fit discarded'
