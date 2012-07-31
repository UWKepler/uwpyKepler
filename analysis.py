import uwpyKepler as kep
import numpy as num
import pylab
from scipy import optimize
import sys
import os
import copy

eDatafileDir = '/astro/store/student-scratch1/johnm26/dFiles/'
name = 'eDataDiscoveries.txt'
name2 = 'eDataFromFits.txt'
name3 = 'binary_eDataDiscoveries.txt'

def writeEDataToFile(kid, period, t0, q, **kwargs):
    fileName = name
    insertInFile = False
    overwriteIdx = 0
    for kw in kwargs:
        # writes to binary file instead of planet file
        if kw == 'binary':
            if kwargs[kw]:
                fileName = name3
        # adds another file entry; not overwrite existing
        elif kw == 'insert':
            insertInFile = kwargs[kw]
        # specifies which file entry to overwrite if multiple
        # 0 overwrites first, 1 second, 2 third, etc.
        elif kw == 'overwrite':
            overwriteIdx = kwargs[kw]
    ofile = open(eDatafileDir + fileName, 'r')
    newline = str(kid) + '\t' \
        + str(period) + '\t' \
        + str(t0) + '\t' \
        + str(q)
    lines = ofile.readlines()
    kids = [line.split()[0] for line in lines]
    if str(kid) in kids:
        idx = kids.index(str(kid))
        if insertInFile:
            lines.insert(idx+1, newline)
        else:
            lines[idx + overwriteIdx] = newline
    else:
        lines.append(newline)
    ofile.close()
    ofile = open(eDatafileDir + fileName, 'w')
    for line in lines:
        print >> ofile, line.strip()

def geteDataFromFile(kid, **kwargs):
    fileName = name
    for kw in kwargs:
        if kw == 'binary':
            if kwargs[kw]:
                fileName = name3
    dFile = open(eDatafileDir + fileName, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            line   = line.split()
            period = float(line[1])
            t0     = float(line[2])
            q      = float(line[3])
            return period, t0, q
    return -1, -1, -1
    
def alleDataFromFile(kid, **kwargs):
    periods = []
    t0s     = []
    qs      = []
    fileName = name
    for kw in kwargs:
        if kw == 'binary':
            if kwargs[kw]:
                fileName = name3
    dFile = open(eDatafileDir + fileName, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            line   = line.split()
            periods.append( float(line[1]) )
            t0s.append( float(line[2]) )
            qs.append( float(line[3]) )
    return periods, t0s, qs
    
def eDataInFile(kid, **kwargs):
    fileName = name
    for kw in kwargs:
        if kw == 'binary':
            if kwargs[kw]:
                fileName = name3
    dFile = open(eDatafileDir + fileName, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            return True
    return False
    
def geteDataFromModelFile(kid):
    periods = []
    t0s     = []
    qs      = []
    incs    = []
    aRss    = []
    RpRss   = []
    u1s     = []
    u2s     = []
    dFile = open(eDatafileDir + name2, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(kid):
            line   = line.split()
            periods.append( float(line[1]) )
            t0s.append(     float(line[2]) )
            qs.append(      float(line[3]) )
            incs.append(    float(line[4]) )
            aRss.append(    float(line[5]) )
            RpRss.append(   float(line[6]) )
            u1s.append(     float(line[7]) )
            u2s.append(     float(line[8]) )
    params = [periods, t0s, qs, incs, aRss, RpRss, u1s, u2s]
    if len(params[0]) == 0:
        return num.zeros(8) - 1
    return params
    
class modelLC:
    def __init__(self, kid, lcData, eData, guess):
        BJDREFI = kep.iodb.getBJDREFI(kid)
        self.kid    = kid
        self.xdata  = lcData['x'] + BJDREFI - 2454900e0
        self.ydata  = lcData['ydt']
        self.yerr   = lcData['yerrdt']
        key = eData['KOI'].keys()[0]
        self.period = eData['KOI'][key]['Period']
        self.t0     = eData['KOI'][key]['T0']
        self.q      = eData['KOI'][key]['Duration']
        self.inc    = guess[0]
        self.aRs    = guess[1]
        self.RpRs   = guess[2]
        self.u1     = guess[3]
        self.u2     = guess[4]
        self.setPhase()
        self.updateModel()
    
    def copy(self):
        return copy.deepcopy(self)
            
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
        self.chiSqr = self.getChiSqr()
            
    def getChiSqr(self):
        return num.sum( ((self.ydata - self.model) / self.yerr)**2 )

    def fitIncArsRprs(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.inc, self.aRs, self.RpRs])
            modelChiSqr = lambda args: \
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    args[0],\
                    args[1],\
                    self.period,\
                    args[2],\
                    self.u1,\
                    self.u2,\
                    self.t0)) / self.yerr)**2)
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
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    args[0],\
                    self.RpRs,\
                    self.u1,\
                    self.u2,\
                    args[1])) / self.yerr)**2)
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
            self.period, self.t0 = guess
            print 'fit two discarded'


    def fitu1(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.u1])
            modelChiSqr = lambda args: \
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    self.period,\
                    self.RpRs,\
                    args[0],\
                    self.u2,\
                    self.t0)) / self.yerr)**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, full_output=True)
            self.u1  = output[0][0]
            warnflag = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.u1 = guess
            print 'u1 fit discarded'


    def fitu2(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.u2])
            modelChiSqr = lambda args: \
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    self.period,\
                    self.RpRs,\
                    self.u1,\
                    args[0],\
                    self.t0)) / self.yerr)**2)
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, full_output=True)
            self.u2  = output[0][0]
            warnflag = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.u2 = guess
            print 'u2 fit discarded'
            
            
    def finalFit(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.inc, self.aRs, self.RpRs, \
                 self.period, self.t0, \
                 self.u1, self.u2])
            modelChiSqr = lambda args: \
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    args[0],\
                    args[1],\
                    args[3],\
                    args[2],\
                    args[5],\
                    args[6],\
                    args[4])) / self.yerr)**2)
            resetPhase = lambda args: self.setPhase(\
                period=args[3], t0=args[4])
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, \
                    callback=resetPhase, full_output=True)
            self.inc    = output[0][0]
            self.aRs    = output[0][1]
            self.RpRs   = output[0][2]
            self.period = output[0][3]
            self.t0     = output[0][4]
            self.u1     = output[0][5]
            self.u2     = output[0][6]
            warnflag    = output[-1]
        self.updateModel()
        
    def fitT0(self):
        warnflag = 1
        while warnflag != 0:
            guess = num.array(\
                [self.t0])
            modelChiSqr = lambda args: \
                num.sum(((self.ydata - kep.tquick.TransitLC(\
                    self.xdata,\
                    1.,\
                    self.inc,\
                    self.aRs,\
                    self.period,\
                    self.RpRs,\
                    self.u1,\
                    self.u2,\
                    args[0])) / self.yerr)**2)
            resetPhase = lambda args: self.setPhase(\
                t0=args[0])
            oldChiSqr = self.getChiSqr()
            output = \
                optimize.fmin(modelChiSqr, guess, \
                    callback=resetPhase, full_output=True)
            self.t0     = output[0][0]
            warnflag    = output[-1]
        self.updateModel()
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.period, self.t0 = guess
            print 'fit two discarded'
        
        newChiSqr = self.getChiSqr()
        if newChiSqr > oldChiSqr:
            self.inc, self.aRs, self.RpRs,
            self.period, self.t0, \
            self.u1, self.u2 = guess
            print 'final fit discarded'
            
    def computeTransitSN(self):
        idx = num.where(self.model < 1)[0]
        N = len(idx)
        signal = num.sum(1 - self.ydata[idx])
        noise = num.sum(self.yerr[idx]**2) / num.sqrt(N)
        return signal / noise

class T0Resetter:
    def __init__(self, t0, period):
        self.t0 = t0
        self.t0_0 = t0
        self.period = period
        self.fig = pylab.gcf()
        self.axes = pylab.gca()
        # connection ID
        self.cid = \
            self.fig.canvas.mpl_connect('button_release_event', self)
    
    def __call__(self, event):
        if event.inaxes != self.axes:
            return
        self.t0 = self.t0_0 + (event.xdata - 0.5) * self.period
        print 't0 will reset to:', self.t0

def resetT0(t0, period):
    print '#-----------------------------#'
    print 'Click on center of transit'
    print 'to reset t0 to the proper value\n'
    print 'TO EXIT: close figure'
    print '#-----------------------------#'
    tReset = T0Resetter(t0, period)
    pylab.show()
    return tReset.t0

class PeriodResolver:
    def __init__(self, p0, stepsize, x, ydt):
        print '#-------------------------------------------#'
        print 'Interactive Period Resolver\n'
        print 'right arrow: increase period'
        print 'left arrow:  decrease period'
        print 'up arrow:    "zoom in"'
        print '  (decrease period stepsize by factor of 10)'
        print 'down arrow:  "zoom out"'
        print '  (increase period stepsize by factor of 10)\n'
        print 'TO EXIT: close figure'
        print '#-------------------------------------------#'
        self.phase = kep.func.foldPhase(x,0,p0)
        self.period = p0
        self.step = stepsize
        self.x = x
        self.ydt = ydt
        pylab.plot(self.phase, self.ydt, 'b.')
        self.fig = pylab.gcf()
        self.cid = \
            self.fig.canvas.mpl_connect('key_press_event', self)
        
    def __call__(self, event):
        if event.key == 'right':
            self.period += self.step
            print 'current period = ' + str(self.period)
            self.phase = kep.func.foldPhase(self.x,0,self.period)
            pylab.ion()
            pylab.cla()
            pylab.plot(self.phase, self.ydt, 'b.')
            pylab.ioff()
        elif event.key == 'left':
            self.period -= self.step
            print 'current period = ' + str(self.period)
            self.phase = kep.func.foldPhase(self.x,0,self.period)
            pylab.ion()
            pylab.cla()
            pylab.plot(self.phase, self.ydt, 'b.')
            pylab.ioff()
        elif event.key == 'up':
            self.step /= 10.
        elif event.key == 'down':
            self.step *= 10.