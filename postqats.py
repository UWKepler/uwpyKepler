import uwpyKepler as kep
import numpy as num
import scipy
import pylab
import os

# returns KID from QATS .data file
# cursor must be at start of file
def getKID(dfile):
    lines = dfile.readline()
    KID = int(lines.split('|')[0].split('#')[1])
    return KID

# change this with new sets of qats runs
def getdFileName(kid, quiet = True):
    kid = str(int(kid))
    sgStr = str(kep.dbinfo.getSkyGroup(kid)).zfill(3)
    basedir = '/astro/store/student-scratch1/johnm26'
    primary_path = os.path.join(basedir, 'SPRING_BREAK_RUNS/SG' + sgStr)
    secondary_path = os.path.join(basedir, 'KOI_RUNS/qatsKOIdata')
    fname = 'signal.SG' + sgStr + '.unflipped.' + kid + '.data'
    dFileName = os.path.join(primary_path, fname)
    if not os.path.isfile(dFileName):
        dFileName = os.path.join(secondary_path, fname)
    
    # quiet defaults to False so as not to interfere with condor
    if not quiet:
        if not os.path.isfile(dFileName):
            print 'WARNING: ' + kid + 'qats data file not found'

    return dFileName

def getBestPeriod(dfile):
    line = dfile.readline()
    return float(line.split('|')[1].split(' ')[1])

def getBestPeriodByKID(kid):
    dfile = open(getdFileName(kid), 'r')
    return getBestPeriod(dfile)

def getQatsData(dfile):
    lines = dfile.readlines()
    if lines[0].startswith('#'):
        del lines[0]
    
    
    #Defining lists to be populated with data
    periods = []
    snrLC=[]
    snrFLAT=[]
    
    #Populating lists with data from file
    for line in lines:
        row = line.split(' | ')
        periods.append(float(row[0]))
        snrLC.append(float(row[1]))
        snrFLAT.append(float(row[2]))

    periods = num.array(periods)
    snrLC = num.array(snrLC)
    snrFLAT = num.array(snrFLAT)
    
    # get actual Signal to Noise Ratio
    snr = snrLC/snrFLAT
    
    return periods, snr, snrLC, snrFLAT
        
# used in flagQats fitting
# note: m=1 and n=1 when p0=p1 or p0=p2
def firstFareyVals(p0,p1,p2):
    #print p1
    #print p2

    m   = num.ones(len(p1))  # start on right side for when p0 = p1
    n   = num.ones(len(p1))
    ln  = num.zeros(len(p1)) # 'left numerator'
    ld  = num.ones(len(p1))  # 'left denominator'
    rn  = num.ones(len(p1))  # 'right numerator'
    rd  = num.ones(len(p1))  # 'right denominator'
    pRatio1 = p1/p0
    pRatio2 = p2/p0

    idx = num.where(pRatio2 > 1)
    temp = pRatio1[idx]
    pRatio1[idx] = 1./pRatio2[idx]
    pRatio2[idx] = 1./temp

    idx = num.where((m/n > pRatio2) | (m/n < pRatio1))[0]
    while len(idx) != 0:
        idxG = num.where(m/n > pRatio2)[0]
        rn[idxG] = m[idxG]
        rd[idxG] = n[idxG]
        m[idxG] += ln[idxG]
        n[idxG] += ld[idxG]

        idxL = num.where(m/n < pRatio1)[0]
        ln[idxL] = m[idxL]
        ld[idxL] = n[idxL]
        m[idxL] += rn[idxL]
        n[idxL] += rd[idxL]

        idx = num.where((m/n > pRatio2) | (m/n < pRatio1))[0]

    return m, n

def sqrtMN(p0, Ps):
    m, n = firstFareyVals(p0, Ps[:-1], Ps[1:])
    return num.sqrt(m * n)

def fitQats(kid, periods, snr, polyOrder, **kwargs):
    order = polyOrder
    allCoeffs = []
    chiSqrs = []
    fits = []
    # periods and snr set to indices [:-1] b/c farey sequence
    # only returns values for periods[:-1]
    snr_cut = snr[:-1]
    periods_cut = periods[:-1]
    for p0 in periods:
        findQs = lambda Ps: 1. / sqrtMN(p0, periods)
        # funcs = [order, order - 1 ... 1, 0, findQs]
        # where the terms before 'findQs' denote the power
        # the data are raised to
        funcs = map(lambda n: lambda x: x[:-1]**n, range(order + 1)[::-1])
        funcs.append(findQs)
        coeffs, fareyVals = \
        kep.linLeastSquares.linLeastSquares\
        (periods, snr_cut, funcs, order + 2, return_func_vals=-1)
        allCoeffs.append(coeffs)
        polynom = scipy.polyval(coeffs[:-1], periods_cut)
        fit = polynom + coeffs[-1] * fareyVals
        fits.append(fit)
        sqrs = (snr_cut - fit)**2 / fit
        chiSqr = sqrs.sum()
        chiSqrs.append(chiSqr)
    minSqr = min(chiSqrs)
    minIdx = chiSqrs.index(minSqr)
    
    bestPeriod = periods[minIdx]
    bestFit = fits[minIdx]
    
    returnList = [(bestPeriod, bestFit, minSqr)]
    for key in kwargs:
        if key == 'qats_best_period_fit':
            if kwargs[key]:
                qatsBestPeriod = getBestPeriodByKID(kid)
                qatsPeriodIdx = periods.tolist().index(qatsBestPeriod)
                qatsPeriodFit = fits[qatsPeriodIdx]
                qatsPeriodSqr = chiSqrs[qatsPeriodIdx]
                returnList.insert(1, (qatsBestPeriod, qatsPeriodFit, qatsPeriodSqr))
        elif key == 'return_coeffs':
            if kwargs[key]:
                returnList.append(tuple(allCoeffs[minIdx]))
        elif key == 'full_output':
            if kwargs[key]:
                out = {'p0s':num.array(periods), \
                'chiSqrs':num.array(chiSqrs), \
                'coeffs':num.array(allCoeffs)}
                return out

    return tuple(returnList)

def getFitFileName(kid, **kwargs):
    path = './'
    fname = str(kid) + 'QatsFit.txt'
    for kw in kwargs:
        if kw == 'path':
            path = kwargs[kw]
        elif kw == 'fname':
            fname = kwargs[kw]
    fname = os.path.join(path, fname)
    return fname

def chiSqr(data, model):
    return num.sum((data - model)**2 / model)

def maxIndices(a):
    aRight = num.roll(a, -1)
    aLeft  = num.roll(a, 1)
    maxIndices = num.where((a > aRight) & (a > aLeft))[0]
    return maxIndices

def getNthMaxIndex(a, n):
    imaxes = maxIndices(a)
    maxes = a[imaxes]
    idx = maxes.argsort()
    return num.where(a == maxes[idx[-1 * n]])[0]

class QatsFeaturesModel:
    #####
    # initialization methods #
    #####
    # pass kid, periods, snr at initialization
    # or use method: setQatsData()
    # or use method: fromFile()
    def __init__(self, *args):
        data = [None, None, None]
        for i in range(len(args)):
            data[i] = args[i]
        self.kid, self.periods, self.snr = data
    
    def setQatsData(self, kid, periods, snr):
        self.kid      = kid
        self.periods  = periods
        self.snr      = snr
    
    #####
    # essential fit methods #
    #####
    def fitQats(self, polyOrder):
        dict = fitQats(self.kid, self.periods, self.snr, \
            polyOrder, full_output=True)
        self.chiSqrs = dict['chiSqrs']
        self.coeffs  = dict['coeffs']
        self.setBestFit()
        
    def setBestFit(self):
        imin = self.chiSqrs.argmin()
        self.polynom = scipy.polyval(self.coeffs[imin][:-1], \
            self.periods[:-1])
        self.bestFit = self.polynom + self.coeffs[imin][-1] * \
            1. / sqrtMN(self.periods[imin], self.periods)
        self.minSqr  = self.chiSqrs[imin]
    
    #####
    # central features method #
    #####
    def returnFeatures(self, iters = 1, **kwargs):
        self.setConvolved()
        #cmax = self.getCMax()
        self.setD_chiSqr()
        imin = self.chiSqrs.argmin()
        #coeffs = self.getPolynomCoeffs()[imin]
        coeffs = self.getCoeffs()[imin]
        # iters defines how many peaks the model will examine
        # for each feature
        n_cmax       = iters
        n_snrmax     = iters
        n_peakwidth  = iters
        n_dchi       = iters
        n_peakperiod = iters
        for kw in kwargs:
            if kw == 'n_cmax':
                n_cmax = kwargs[kw]
            elif kw == 'n_snrmax':
                n_snrmax = kwargs[kw]
            elif kw == 'n_peakwidth':
                n_peakwidth = kwargs[kw]
            elif kw == 'n_dchi':
                n_dchi = kwargs[kw]
            elif kw == 'n_peakperiod':
                n_peakperiod = kwargs[kw]
        cmaxes      = num.zeros(n_cmax)
        snrMaxes    = num.zeros(n_snrmax)
        peakWidths  = num.zeros(n_peakwidth)
        dchis       = num.zeros(n_dchi)
        peakPeriods = num.zeros(n_peakperiod)
        for i in num.arange(n_cmax) + 1:
            cmaxes[i - 1] = self.getCMax(which_peak=i)
        for i in num.arange(n_snrmax) + 1:
            snrMaxes[i - 1] = self.getSnrMax(which_peak=i)
        for i in num.arange(n_peakwidth) + 1:
            peakWidths[i - 1] = self.getPeakWidth(which_peak=i)
        for i in num.arange(n_dchi) + 1:
            dchis[i - 1] = self.getMaxD_chiSqr(which_peak=i)
        for i in num.arange(n_peakperiod) + 1:
            peakPeriods[i - 1] = self.getPeakPeriod(which_peak=i)
        
        features = num.array([])
        features = num.hstack((features, cmaxes))
        #print features
        features = num.hstack((features, coeffs))
        #print features
        features = num.hstack((features, snrMaxes))
        #print features
        features = num.hstack((features, peakWidths))
        #print features
        features = num.hstack((features, dchis))
        #print features
        #exit()
        #features = num.hstack((features, peakPeriods))
        return features
    
    #####
    # 'SET' features methods #
    #####
    def setConvolved(self):
        imin = self.chiSqrs.argmin()
        baseline = self.polynom
        dataPeaks = self.snr[:-1] - baseline
        modelPeaks = self.bestFit - baseline
        #cmax = max(dataPeaks * modelPeaks)
        self.convolved = dataPeaks * modelPeaks
        #self.cmax = max(self.convolved)
    
    def setD_chiSqr(self):
        baseline_chiSqr = chiSqr(self.snr[:-1], self.polynom)
        self.d_chiSqr = baseline_chiSqr - self.chiSqrs
    
    #####
    # 'GET' features methods #
    #####
    def getSnrMax(self, **kwargs):
        # wp = 'which peak';
        # an integer defining which peak is chosen
        # 1 = first highest, 2 = second highest, etc...
        wp = 1
        for kw in kwargs:
            if kw == 'which_peak':
                wp = kwargs[kw]
        dataPeaks = self.snr[:-1] - self.polynom
        imaxes = maxIndices(dataPeaks)
        maxes = dataPeaks[imaxes]
        maxes.sort()
        return maxes[-1 * wp]
        #dataPeaks.sort()
        #return dataPeaks[maxes[-1 * wp]]
    
    def getCMax(self, **kwargs):
        wp = 1
        for kw in kwargs:
            if kw == 'which_peak':
                wp = kwargs[kw]
        imaxes = maxIndices(self.convolved)
        maxes = self.convolved[imaxes]
        maxes.sort()
        return maxes[-1 * wp]
    
    def getPeakWidth(self, **kwargs):
        wp = 1
        for kw in kwargs:
            if kw == 'which_peak':
                wp = kwargs[kw]
        snr = self.snr[:-1]
        dataPeaks = snr - self.polynom
        imax = getNthMaxIndex(dataPeaks, wp)
        snrRight = num.roll(snr, -1)
        snrLeft  = num.roll(snr, 1)
        imins = num.where((snr < snrRight) & (snr < snrLeft) \
            & (dataPeaks < 0.5 * dataPeaks[imax]))[0]
        nearestMin = num.argmin(num.abs(imins - imax))
        try:
            if imins[nearestMin] > imax:
                width = imins[nearestMin] - imins[nearestMin - 1]
            else:
                width = imins[nearestMin + 1] - imins[nearestMin]
        # for edge cases; rough guess hopefully won't 
        # disrupt machine learning categorization
        except:
            width = 4

        return width
    
    def getPeakPeriod(self, **kwargs):
        wp = 1
        for kw in kwargs:
            if kw == 'which_peak':
                wp = kwargs[kw]
        dataPeaks = self.snr[:-1] - self.polynom
        imax = getNthMaxIndex(dataPeaks, wp)
        return self.periods[imax][0]
    
    def getMaxD_chiSqr(self, **kwargs):
        wp = 1
        for kw in kwargs:
            if kw == 'which_peak':
                wp = kwargs[kw]
        imax = getNthMaxIndex(self.d_chiSqr, wp)
        # sometimes the best chi square improvement
        # does not occur at the snr peak period(s)
        # the following lines ensure chi square improvements
        # from these periods are chosen
        #dataPeaks = self.snr - self.polynom
        #imax = getNthMaxIndex(dataPeaks, wp)
        return self.d_chiSqr[imax]
    
    def getPolynomCoeffs(self):
        return self.coeffs[:, :-1]
    
    def getCoeffs(self):
        return self.coeffs
    
    #####
    # file write/read methods #
    #####
    def toFile(self, **kwargs):
        fname = getFitFileName(self.kid, **kwargs)
        ofile = open(fname, 'w')
        # column 1: periods
        periods = self.periods.reshape(len(self.periods), -1)
        # column 2: snr
        snr     = self.snr.reshape(len(self.snr), -1)
        # column 3: chiSqrs
        chiSqrs = self.chiSqrs.reshape(len(self.chiSqrs), -1)
        # column 4+: coeffs
        out = num.hstack((periods, snr, chiSqrs, self.coeffs))
        num.savetxt(ofile, out)
        ofile.close()
    
    def fromFile(self, **kwargs):
        fname = getFitFileName(self.kid, **kwargs)
        ifile = open(fname, 'r')
        indata = num.loadtxt(ifile)
        self.periods = indata[:, 0]
        self.snr     = indata[:, 1]
        self.chiSqrs = indata[:, 2]
        self.coeffs  = indata[:, 3:]
        self.setBestFit()
        ifile.close()
    
    #####
    # plot methods #
    #####
    def plot(self):
        pylab.figure()
        ax = pylab.gca()
        ax.semilogx()
        imin = num.argmin(self.chiSqrs)
        qatsBestPeriod = getBestPeriodByKID(self.kid)
        iqats = num.where(self.periods == qatsBestPeriod)[0]
        ax.plot(self.periods, self.snr, 'k-',\
        self.periods[:-1], self.bestFit, 'c-',\
        self.periods[imin], self.snr[imin], 'ro',\
        self.periods[iqats], self.snr[iqats], 'bo',)
        ax.legend(('qats SNR', 'fit to qats', \
                'fit best period', 'qats best period'))
        pylab.title('KID: %s' % self.kid)
        pylab.ylabel('SNR')
        pylab.xlabel('periods')
        pylab.show()
    
    def plot_convolution(self):
        pylab.figure()
        ax = pylab.gca()
        ax.semilogx()
        pylab.title('KID: %s' % self.kid)
        pylab.ylabel('convolved')
        pylab.xlabel('periods')
        imax = self.convolved.argmax()
        pylab.plot(self.periods[:-1], self.convolved, 'k-')
        pylab.plot(self.periods[imax], self.convolved[imax], 'bo')
        pylab.show()
    
    def plot_snrMaxes(self):
        pylab.figure()
        ax = pylab.gca()
        ax.semilogx()
        pylab.title('KID: %s' % self.kid)
        pylab.ylabel('SNR')
        pylab.xlabel('periods')
        idx = []
        for n in num.arange(3) + 1:
            m = self.getSnrMax(which_peak=n)
            idx.append(num.where(self.snr[:-1] - self.polynom == m)[0])
        idx = num.array(idx)
        pylab.plot(self.periods, self.snr, 'k-')
        pylab.plot(self.periods[idx], self.snr[idx], 'bo')
        pylab.show()
        
    def print_peakWidth(self):
        for n in num.arange(3) + 1:
            width = self.getPeakWidth(which_peak=n)
            height = self.getSnrMax(which_peak=n)
            period = self.getPeakPeriod(which_peak=n)
            print 'Peak ' + str(n) + ' width  = ' + str(width)
            print 'Peak ' + str(n) + ' height = ' + str(height)
            print 'Peak ' + str(n) + ' period = ' + str(period)
    
    def plot_dchiMaxes(self):
        pylab.figure()
        ax = pylab.gca()
        ax.semilogx()
        pylab.title('KID: %s' % self.kid)
        pylab.ylabel('delta chi')
        pylab.xlabel('periods')
        idx = []
        for n in num.arange(3) + 1:
            m = self.getMaxD_chiSqr(which_peak=n)
            idx.append(num.where(self.d_chiSqr == m)[0])
        idx = num.array(idx)
        pylab.plot(self.periods, self.d_chiSqr, 'k-')
        pylab.plot(self.periods[idx], self.d_chiSqr[idx], 'bo')
        pylab.show()
    
    