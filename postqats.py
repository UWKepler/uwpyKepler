import uwpyKepler as kep
import numpy as num
import scipy

# returns KID from QATS .data file
# cursor must be at start of file
def getKID(dfile):
    lines = dfile.readline()
    KID = int(lines.split('|')[0].split('#')[1])
    return KID

# change this with new sets of qats runs
def getdFileName(kid):
    kid = str(kid)
    sgStr = str(kep.dbinfo.getSkyGroup(kid)).zfill(3)
    dFileName = '/astro/store/student-scratch1/'\
                'johnm26/SPRING_BREAK_RUNS/SG' + sgStr + \
                '/signal.SG' + sgStr + \
                '.unflipped.' + kid + '.data'
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
def firstFareyVal(p0,p1,p2):
    m   = 1. # start on right side for when p0 = p1
    n   = 1.
    ln  = 0. # 'left numerator'
    ld  = 1. # 'left denominator'
    rn  = 1. # 'right numerator'
    rd  = 1. # 'right denominator'
    pRatio1 = p1/p0
    pRatio2 = p2/p0
    if pRatio2 > 1:
        temp = pRatio1
        pRatio1 = 1./pRatio2
        pRatio2 = 1./temp
    while m/n < pRatio1 or m/n > pRatio2:
        if m/n < pRatio1:
            ln = m
            ld = n
            m += rn
            n += rd
        else:
            rn = m
            rd = n
            m += ln
            n += ld
    return m, n

def firstFareyVals(periods, p0):
    Qs = []
    for i in range(len(periods) - 1):
        m, n = kep.postqats.firstFareyVal(p0, periods[i], periods[i + 1])
        Qs.append(1. / num.sqrt(m * n))
    Qs.append(Qs[-1])
    return num.array(Qs)

def fitQats(kid, periods, snr, polyOrder, **kwargs):
    order = polyOrder
    allCoeffs = []
    chiSqrs = []
    fits = []
    for p0 in periods:
        findQs = lambda Ps: firstFareyVals(Ps, p0)
        funcs = map(lambda n: lambda x: x**n, range(order + 1)[::-1])
        funcs.append(findQs)
        #funcs = [order, order - 1 ... 1, 0, findQs]
        # where the terms before 'findQs' denote the power
        # the data are raised to
        coeffs, fareyVals = \
        kep.linLeastSquares.linLeastSquares\
        (periods, snr, funcs, order + 2)
        allCoeffs.append(coeffs)
        polynom = scipy.polyval(coeffs[:-1], periods)
        fit = polynom + coeffs[-1] * fareyVals
        fits.append(fit)
        sqrs = (snr - fit)**2 / fit
        chiSqr = sqrs.sum()
        chiSqrs.append(chiSqr)
    minSqr = min(chiSqrs)
    minIdx = chiSqrs.index(minSqr)
    
    bestPeriod = periods[minIdx]
    bestFit = fits[minIdx]
    
    qatsBestPeriod = kep.postqats.getBestPeriodByKID(kid)
    qatsPeriodIdx = periods.tolist().index(qatsBestPeriod)
    qatsPeriodFit = fits[qatsPeriodIdx]
    qatsPeriodSqr = chiSqrs[qatsPeriodIdx]
    
    for key in kwargs:
        if key == 'qats_best_period_fit':
            if kwargs[key] == True:
                return bestPeriod, bestFit, minSqr,\
                       qatsBestPeriod, qatsPeriodFit, qatsPeriodSqr
    return bestPeriod, bestFit, minSqr
