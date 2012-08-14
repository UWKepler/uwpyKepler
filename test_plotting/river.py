#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import numpy as num
import pylab
import sys
#import copy
import optparse
import pdb
#from mpl_toolkits.axes_grid1.axes_rgb import make_rgb_axes, RGBAxes

# read data,
#   need lc, eData from fit
#   optional argument to select planet "number" to analyze
#   if argument not given, prompt user based on data in file
# make lightcurve
#   remove other transit data points

kid = sys.argv[1]
FRACTION = 0.8

class KOI:
    def __init__(self, kid, period, t0, q):
        self.kid    = kid
        self.period = period
        self.t0     = t0
        self.q      = q
        
    def toString(self):
        return 'period=%s, t0=%s, q=%s' % \
            (self.period, self.t0, self.q)

# prompts user for which eData to use
def getPlanets(ed):
    count = 0
    koiTups = []
    print 'select an object to analyze:'
    for key in ed.keys():
        if len(ed[key]) > 0:
            print 'from ' + key
        for koi in ed[key]:
            count += 1
            print '\t' + str(count) + ': ' + koi.toString()
            koiTups.append((count, koi, key))
    prompt = 'input '
    choiceList = num.arange(count) + 1
    for choice in choiceList:
        if choice != choiceList[-1]:
            prompt += str(choice) + ' '
        else:
            prompt += 'or ' + str(choice)
    objChoice = int(raw_input(prompt + '\n> '))
    maskThese = []
    for koi in koiTups:
        if koi[0] == objChoice:
            choice = koi[1]
        elif koi[2] == 'database' \
        or koi[2] == 'transitFit file':
            maskThese.append(koi[1])

    return choice, maskThese

# retrieves all archived eData
def geteDataDict(lc):
    eds = {'database':[], 'eDataDiscovery file':[], \
        'transitFit file': []}
    eDataDiscovCount = 0
    transitFitFileCount = 0
    eData1 = lc.eData
    if eData1['eDataExists']:
        for key in eData1['KOI'].keys():
            p = eData1['KOI'][key]['Period']
            t = eData1['KOI'][key]['T0']
            q = eData1['KOI'][key]['Duration']
            koi = KOI(kid, p, t, q)
            eds['database'].append(koi)
    eData2 = kep.analysis.alleDataFromFile(lc.KID)
    for i in range(len(eData2[0])):
        koi = KOI(kid, eData2[0][i], eData2[1][i], eData2[2][i])
        eds['eDataDiscovery file'].append(koi)
    eData3 = kep.analysis.geteDataFromModelFile(lc.KID)
    for i in range(len(eData3[0])):
        koi = KOI(kid, eData3[0][i], eData3[1][i], eData3[2][i])
        eds['transitFit file'].append(koi)

    return eds

# mask transits that would otherwise disrupt plot color scheme
# by their depth
def removeOtherTransits(lc, maskThese):
    for i in range(len(maskThese)):
        lc.eData['eDataExists'] = True
        lc.eData['KOI']['fake' + str(i + 1)] = { \
            'Period':maskThese[i].period, \
            'T0':maskThese[i].t0, \
            'Duration':maskThese[i].q}
    kw = kep.quicklc.quickKW()
    lc.runPipeline(kw)
    kep.lcmod.ApplyMask(lc.lcFinal, 'eMask')

# a is "anchor" (center point of reshaping data)
# returns x and ydt reshaped according to phase cycles
def cycleData(lc, a, p):
    xdata = num.array(lc.lcFinal['x']) - a
    ydata = lc.lcFinal['ydt']
    nCycles = num.floor((xdata.max() - xdata.min()) / p)
    xdata = xdata / p - xdata // p

    diffs = xdata - num.roll(xdata, 1)
    cycles = num.where(diffs < 0)[0]
    cycles = num.hstack( (cycles, -1) )
    
    lightCycles = []

    for i in range(len(cycles) - 1):
        xCycle = xdata[cycles[i]:cycles[i + 1]]
        yCycle = ydata[cycles[i]:cycles[i + 1]]
        lightCycles.append((xCycle, yCycle))
    
    return lightCycles

# determines binsize for refine shape:
# takes the mean of all phase cycle lengths 
# after rejecting 1-sigma outliers
def findNBoxes(data_phaseCycled):
    lengths = num.array([len(data_phaseCycled[i][0]) \
        for i in range(len(data_phaseCycled))])
    idx = num.where(num.abs(lengths-num.mean(lengths)) < num.std(lengths))[0]
    lengths = lengths[idx]
    return num.round(num.mean(lengths))

# bins data_phasedCycled so that each array has the same length
def refineShape(phasedData, nBoxes):
    # binned phase arrays
    bpx = num.zeros([len(phasedData), nBoxes])
    bpy = num.zeros([len(phasedData), nBoxes])
    for i in range(len(phasedData)):
        bx, by = kep.binningdetect.binTime_nBins( \
            phasedData[i][0], phasedData[i][1], nBoxes)
        if len(bx) != len(bpx[i]):
            bx = num.hstack( (num.zeros(len(bpx[i]) - len(bx)), bx) )
            by = num.hstack( (num.zeros(len(bpy[i]) - len(by)), by) )
        bpx[i] = bx
        bpy[i] = by
    bpx = num.ma.MaskedArray(bpx, num.isnan(bpx))
    bpy = num.ma.MaskedArray(bpy, num.isnan(bpy))
    idx1 = num.where(bpx == 0)
    idx2 = num.where(bpy == 0)
    bpx.mask[idx1] = True
    bpy.mask[idx2] = True
    
    return bpx, bpy

# uses binned data to produce the "river" image array
def getImageArray(x, y):
    # determine horizontal length of image array
    phaseDiffs = num.abs(x[1:] - x[:-1])
    phaseDiffs = phaseDiffs.ravel()
    sortIdx = phaseDiffs.argsort()
    dx = phaseDiffs[sortIdx[num.floor(0.1585*len(phaseDiffs))]]
    db = FRACTION*dx
    xrng = x.max() - x.min()
    L = num.ceil(xrng / db)
    print 'len(x[0]) , len(image[0]):', len(x[0]),',',L
    initVals = num.zeros( (len(x), L) )
    image = num.ma.MaskedArray(initVals, initVals) # second initVals is the mask
    # get spacing between adjacent x to reject bad xvals in nested loop
    dx2 = num.median(x[:,1:] - x[:,:-1])
    reject_thresh = 1.5*dx2
    # fill image array with y values
    # values are assigned based on how each horizontal index
    #   corresponds to a given x value
    for i in range(len(image)):
        for j in range(len(image[i])):
            bestXVal = xrng * (float(j) / L)
            bestYIdx = (num.abs(x[i] - bestXVal)).argmin()
            if num.abs(x[i][bestYIdx] - bestXVal) > reject_thresh:
                image.mask[i][j] = True
            else:
                image[i][j] = y[i][bestYIdx]

    return image

# used to convert x-axis labels from array indices to phase units
def imageIndexToTimestamp(x, image, idx):
    xrng = x.max() - x.min()
    L = len(image[0])
    bestXVal = xrng * (float(idx) / L)
    closestXIdx = (num.abs(x - bestXVal)).argmin()
    return x[closestXIdx]

# converts x-axis labels from array indices to phase units
def setTimeTicks(ax, x, image):
    ticks = ax.get_xticks()
    #argx = num.round(0.5 * len(x))
    maskCounts = num.zeros(len(x))
    for i in range(len(x)):
        maskCounts[i] = len(num.where(x[i].mask)[0])
    argx = maskCounts.argmin()
    getTimestamps = \
        lambda tick: imageIndexToTimestamp(x[argx], image, tick)
    newticks = num.array(map(getTimestamps, list(ticks)))
    newticks = num.round(newticks, decimals=3)
    ax.set_xticklabels(newticks)

# main method: generates plot and specifies visual parameters
def river(lc, koi):
    p = koi.period
    t = koi.t0
    q = koi.q
    t += 2454900e0 - lc.BJDREFI[0]
    data_phaseCycled = cycleData(lc, t + 0.5 * p, p)
    nBoxes = findNBoxes(data_phaseCycled)
    x, y = refineShape(data_phaseCycled, nBoxes)
    image = getImageArray(x, y)
    #left = num.floor(len(image[0])*0.49)
    #right = num.ceil(len(image[0])*0.60)
    left = num.floor( len(image[0]) * (0.5 - (q/24.)/p) )
    right = num.ceil( len(image[0]) * (0.5 + (q/24.)/p) )
    flat = image.ravel()
    flat = flat[num.where(flat.mask == False)[0]]
    vmin = flat[flat.argsort()[num.floor(0.05*len(flat))]]
    vmax = flat[flat.argsort()[num.ceil(0.95*len(flat))]]
    pylab.imshow(image, aspect='auto', \
        cmap=pylab.matplotlib.cm.GnBu_r, interpolation='nearest', vmin=vmin, vmax=vmax)
    #pylab.imshow(image, aspect='auto', \
        #cmap=pylab.matplotlib.cm.Accent, interpolation='nearest', vmin=vmin)
    #pylab.imshow(image, aspect='auto', \
        #cmap=pylab.matplotlib.cm.Accent, interpolation='nearest')
    ax = pylab.gca()
    ax.set_xlim((left, right))
    setTimeTicks(ax, x, image)
    pylab.title('River: %s' % lc.KID)
    pylab.xlabel('phase')
    pylab.ylabel('transit number')
    pylab.colorbar()
        

if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog\nUse this script to make a 'river' plot (visual representation of TTVs).\nScript requires a Kepler ID as a system argument.")
    parser.add_option('-n','--planetNumber',\
                        type=int,\
                        dest='pn',\
                        default=None,\
                        help='specify which planet under KID'\
                        ' in file to use;\n0 is first, ' \
                        '1 is second, etc.')
    opts, args = parser.parse_args()

    lc = kep.keplc.keplc(kid)
    eds = geteDataDict(lc)
    koi, toMask = getPlanets(eds)
    removeOtherTransits(lc, toMask)
    river(lc, koi)

    pylab.show()


