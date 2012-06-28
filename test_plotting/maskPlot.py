#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse
import os
import pdb

def geteData(kid):
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

def sharexPlot(lc):
    lcData = lc.lcFinal
    idx = num.where(lcData['eMask'] == 1)[0]
    tx = lcData['x'][idx]
    ty = lcData['y'][idx]
    tydt = lcData['ydt'][idx]
    ax1 = pylab.subplot(211)
    pylab.title('KID: %s' % kid)
    pylab.ylabel('Flux')
    pylab.plot(lcData['x'], lcData['y'], 'b.')
    pylab.plot(lcData['x'], lcData['correction'], 'k.')
    pylab.plot(tx, ty, 'ro')
    ax2 = pylab.subplot(212, sharex=ax1)
    pylab.xlabel('BJD')
    pylab.ylabel('CorrFlux')
    pylab.plot(lcData['x'], lcData['ydt'], 'b.')
    pylab.plot(tx, tydt, 'ro')
    
def foldPlot(lc):
    period = lc.eData['KOI']['fake']['Period']
    t0 = lc.eData['KOI']['fake']['T0']
    lcData = lc.lcFinal
    phase = kep.func.foldPhase(lcData, t0 + 0.5 * period, period)
    pylab.title('KID: %s' % kid)
    pylab.ylabel('CorrFlux')
    pylab.xlabel('Phase')
    idx = num.where(lcData['eMask'] == 1)[0]
    tphase = phase[idx]
    ty = lcData['ydt'][idx]
    pylab.plot(phase, lcData['ydt'], 'b.')
    pylab.plot(tphase, ty, 'ro')

kid = sys.argv[1]
if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog\nUse this script to make various lightcurve plots,\nmasking transits using eclipse data from file:\n\t/astro/store/student-scratch1/johnm26/dFiles/eDataDiscoveries.txt\nScript requires a Kepler ID as a system argument.")
    parser.add_option('-f','--fold',\
                        action='store_true',\
                        dest='fold',\
                        default=False,\
                        help='fold the lightcurve')
    cChoice = ('LC','SC','')
    parser.add_option('-c','--ctype',\
                        choices=cChoice,\
                        type='choice',\
                        dest='ctype',\
                        default='LC',\
                        help='The expected cadence type('\
                        +', '.join(cChoice)\
                        +') [default: %default]')

    opts, args = parser.parse_args()
    period, t0, q = geteData(kid)
    lc = kep.keplc.keplc(kid)
    lc.eData['eDataExists'] = True
    lc.eData['KOI'] = {'fake':\
                      {'Period':period, 'T0':t0, 'Duration':q}}
    kw = kep.quicklc.quickKW(ctype=opts.ctype)
    lc.runPipeline(kw)
    lc.lcFinal['x'] = lc.lcFinal['x'] + lc.BJDREFI
    
    if opts.fold:
        foldPlot(lc)
    else:
        sharexPlot(lc)

    pylab.show()