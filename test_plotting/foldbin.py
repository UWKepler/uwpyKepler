#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse
import os

def getphase(eData,kid,auto,man):
    if not man:
        try:
            eKeys = eData['KOI'].keys()
            if len(eKeys) > 1:
                print 'multiple catalogued periods:'
                for key in eKeys:
                    print '\t' + key + ' period: ' + str(eData['KOI'][key]['Period'])
                period = eData['KOI'][eKeys[0]]['Period']
                t0 = eData['KOI'][eKeys[0]]['T0']
                print 'using period for ' + eKeys[0] + ': ' + str(period)
            else:
                eData_idx = eKeys[0]
                eData = eData['KOI'][eData_idx]
                period = eData['Period']
                t0 = eData['T0']
                print 'catalogued period:', period
            phase = \
            kep.func.foldPhase(lcData['x'],t0 + 0.5*period,period)
        except:
            try:
                period = kep.postqats.getBestPeriodByKID(kid)
                phase = kep.func.foldPhase(lcData['x'],0,period)
                print 'QATS best period:', period
            except:
                print('no period data found; exiting...\n')
                sys.exit()
    else:
        phase = kep.func.foldPhase(lcData['x'],0,man)

    return phase

def bin(binsize,phase):
    x = num.array(phase)
    binnedx    = []
    binnedy    = []
    binnedyerr = []
    for bin in num.arange(0,1,binsize):
        idx   = num.where((x >= bin) \
              & (x < bin+binsize))[0]
        xbin  = num.mean(x[idx])
        ybin  = num.mean(lcData['ydt'][idx])
        ybstd = num.std(lcData['ydt'][idx])
        binnedx.append(xbin)
        binnedy.append(ybin)
        binnedyerr.append(ybstd)
        
    return binnedx, binnedy, binnedyerr

KID = sys.argv[1]
if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog\nUse this script to fold or bin a given lightcurve.\nScript takes a KID as a system argument.")
    parser.add_option('-a','--auto',\
                        action='store_true',\
                        dest='auto',\
                        default=True,\
                        help='query the database and qats file '\
                        'for a known period')
    parser.add_option('-m','--manual',\
                        type=float,\
                        dest='manual',\
                        default=None,\
                        help='input a manual period')
    parser.add_option('-b','--bin',\
                        type=float,\
                        dest='bin',\
                        default=False,\
                        help='bin the lightcurve according '\
                        'to a specified binsize')
    parser.add_option('-e','--errbars',\
                        action='store_true',\
                        dest='errbars',\
                        default=False,\
                        help='if -b was chosen, '\
                        'bin with errorbars')
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

keplc = kep.keplc.keplc(KID)
eData = keplc.eData
BJDREFI = keplc.BJDREFI

kw = kep.quicklc.quickKW(ctype=opts.ctype)

lcData = kep.keplc.lcData(KID,eData,BJDREFI,kw).lcData
lcData['x'] = lcData['x'] + BJDREFI

phase = getphase(eData,KID,opts.auto,opts.manual)

if opts.bin:
    bx, by, byerr = bin(opts.bin,phase)
    if not opts.errbars:
        byerr = None
    pylab.errorbar(bx, by, yerr=byerr, fmt='.')
    pylab.title('Binsize: %s' % opts.bin)
    pylab.show()
    sys.exit()

pylab.title('KID: %s' % KID)
pylab.xlabel('Phase')
pylab.ylabel('Flux')
pylab.plot(phase, lcData['ydt'], 'b.')
pylab.show()
