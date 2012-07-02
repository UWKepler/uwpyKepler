#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse
import os

KID = sys.argv[1]
if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "Use this script to find a more precise period.\nScript produces 10 plots per run.\nScript takes a KID as a system argument.")
    parser.add_option('-m','--manual',\
                        type=float,\
                        dest='manual',\
                        default=None,\
                        help='input a manual period')
    parser.add_option('-z','--zoomfactor',\
                        type=int,\
                        dest='zoomfac',\
                        default=1,\
                        help='input a zoomfactor')
    parser.add_option('-p','--plot',\
                        action='store_true',\
                        dest='plot',\
                        default=False,\
                        help='plots rather than saves figures')
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

kw = kep.keplc.kw(\
ctype=opts.ctype,\
gsize=2,\
owin=15,\
othresh=5,\
dwin=50,\
polyorder=6,\
agap=1,\
durfac=2)

lcData = kep.keplc.lcData(KID,eData,BJDREFI,kw).lcData

if opts.manual:
    p0 = opts.manual
else:
    try:
        eData_idx = eData['KOI'].keys()[0]
        eData = eData['KOI'][eData_idx]
        p0 = eData['Period']
        print 'catalogued period:', eData['Period']
    except:
        try:
            p0 = kep.postqats.getBestPeriodByKID(KID)
            print 'QATS best period:', p0
        except:
            print('no period data found; exiting...\n')
            sys.exit()

zf = opts.zoomfac
t_dur0 = kep.qats.tdur(1.4,0,p0)
t_dur  = t_dur0/(10.**(zf-1))

testPs = num.arange(p0-t_dur,p0+t_dur,2*t_dur/10.)
for plot in range(len(testPs)):
    phase = kep.func.foldPhase(lcData['x'],0,testPs[plot])
    pylab.figure(plot)
    pylab.title('KID: %s, P0: %s, ZoomFac: %s' % \
                (KID, testPs[plot], zf))
    pylab.xlabel('Phase')
    pylab.ylabel('Flux')
    pylab.plot(phase, lcData['ydt'], 'b.')
    if not opts.plot:
        pylab.savefig('%s_%s_%s.png' %\
                      (KID, testPs[plot], zf))
if opts.plot:
    pylab.show()


