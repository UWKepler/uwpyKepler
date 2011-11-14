#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse

def rawplot(kid,ctype):
    lcData = kep.iodb.ReadLightCurve(kid,selection=ctype)
    #Look at difference between plot['x'] and bjdtime#
    BJDREFI = kep.iodb.getBJDREFI(kid)
    pylab.title('KID: %s' % kid)
    pylab.xlabel('BJD')
    pylab.ylabel('Flux')
    pylab.plot(lcData['x']+BJDREFI, lcData['y'], 'b.')
    pylab.show()

def pipelineplot(kid,ctype,pltlcFinal,dtopt,sharex):
    keplc = kep.keplc.keplc(kid)
    eData = keplc.eData
    BJDREFI = keplc.BJDREFI
    
    kw = kep.keplc.kw(\
    ctype=ctype,\
    gsize=2,\
    owin=15,\
    othresh=5,\
    dwin=50,\
    polyorder=6,\
    agap=1,\
    durfac=2)
    
    lcData = kep.keplc.lcData(kid,eData,BJDREFI, kw).lcData
    
    if sharex:
        ugly = pylab.subplot(211)
        pylab.title('KID: %s' % kid)
        pylab.ylabel('Flux')
        pylab.plot(lcData['x'], lcData['y'], 'b.')
        pylab.plot(lcData['x'], lcData['correction'], 'k.')
        pretty = pylab.subplot(212, sharex=ugly)
        pylab.xlabel('BJD')
        pylab.ylabel('CorrFlux')
        pylab.plot(lcData['x'], lcData['ydt'], 'b.')
        pylab.show()
        exit()
    
    i=1
    if pltlcFinal:
        pylab.figure(i)
        pylab.title('KID: %s' % kid)
        pylab.xlabel('BJD')
        pylab.ylabel('CorrFlux')
        pylab.plot(lcData['x']+BJDREFI, lcData['ydt'], 'b.')
        i+=1
    if dtopt:
        pylab.figure(i)
        pylab.title('KID: %s' % kid)
        pylab.xlabel('BJD')
        pylab.ylabel('Flux')
        pylab.plot(lcData['x']+BJDREFI, lcData['y'], 'b.')
        pylab.plot(lcData['x']+BJDREFI, lcData['correction'], 'k-')
    pylab.show()

def runpwn(kid,ctype,rawplt,pipelineplt,showDetrend,sharex,pipebool):
    if rawplt:
        rawplot(kid,ctype)
    if pipebool:
        pipelineplot(kid,ctype,pipelineplt,showDetrend,sharex)
    if not rawplt and not pipebool:
        rawplot(kid,ctype)

KID = sys.argv[1]
if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog: Use this script to make various lightcurve plots.")
    parser.add_option('-r','--rawplot',\
                        action='store_true',\
                        dest='rawplot',\
                        default=False,\
                        help='make a plot'\
                        ' showing raw lc data')
    parser.add_option('-p','--pipelineplot',\
                        action='store_true',\
                        dest='pipelineplot',\
                        default=False,\
                        help='make a plot'\
                        ' showing a detrended lc')
    parser.add_option('-d','--showDetrend',\
                        action='store_true',\
                        dest='showDetrend',\
                        default=False,\
                        help='make a plot showing'\
                        ' the detrending function')
    parser.add_option('-s','--sharex',\
                        action='store_true',\
                        dest='sharex',\
                        default=False,\
                        help='plot the detrending function'\
                        ' and detrended data with a shared x-axis')
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
    pipebool = False
    if opts.pipelineplot or opts.showDetrend or opts.sharex:
        pipebool = True
    runpwn(KID,opts.ctype,opts.rawplot,opts.pipelineplot,opts.showDetrend,opts.sharex,pipebool)