#!/astro/apps/pkg/python64/bin//python

import sys
import uwpyKepler as kep
import numpy as num
import pylab
import optparse

def runbin(kid,ctype,share,one,man):
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
    
    lcData = kep.keplc.lcData(kid,eData,BJDREFI,kw).lcData
    lcData['x'] = lcData['x'] + BJDREFI
    
    if not man:
        lst = [1, 5, 10, 20, 50, 100]
        bind = kep.binningdetect.bin(lcData, lst)
        colors = ['r','y','b','c','g','m','k']
        if share:
            nplots = len(lst)
            for binidx in range(len(bind.keys())):
                if binidx == 0:
                    ax1 = pylab.subplot(nplots,1,binidx+1)
                    pylab.plot(bind[bind.keys()[binidx]]['x'],\
                    bind[bind.keys()[binidx]]['y'],colors[binidx]+'.')
                else:
                    pylab.subplot(nplots,1,binidx+1,sharex=ax1)
                    pylab.plot(bind[bind.keys()[binidx]]['x'],\
                    bind[bind.keys()[binidx]]['y'],colors[binidx]+'.')
                pylab.title('Binsize: %s' % lst[binidx])
        elif one:
            for binidx in range(len(bind.keys())):
                pylab.plot(bind[bind.keys()[binidx]]['x'],\
                bind[bind.keys()[binidx]]['y'],colors[binidx]+'o')
    else:
        lst = [man]
        bind = kep.binningdetect.bin(lcData, lst)
        pylab.plot(bind['binsize'+str(man)]['x'],\
        bind['binsize'+str(man)]['y'], 'b.')
        pylab.title('Binsize: %s' % man)
    pylab.show()

KID = sys.argv[1]
if __name__ == '__main__':

    parser = optparse.OptionParser(usage=\
    "%prog\nUse this script to fold a given lightcurve./nScript requires a Kepler ID as a system argument.")

    parser.add_option('-s','--sharex',\
                        action='store_true',\
                        dest='sharex',\
                        default=False,\
                        help='divide bins into subplots')
    parser.add_option('-o','--oneplot',\
                        action='store_true',\
                        dest='oneplot',\
                        default=True,\
                        help='display bins on same plot')
    parser.add_option('-m','--manual',\
                        type=int,\
                        dest='manual',\
                        default=None,\
                        help='input a manual binning number')
    cChoice = ('LC','SC','')
    parser.add_option('-c','--ctype',\
                        choices=cChoice,\
                        type='choice',\
                        dest='ctype',\
                        default='LC',\
                        help='The expected cadence type'\
                        +', '.join(cChoice)\
                        +') [default: %default]')
    opts, args = parser.parse_args()

runbin(KID,opts.ctype,opts.sharex,opts.oneplot,opts.manual)

