#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse
import os

def findInSubdirectory(filename, subdirectory):
    path = subdirectory
    for root, dirs, names in os.walk(path):
        if filename in names:
            return os.path.join(root, filename)
    raise 'File not found'

def runfold(kid,ctype,auto,hint,man):
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

    if not man:
        try:
            eData_idx = eData['KOI'].keys()[0]
            eData = eData['KOI'][eData_idx]
            x = \
            kep.func.foldPhase(lcData,eData['T0'],eData['Period'])
        except:
            print 'checking for qats file...'
            try:
                qtsfile = open('signal.'+kid+'.data','r')
                top = qtsfile.readlines()[0]
                split = top.split(' ')
                period = eval(split[3])
                x = kep.func.foldPhase(lcData,0,period)
            except:
                try:
                    if hint[-4:] == 'data':
                        qtsfile = open(hint, 'r')
                    else:
                        qtsfilepath =\
                        findInSubdirectory('signal.'+kid+'.data',hint)
                        qtsfile = open(qtsfilepath,'r')
                    top = qtsfile.readlines()[0]
                    split = top.split(' ')
                    period = eval(split[3])
                    x = kep.func.foldPhase(lcData,0,period)
                except:
                    print('no period data found; exiting...\n')
                    sys.exit()
    else:
        x = kep.func.foldPhase(lcData,0,man)
    pylab.title('KID: %s' % kid)
    pylab.xlabel('Phase')
    pylab.ylabel('Flux')
    pylab.plot(x, lcData['ydt'], 'b.')
    pylab.show()

KID = sys.argv[1]
if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog\nUse this script to fold a given lightcurve.\nScript requires a Kepler ID as a system argument.")
    parser.add_option('-a','--auto',\
                        action='store_true',\
                        dest='auto',\
                        default=True,\
                        help='query the database and qats files '\
                        'for a known period')
    parser.add_option('-p','--pathhint',\
                        type=str,\
                        dest='pathhint',\
                        default='/astro/store/student-scratch1',\
                        help='supply a path to period file '\
                        'if outside current directory to '\
                        'dramatically increase speed')
    parser.add_option('-m','--manual',\
                        type=float,\
                        dest='manual',\
                        default=None,\
                        help='input a manual period')
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

runfold(KID,opts.ctype,opts.auto,opts.pathhint,opts.manual)
