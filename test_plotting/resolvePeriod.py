#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import optparse
import os

def makeMovie(fileNames):
    kid = fileNames[0].split('_')[0]
    p0 = fileNames[0].split('_')[1]
    pf = fileNames[-1].split('_')[1]
    movieName = kid + '_' + p0 + '_to_' + pf + '_foldAnimation'
    os.system("mencoder 'mf://" + kid + "_*.png' -mf type=png:fps=10 " \
        "-ovc lavc -lavcopts vcodec=wmv2 -oac copy " \
        "-o " + movieName + ".avi")
    for file in fileNames:
        os.remove(file)

class interactivePeriodResolver:
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
    parser.add_option('-a','--animate',\
                        type=float,\
                        dest='animate',\
                        default=True,\
                        help='input a frame number to make\n' \
                        'a movie with one period per frame;\n' \
                        'period range is still defined by' \
                        'zoomfactor')
    parser.add_option('-i','--interactive',\
                        action='store_true',\
                        dest='interactive',\
                        default=False,\
                        help='enter interactive mode')
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

if opts.manual:
    p0 = opts.manual
else:
    try:
        eKeys = eData['KOI'].keys()
        if len(eKeys) > 1:
            print 'multiple catalogued periods:'
            for key in eKeys:
                print '\t' + key + ' period: ' + str(eData['KOI'][key]['Period'])
            p0 = eData['KOI'][eKeys[0]]['Period']
            print 'using period for ' + eKeys[0] + ': ' + str(p0)
        else:
            eData_idx = eKeys[0]
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
if opts.interactive:
    phase = kep.func.foldPhase(lcData['x'],0,p0)
    pylab.plot(phase, lcData['ydt'], 'b.')
    pFinder = interactivePeriodResolver(p0, 2*t_dur/10., \
        lcData['x'], lcData['ydt'])
    pylab.show()
    sys.exit()
if opts.animate:
    nFrames = opts.animate
    testPs = num.arange(p0-t_dur,p0+t_dur,2*t_dur/nFrames)
    fileNames = []
else:
    testPs = num.arange(p0-t_dur,p0+t_dur,2*t_dur/10.)

for plot in range(len(testPs)):
    phase = kep.func.foldPhase(lcData['x'],0,testPs[plot])
    pylab.figure(plot)
    pylab.title('KID: %s, P0: %s, ZoomFac: %s' % \
                (KID, testPs[plot], zf))
    pylab.xlabel('Phase')
    pylab.ylabel('Flux')
    pylab.plot(phase, lcData['ydt'], 'b.')
    if opts.animate or not opts.plot:
        fileName = '%s_%s_%s.png' % (KID, testPs[plot], zf)
        pylab.savefig(fileName)
        if opts.animate:
            fileNames.append(fileName)

if opts.plot:
    pylab.show()
elif opts.animate:
    makeMovie(fileNames)


