import sys
import uwpyKepler as kep
import numpy as num
import pylab
import optparse

kid = sys.argv[1]
dwindowsize = 50
dorder = 3
usage = """%prog: Use to detrend python data.\nInput a KID first.\nDetrend function properties...\n-Window Size\n\tDefault: 50 data points\n\tSet with: -w\n-Polynomial Order\n\tDefault: 3\n\tSet with: -o"""
parser = optparse.OptionParser(usage)
parser.add_option('-w', '--window', type=int, default = dwindowsize, help='Detrending window size.')
parser.add_option('-o', '--order', type=int, default = dorder, help='Detrending function polynomial order.')
options, kid = parser.parse_args()
kid = kid[0]
#kid = '10341831'
#kid = ''

d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.05,2,0.05)
d3 = kep.io.FlagOutliers(pd,10,4)
d4 = kep.proc.detrendData(d3,options.window,options.order)

pylab.plot(d1['x'],d1['y'],'b.')

for portion in d4.keys():
    #pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    pylab.plot(d4[portion]['x'], d4[portion]['Correction'], 'k-')

pylab.title('Object ' + kid)
pylab.xlabel('Time in BJD')
pylab.ylabel('Flux')
#pylab.legend(('Data', 'Correction Function'),'upper right', shadow=True)
pylab.show()
#pylab.savefig('testplota.png')

pylab.figure(2)
#pylab.plot(d4[portion]['x'],d4[portion]['y'],'y.')

d7 = kep.proc.stackPortions(d4)
#d8 = kep.proc.cutOT(d7)
#pylab.plot(d8['x'],d8['y'],'b.')

pylab.plot(d7['x'],d7['y'],'y.')

d8 = kep.proc.onlyOutliers(d7)

pylab.plot(d8['x'],d8['y'],'b.')

d9 = kep.proc.stackPortions(d4)
d10 = kep.proc.onlyTransits(d9)

pylab.plot(d10['x'],d10['y'],'r.')

pylab.title('Object ' + kid)
pylab.xlabel('Time in BJD')
pylab.ylabel('Flux')
#pylab.legend(('Data', 'Outliers','Transits'),'lower right', shadow=True)

pylab.show()
#pylab.savefig('testplotb.png')
