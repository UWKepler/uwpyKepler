import sys
import uwpyKepler as kep
import numpy as num
import pylab

def mov_avg(data,window):
    
    mvavg = []
    npoints = len(data)
    for i in range(npoints-1):
        j1  = max(0,i-window)
        j2 = min(i+window,npoints-1)
        mvavg.append(num.median(data[j1:j2]))

    return num.array(mvavg)

def smoothing(data,window):
    
    mvavg = []
    npoints = len(data['x'])
    data['x'].mask = data['UnMasked']
    data['y'].mask = data['UnMasked']
    data['yerr'].mask = data['UnMasked']
    
    mvavg1 = mov_avg(data['y'],window)
    mvavg1 = num.append(mvavg1,mvavg1[-1])
    
    #mvavg1 = data['y']

    diff1 = num.diff(mvavg1)
    diff1 = num.hstack((diff1,diff1[-1]))
    diff2 = num.diff(diff1)
    diff2 = num.hstack((diff2[-1],diff2))
    tgap = []
    fgap = []
    print type(diff1)
    sig = kep.io.compute1Sigma(diff2)
    for i in range(npoints-1):
        dt = data['x'][i+1] - data['x'][i]
        dy1 = diff1[i]
        dy2 = diff2[i]
        if dt > 0.1:
            tgap.append(i)
        if num.abs(dy2) > 7*sig:
            try:
                if i-2 == fgap[-1] or i-1 == fgap[-1]:
                    pass
                elif num.abs(diff[i+1]) > 7*sig:
                    pass
                else:
                    fgap.append(i)
                    
            except:
                fgap.append(i)
    pylab.plot(data['x'],data['y'],'b.')
    for el in tgap:
        print el, 'time gap', data['x'][el-1],data['x'][el]
        pylab.plot(data['x'][el],data['y'][el],'gs',markersize=10,linewidth=5)
    for el in [el for el in fgap if el not in tgap]:
        #print el, 'flux gap', data['y'][el-1],data['y'][el]
        print el, 'flux gap', diff1[el-2], diff1[el-1],diff1[el], sig, num.std(diff2)
        pylab.plot(data['x'][el],data['y'][el],'rs',markersize=10,linewidth=5)
    
    #pylab.plot(data['x'],mvavg,'k-',linewidth=3)
    #pylab.plot(data['x'],mvavg1,'g-',linewidth=3)
    #pylab.plot(data['x'],mvavg2,'r-',linewidth=3)
    #pylab.plot(data['x'],1.92e7+num.abs(diff1),'r.')
    #pylab.plot(data['x'],1.92e7+num.abs(diff2),'g.')
    #pylab.plot(data['x'],1.92e7+num.zeros(len(data['x']))+num.median(num.abs(diff2))+6*sig,'g-')
    
    pylab.show()

kid = sys.argv[1]
#kid = '9663113'
kid = '10341831'
#kid = ''

d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)

smoothing(d2,5)

