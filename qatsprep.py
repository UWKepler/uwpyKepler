import numpy as num
import pylab

def padLC(lcData):
    """
    Pads missing cadences with ones.
    Returns the padded lightcurve.
    """
    newcadence = ['Placeholder']
    newtstamp = ['Placeholder']
    gaps = True
    efficiencynum = 0
    diffs = []
    
    for i in range(len(lcData['cadence'])):
        if i < len(lcData['cadence'])-1:
            if lcData['cadence'][i+1]-lcData['cadence'][i] > 1:
                diffs.append(i)
                
    gapnum = len(diffs)
    
    while gapnum > 0:
        for i in range(len(lcData['cadence'])):
            if i < len(lcData['cadence'])-1:
                diffcadence = lcData['cadence'][i+1]-lcData['cadence'][i]
                if diffcadence > 1:
                    gapnum -= 1
                    newcadence = []
                    newtstamp = []
                    diff = lcData['cadence'][i+1]-lcData['cadence'][i]
                    cnewidx = num.arange(diff-1)+i+1 #Does not work after first gap: indices off. Problem addressed below
                    tnewidx = num.arange(diff-1)
                    for j in range(len(cnewidx)):
                        if j == 0:
                            newcadence.append(cnewidx[j])
                            newtstamp.append(lcData['x'][i]+difftstamp)
                        else: #Prevents appending same timestamp when gap is larger than 1 point
                            newcadence.append(cnewidx[j])
                            newtstamp.append(newtstamp[j-1]+difftstamp)
                    ### Put inside for loop to reset lcData at each iteration, solving index problem above ###
                    newcadence = num.array(newcadence, dtype='int64')
                    newtstamp = num.array(newtstamp, dtype='float')
                    ones = num.zeros(len(newcadence))+1e0
                    lcData['cadence'] = num.hstack((lcData['cadence'],newcadence))
                    lcData['x'] = num.hstack((lcData['x'],newtstamp))
                    sortindex = lcData['cadence'].argsort()
                    lcData['ydt'] = num.hstack((lcData['ydt'],ones)) 
                    lcData['yerrdt'] = num.hstack((lcData['yerrdt'],ones)) 
                    lcData['y'] = num.hstack((lcData['y'],ones)) 
                    lcData['yerr'] = num.hstack((lcData['yerr'],ones)) 
                    
                    lcData['ydt'] = lcData['ydt'][sortindex]
                    lcData['yerrdt'] = lcData['yerrdt'][sortindex]
                    lcData['y'] = lcData['y'][sortindex]
                    lcData['yerr'] = lcData['yerr'][sortindex]
                    lcData['cadence'] = lcData['cadence'][sortindex]
                    lcData['x'] = lcData['x'][sortindex]
                    break
                else:
                    if efficiencynum == 0:
                        difftstamp = lcData['x'][i+1]-lcData['x'][i]
                        efficiencynum += 1
                        
    lcData['flagflat'] = num.zeros(len(lcData['y']))
    lcData['flagflat'][num.where(lcData['y'] == 1)[0]] = 1
    
    #x1 = pylab.subplot(211)
    #pylab.plot(lcData['x'], lcData['y'], 'bo')
    #x2 = pylab.subplot(212, sharex=x1)
    #pylab.plot(lcData['x'], lcData['flagflat'], 'ro')
    #pylab.show()
            
    return lcData
    
    
    
    
    