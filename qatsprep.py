import numpy as num

def padLC(lcData):
    """
    Pads missing cadences with ones.
    Returns the padded lightcurve.
    """
                
    newcadence = []
    newtstamp = []
    oldcadence = lcData['cadence']
    efficiencynum = 1
    iteration = 1
    for i in range(len(oldcadence)):
        if i < len(oldcadence)-1:
            diffcadence = lcData['cadence'][i+1]-lcData['cadence'][i]
            if diffcadence > 1:
                diff = lcData['cadence'][i+1]-lcData['cadence'][i]
                cnewidx = num.arange(diff-1)+i+1
                tnewidx = num.arange(diff-1)+1
                for j in cnewidx:
                    if iteration == 0:
                        newcadence.append(j)
                        newtstamp.append(lcData['x'][i]+difftstamp)
                        iteration =+ 1
                    else: #Code in progress to not keep appending same timestamp when gap is larger than 1 point
                        newcadence.append(j)
                        newtstamp.append(lcData['x'][i+j]+difftstamp)
            else:
                while efficiencynum > 0:
                    difftstamp = lcData['x'][i+1]-lcData['x'][i]
                    efficiencynum -= 1

    newcadence = num.array(newcadence, dtype='int64')
    newtstamp = num.array(newtstamp, dtype='int64')
    ones = num.zeros(len(newcadence))+1e0
    #lcData['x'][newcadence] = 1
    lcData['cadence'] = num.hstack((oldcadence,newcadence))
    lcData['x'] = num.hstack((lcData['x'],newtstamp))
    sortindex = lcData['cadence'].argsort()
    #lcData['cadence'] = lcData['cadence'][sortindex]
    #lcData['x'] = lcData['x'][sortindex]
    for key in lcData.keys():
        if key == ['ydt', 'yerrdt', 'y', 'yerr']:
            lcData[key] = num.hstack((lcData[key],ones)) 
            lcData[key] = lcData[key][sortindex]
        if key == ['cadence', 'x']:
            lcData[key] = lcData[key][sortindex]
    
    return lcData
    
    
    
    
    