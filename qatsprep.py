import numpy as num

def padLC(lcData):
    
    """
    Pads missing cadences with ones.
    Returns the padded lightcurve.
    """
    
    minCad = min(lcData['cadence'])
    maxCad = max(lcData['cadence'])
    CadSet = num.linspace(minCad,maxCad,1)
    
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
                    #for key in lcData.keys():
                        #if key == 'ydt' or 'yerrdt' or 'y' or 'yerr':
                            #print key
                            #lcData[key] = num.hstack((lcData[key],ones)) 
                            #lcData[key] = lcData[key][sortindex]
                        #if key == 'cadence' or 'x':
                            #lcData[key] = lcData[key][sortindex]
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
    
    #delindex = []
    #for i in lcData['x']:
        #if len(num.where(lcData['x'][i] == lcData['x'][i])[0]) > 1:
            #print len(num.where(lcData['x'][i] == lcData['x'][i])[0])
            #delindex = num.where(lcData['x'][i] == lcData['x'][i])[0][1]
            #lcData['cadence'] = lcData['cadence'].tolist()
            #lcData['x'] = lcData['x'].tolist()
            #lcData['ydt'] = lcData['ydt'].tolist()
            #lcData['yerrdt'] = lcData['yerrdt'].tolist()
            #lcData['y'] = lcData['y'].tolist()
            #lcData['yerr'] = lcData['yerr'].tolist()
            #del(lcData['cadence'][delindex])
            #del(lcData['x'][delindex])
            #del(lcData['ydt'][delindex])
            #del(lcData['yerrdt'][delindex])
            #del(lcData['y'][delindex])
            #del(lcData['yerr'][delindex])
            #lcData['cadence'] = num.array(lcData['cadence'])
            #lcData['x'] = num.array(lcData['x'])
            #lcData['ydt'] = num.array(lcData['ydt'])
            #lcData['yerrdt'] = num.array(lcData['yerrdt'])
            #lcData['y'] = num.array(lcData['y'])
            #lcData['yerr'] = num.array(lcData['yerr'])
            
    return lcData
    
    
    
    
    