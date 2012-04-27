import numpy as num
import math


def Masking(lc):
    #masking out known planet transits
    time = num.array(lc.data['x'])
    
    t0 = []
    period = []
    td = []
    for koi in lc.eData['KOI'].keys():
        t0.append(lc.eData['KOI'][koi]['T0'])
        period.append(lc.eData['KOI'][koi]['Period'])
        td.append(lc.eData['KOI'][koi]['Duration'])
    
    t0 = num.array(t0)
    period = num.array(period)
    td = num.array(td)
    nplanet = len(period)
    #pdpt = num.array([]) currently ignoring quadratic trends in periodicity
    
    #n0 - the (neg) phase of known planet
    n0 = num.floor((num.amin(time)-t0) / period)
    
    #ntrans - max number of transits for this planet within our time range
    ntrans = num.ceil((num.amax(time)-num.amin(time)) / period) + 1
    ntime = len(time)
    #import pdb; pdb.set_trace()
    
    #mask - 2-d array of masks for all planets
    mask = num.zeros((nplanet,ntime))
    for ip in range(nplanet):
        nnp = n0[ip] + num.arange(ntrans[ip])
        ttp = t0[ip] + period[ip]*nnp #+ pdpt[ip]*nnp**2
        for it in range(ntrans[ip]):
            in_transit = num.where(abs(time-ttp[it]) < td[ip]/2.)
            mask[ip][in_transit] = 1
    #mask = mask > 0 #converted into a boolean
    masksum = num.sum(mask,axis=0)
    lc.data['x'] = lc.data['x'][num.where(masksum==0]
    lc.data['y'] = lc.data['y'][num.where(masksum==0]
    return lc