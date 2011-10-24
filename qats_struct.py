
import numpy as num
import qats

def gamma_m_list(mlist,tmin,tmax,omDict,q):
    """
    """
    
    outDict = {}
    for m in mlist:
        tempDict = {}
        for MM in omDict[m].keys():
            tempDict[MM] = {}
            for n in range(omDict[m][MM]['k0'],omDict[m][MM]['k1']+1,1):
                tempDict[MM][n] =\
                {'j0':num.int_(max([(m-2L)*tmin,n-tmax])),
                 'j1':num.int_(min([tmax-q+(m-2L)*tmax,n-tmin]))}
        outDict[m] = tempDict
        
    return outDict

def omega_m_list(mlist,MMin,MMax,tmin,tmax,n,q):
    """
    """
    
    outDict = {}
    for m in mlist:
        tempDict = {}
        for MM in range(MMin,MMax+1,1):
            tempDict[MM] = {'k0':num.int_(max([(m-1L)*tmin,n-(MM-m+1)*tmax])),
                            'k1':num.int_(min([-q+m*tmax,n-q-(MM-m)*tmin]))}
                    
        outDict[m] = tempDict
        
    return outDict
    
def qpt_convolve(data,q):
    """
    """
    
    tmp=num.convolve(data,num.ones(q),'valid')
    return tmp
    
def qpt_convolve_multi(data,qlist):
    
    unique_q = list(set(qlist))
    dcDict = {}
    for uq in unique_q:
        dcDict[uq] = qpt_convolve(data,uq)
    
    return dcDict

def fmax_multi(MMin,MMax,Ndata,tmin,tmax,q,dc):
    """
    """
    
    mList = range(1,MMax+1,1)
    
    # Compute omega:
    omDict = omega_m_list(mList,MMin,MMax,tmin,tmax,Ndata,q)
    gmDict = gamma_m_list(mList,tmin,tmax,omDict,q)
    fmaxL = []
    
    for MM in range(MMin,MMax+1,1):
        shape = (MM,Ndata)
        fmnMM = num.zeros(shape)
        for m in mList:
            #om = qats.omega_m(m,MM,tmin,tmax,Ndata,q)
            if m <= MM:
                if m == 1:
                    fmnMM[m-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1] =\
                    dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]
                else:
                    maxArr = []
                    #i = 0
                    for k in range(omDict[m][MM]['k0'],omDict[m][MM]['k1']+1,1):
                        #gm = qats.gamma_m(m,tmin,tmax,om[i],q)
                        #i += 1
                        test = fmnMM[m-2][gmDict[m][MM][k]['j0']:gmDict[m][MM][k]['j1']+1]
                        maxArr.append(max(test))
                        #print gm, gmDict[m][MM][k]['j0'],gmDict[m][MM][k]['j1'], len(test)
                    #print len(om), len(dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]), len(maxArr)
                    #print num.shape(dc), omDict[m][MM]['k0'],omDict[m][MM]['k1'], m, MM
                    fmnMM[m-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1] =\
                    dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]+\
                    num.array(maxArr)
        
                fmax0 = max(fmnMM[MM-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1].ravel())
                #fmax1 = max(fmnMM[MM-1][om].ravel())
                #print om, omDict[m][MM]['k0'],omDict[m][MM]['k1'], num.shape(fmnMM), m, MM
                #print fmax0, fmax1
                fmaxL.append(fmax0)
        
    return fmaxL

def qpt_detect_multi(dcData,pData,Ndata):
    
    spmax = 0
    snr = []
    for ip in pData.keys():
        print ip,pData[ip]['period'],pData[ip]['q']
        fMax =  fmax_multi(pData[ip]['mmin'],\
                           pData[ip]['mmax'],\
                           Ndata,\
                           pData[ip]['tmin'],\
                           pData[ip]['tmax'],\
                           pData[ip]['q'],\
                           dcData[pData[ip]['q']])
        smax = max(fMax)
        imax = smax.argmax()
        MM = pData[ip]['mmin']+imax

        if smax > spmax:
            snr.append(smax/(MM*pData[ip]['q']))
            spmax = smax
            
    return snr