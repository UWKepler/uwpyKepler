import numpy as num
from lcmod import returnData
from dbinfo import returnLOGG, returnRstar
import qats
import pylab

def getListIndicies(Array,ListValues):
    """
    Return a list of indices where a list of values exists
    in a given array
    """ 
    
    lArray = Array.tolist()
    lLV = ListValues.tolist()
    
    indX = [lArray.index(x) if x in lLV else None for x in lLV]
    
    return indX

def stellar_dens(KID):
    """
    compute stellar density
    """
    
    Logg = returnLOGG(KID)
    Rstar = returnRstar(KID)
    Grav = 6.674e-8
    Msol = 1.99e33
    Rsol = 6.96e10
    rho_s = 4e0*num.pi*(10**(Logg))/(Grav*Rstar*Rsol)
    return rho_s

def tdur(rho_s,b,period):
    """
    compute transit duration given stellar density,
    impact parameter and period
    """
    
    p_sec = period*86400e0
    Grav = 6.674e-8
    tdur = ((3e0*num.pi/Grav)**(1e0/3e0))*\
         ((1e0/num.pi)*(p_sec/rho_s)**(1e0/3e0))*\
         (num.sqrt(1e0-b**2))
         
    tdur = tdur/86400e0
    
    return tdur

class qatslc:

    def __init__(self,lcData,KID):
        
        self.lcData = lcData
        self.kid = KID
        self.status = 'Initial Input'

    def padLC(self,**kwargs):
        """
        Pads missing cadences with ones.
        Returns the padded lightcurve.
        
        FlagIDs must corresepond to input lcData indices
        i.e. indices of the array from the pipline
        """
        
        #Default = no flaged points
        FlagIDs = []
        for key in kwargs:
            if key.lower() == 'flagids':
                FlagIDs = kwargs[key]

        # the eclipse lightcurve data are the best place to start
        x,y,yerr,cad = returnData(self.lcData,'elc')
        
        if len(FlagIDs) > 0:
            if len(x)-1 < max(FlagIDs):
                raise NameError("len(x) < max(FlagIds) "+\
                str(len(x))+" < "+str(len(FlagIds)) )
            for el in FlagIDs:
                x.pop(el)
                y.pop(el)
                yerr.pop(el)
                cad.pop(el)
        else:
            pass
        
        x0 = min(x)
        x1 = max(x)
        c0 = min(cad)
        c1 = max(cad)
        dx = num.median(x[1L:len(x)]\
                            -x[0L:len(x)-1])
        Tcad = num.median(dx)
        xcomplete = num.arange(x0,x1,Tcad)
        cadcomplete = num.arange(c0,c1+1,1)
        
        if len(xcomplete) != len(cadcomplete):
            print len(xcomplete), len(cadcomplete)
            raise NameError("Mismatch in Cadence and Time Intervals")
        
        #missing = num.array(list(set.difference(set(cadcomplete),set(cad))))
        #missingIDX = getListIndicies(cadcomplete,missing)
        existingIDX = getListIndicies(cadcomplete,cad)
    
        zeros = num.zeros(len(xcomplete))
        yerrcomplete = zeros                #padding errors with 0
        padflag = zeros
        ycomplete = zeros+1e0               #padding lc with 1
        
        # padded datasets
        # re-using original time-stamps for existing cadences
        xcomplete[existingIDX] = x
        ycomplete[existingIDX] = y
        yerrcomplete[existingIDX] = yerr
        qflag = zeros[existingIDX] = 1
        self.lcData = {'x':xcomplete,\
                       'y':ycomplete,\
                       'yerr':yerrcomplete,\
                       'padflag':padflag}
        self.status = 'Padded Lightcurve'

    def addNoise(self,**kwargs):
        """
        add noise to padded data
        """
        
        # Default Sigma
        Sigma = num.std(self.lcData['y'])
        
        for key in kwargs:
            if key.lower() == 'noise':
                Sigma = kwargs[key]
            else:
                continue

        NoiseIDs = num.where(self.lcData['yerr'] == 0e0)[0]
        self.lcData['y'][NoiseIDs] += Sigma*num.random.randn(len(NoiseIDs))
        self.status = 'QATS Ready'
        
    def runQATS(self, **kwargs):
        
        NPoints = len(self.lcData['x'])
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])

        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        f = 0.005e0
        b = 0e0
        rho_s = stellar_dens(self.kid)

        self.pmin = pmin
        self.pmax = pmax
        self.f = f
        self.nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))
        period0 = pmin/(1+f/2)*(1-f/2)

        period = []
        speriod = []
        mbest = []
        qbest = []

        for i in range(self.nperiod):
            #print period0
            spmax = 0e0
            period0 = (period0*(1+f/2)/(1-f/2))
            period.append(period0)
            tmin = num.floor(period0*(1e0-f/2e0))
            tmax = num.ceil(period0*(1e0+f/2e0))
            q = num.floor(tdur(rho_s,b,period0)/self.dt)
            MM, nhat, smax, dc = \
            qats.qpt_detect(self.lcData['y'],tmin,tmax,q)
            speriod.append(smax)
            if smax > spmax:
                mbest.append(MM)
                qbest.append(q)
                spmax = smax

        period = num.array(period)
        speriod = num.array(speriod)
        mbest = num.array(mbest)
        qbest = num.array(qbest)
        snr = speriod/num.sqrt(mbest*qbest)
        self.period = period
        self.snr = snr
        
        
