#!/astro/apps/pkg/python64/bin/python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import uwpyKepler as kep
import sys
import numpy as num
import scipy

def main(FileName):
    
    for trial in kep.iofiles.passFileToTrialList(FileName):
        tr = kep.pipelinepars.keptrial(trial)
        kw = kep.keplc.kw(**tr.kw)
        X = kep.keplc.keplc(tr.kid)
        X.runPipeline(kw)
        lc = X.lcData.lcData
        Y = kep.qats.qatslc(X.lcFinal,X.KID,maske=kw.maske)
        
        #OLD MASKING METHOD
        #idx = num.where((X.lcFinal['eMask'] == False))[0]
        #Y.padLC(flagids=idx)
        
        Y.padLC()
        Y.addNoise()
        Y.runQATS(f=0.01)
        
        P,ignored = calcPeriodogram(Y.periods, lc['x'], lc['ydt'], lc['yerr']) ### ADDED ###
                
        coeff = num.polyfit(num.log10(Y.periods),num.log10(Y.snrLC),1)
        outy = scipy.polyval(coeff,num.log10(Y.periods))
        normalizedPower = 10**(outy)
                
        coeff2 = num.polyfit(num.log10(Y.periods),num.log10(Y.snrFLAT),1) ### ADDED ###
        outy2 = scipy.polyval(coeff2,num.log10(Y.periods)) ### ADDED ###
        normalizedPower2 = 10**(outy2) ### ADDED ###
                
        fittedr=Y.snrFLAT/normalizedPower2 ### ADDED ###
        fittedb=Y.snrLC/normalizedPower2 ### ADDED ###
        normalizedPower3 = normalizedPower2/normalizedPower2 ### ADDED ###
                
        plt.subplot(211) ### ADDED ###
        plt.title('KID = '+X.KID)
        plt.plot(Y.periods,Y.snrFLAT,'r-') 
        plt.plot(Y.periods,Y.snrLC,'b-')
        plt.plot(Y.periods,normalizedPower,'k-')
                #plt.plot(Y.periods,normalizedPower2, 'k-') ### ADDED ###
        plt.setp(plt.gca().set_xscale('log'))
        plt.ylabel('Signal Power')
        plt.subplot(212)
        plt.plot(Y.periods,P,'g-') ### ADDED ###
        plt.setp(plt.gca().set_xscale('log')) ### ADDED ###
        plt.xlabel('Period (days)') ### MOVED ###
        plt.ylabel('F Transform')
        plt.savefig('sn.test.'+X.KID+'_'+kw.maske+'.png') ### MOVED ###
        #coeff = num.polyfit(num.log10(Y.periods),num.log10(Y.snrLC),1)
        #outy = scipy.polyval(coeff,num.log10(Y.periods))
        #normalizedPower = 10**(outy)
        
        #plt.plot(Y.periods,Y.snrFLAT,'r-') 
        #plt.plot(Y.periods,Y.snrLC,'b-')
        #plt.plot(Y.periods,normalizedPower,'k-')
        #plt.setp(plt.gca().set_xscale('log'))
        #plt.savefig('sn.'+X.KID+'.png')
        
        dfile = open('signal.'+X.KID+'.data','w')
        print >> dfile,'#',X.KID,'|', Y.periods[num.argmax(Y.snrLC/normalizedPower)],\
                Y.periods[num.argmax(Y.snrLC)], max(Y.snrLC/normalizedPower), max(Y.snrLC)
        for i in range(len(Y.SignalPower)):
            print >> dfile, Y.periods[i],'|',Y.snrLC[i],'|',\
                            Y.snrFLAT[i],'|',normalizedPower[i]
        dfile.close()


def calcPeriodogram(periods, times, data, errors):
    # weights
    W      = num.sum(1. / errors**2)
    w      = 1. / W / errors**2
    wy     = w * data
   
    # makes times a phase, 0..1
    phases = num.outer(1./periods, times) - num.outer(1./periods, times).astype(num.int)

    # makes wt which runs from 0..2pi
    freqt  = 2 * num.pi * phases
   
    # solve for tau
    cfreqt  = num.cos(freqt)
    sfreqt  = num.sin(freqt)
    csfreqt = cfreqt * sfreqt
    ccfreqt = cfreqt * cfreqt
    ssfreqt = 1 - ccfreqt

    C       = num.sum(w * cfreqt, axis = 1)
    S       = num.sum(w * sfreqt, axis = 1)
    CCh     = num.sum(w * ccfreqt, axis = 1)
    SSh     = num.sum(w * ssfreqt, axis = 1)
    CSh     = num.sum(w * csfreqt, axis = 1)

    CS      = CSh - C * S
    CC      = CCh - C * C
    SS      = SSh - S * S
    wtau    = num.arctan(2 * CS / (CC - SS))

    freqt  -= wtau[:,num.newaxis]

    # resolve for matrices
    cfreqt  = num.cos(freqt)
    sfreqt  = num.sin(freqt)
    csfreqt = cfreqt * sfreqt
    ccfreqt = cfreqt * cfreqt
    ssfreqt = 1 - ccfreqt

    Y     = num.sum(wy)
    C     = num.sum(w * cfreqt, axis = 1)
    S     = num.sum(w * sfreqt, axis = 1)

    YYh   = num.sum(w * data**2)
    YCh   = num.sum(wy * cfreqt, axis = 1)
    YSh   = num.sum(wy * sfreqt, axis = 1)
    CCh   = num.sum(w * ccfreqt, axis = 1)
    SSh   = num.sum(w * ssfreqt, axis = 1)
    CSh   = num.sum(w * csfreqt, axis = 1)

    YY    = YYh - Y * Y
    YC    = YCh - Y * C
    YS    = YSh - Y * S
    CC    = CCh - C * C
    SS    = SSh - S * S
    CS    = CSh - C * S

    p     = (YC * YC / CC + YS * YS / SS) / YY
    return p, wtau

        
if __name__ == '__main__':
    
    logFile= sys.argv[1]+'.log'
    main(logFile)
