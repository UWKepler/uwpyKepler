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
        Y = kep.qats.qatslc(X.lcFinal,X.KID)
        Y.padLC()
        Y.addNoise()
        Y.runQATS(f=0.01)
        
        coeff = num.polyfit(num.log10(Y.periods),num.log10(Y.snrLC),1)
        outy = scipy.polyval(coeff,num.log10(Y.periods))
        normalizedPower = 10**(outy)
        
        plt.plot(Y.periods,Y.snrFLAT,'r-') 
        plt.plot(Y.periods,Y.snrLC,'b-')
        plt.plot(Y.periods,normalizedPower,'k-')
        plt.setp(plt.gca().set_xscale('log'))
        plt.savefig('sn.'+X.KID+'.png')
        
        dfile = open('signal.'+X.KID+'.data','w')
        print >> dfile,'#',X.KID,'|', Y.periods[num.argmax(Y.snrLC/normalizedPower)],\
                Y.periods[num.argmax(Y.snrLC)], max(Y.snrLC/normalizedPower), max(Y.snrLC)
        for i in range(len(Y.SignalPower)):
            print >> dfile, Y.periods[i],'|',Y.snrLC[i],'|',\
                            Y.snrFLAT[i],'|',normalizedPower[i]
        dfile.close()
        
if __name__ == '__main__':
    
    logFile= sys.argv[1]+'.log'
    main(logFile)