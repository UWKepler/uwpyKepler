#!/astro/apps/pkg/python64/bin/python

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import uwpyKepler as kep
import sys
import numpy as num
import scipy

def main(FileName, kid):
    
    for trial in kep.iofiles.passFileToTrialList(FileName):
        tr = kep.pipelinepars.keptrial(trial)
        kw = kep.keplc.kw(**tr.kw)
        X = kep.keplc.keplc(tr.kid)
        X.runPipeline(kw)
        Y = kep.qats.qatslc(X.lcFinal,X.KID)
        print 'have lightcurve data'
        Y.padLC()
        Y.addNoise()
        print 'have padded and added noise'
        print 'now running Qats'
        Y.runQATS(f=0.01)
        print 'man that took forever'
        coeff = num.polyfit(num.log10(Y.periods),num.log10(Y.snrLC),1)
        outy = scipy.polyval(coeff,num.log10(Y.periods))
        normalizedPower = 10**(outy)
        
        plt.plot(Y.periods,Y.snrFLAT,'r-') 
        plt.plot(Y.periods,Y.snrLC,'b-')
        plt.plot(Y.periods,normalizedPower,'k-')
        plt.setp(plt.gca().set_xscale('log'))
        plt.savefig('sn.'+X.KID+'.png')
        
        print 'initial Qats done, saved figure, now running mainHD'
        nhat0 = kep.nhat_qats_struct.mainHD(Y.tmin[num.argmax(Y.snrLC/normalizedPower)],Y.tmax[num.argmax(Y.snrLC/normalizedPower)], kid)

        
        dfile = open('signal.'+X.KID+'.data','w')
        print >> dfile,'#','KID','|', 'Period at max(Y.snrLC/normalizedPower)',\
                'period at max(Y.snrLC)', 'max(Y.snrLC/normalizedPower)', 'max(Y.snrLC)', '|', 'nhat0'
        print >> dfile,'#',X.KID,'|', Y.periods[num.argmax(Y.snrLC/normalizedPower)],\
                Y.periods[num.argmax(Y.snrLC)], max(Y.snrLC/normalizedPower), max(Y.snrLC),'|',
        for t0 in nhat0:
            print>> dfile, t0, 
        print>> dfile, '\n#'
        print >> dfile,'#','Period ','|','snrLC','|',\
                'snrFlat','|','normalizedPower','|',\
                'tmin','|','tmax','|','q'
        for i in range(len(Y.SignalPower)):
            print >> dfile, Y.periods[i],'|',Y.snrLC[i],'|',\
                            Y.snrFLAT[i],'|',normalizedPower[i],'|',\
                            Y.tmin[i],'|',Y.tmax[i],'|',Y.q[i]
        dfile.close()
        
if __name__ == '__main__':
    
    file = sys.argv[1]
    logFile= file+'.log'
    kid = file[5:]
    main(logFile, kid)