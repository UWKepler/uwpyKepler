
import uwpyKepler as kep
import sys
import pylab
import numpy as num
import scipy
import cPickle as pickle

LogFile = open('testQATS.log','w')
LogFile = open('testQATS.log','a')

for trial in kep.iofiles.passFileToTrialList('test10.log'):
    print trial
    tr = kep.pipelinepars.keptrial(trial)
    kw = kep.keplc.kw(**tr.kw)
    X = kep.keplc.keplc(tr.kid)
    X.runPipeline(kw)
    
    Y = kep.qats.qatslc(X.lcFinal,X.KID)
    Y.padLC()
    Y.addNoise()
    Y.runQATS(f=0.01)
    #OutFile = open('ydata','wb')
    #pickle.dump(Y,OutFile)
    #OutFile.close()
    
    coeff = num.polyfit(num.log10(Y.periods),num.log10(Y.snrLC),1)
    outy = scipy.polyval(coeff,num.log10(Y.periods))
    normalizedPower = 10**(outy)
    
    pylab.plot(Y.periods,Y.snrFLAT,'r-') 
    pylab.plot(Y.periods,Y.snrLC,'b-')
    pylab.plot(Y.periods,normalizedPower,'k-')
    pylab.setp(pylab.gca().set_xscale('log'))
    #pylab.show()
    pylab.savefig('sn.'+X.KID+'.png')
    dfile = open('signal.'+X.KID+'.data','w')
    print >> LogFile, X.KID,'|', Y.periods[num.argmax(Y.snrLC/normalizedPower)],\
             Y.periods[num.argmax(Y.snrLC)], max(Y.snrLC/normalizedPower), max(Y.snrLC)
    for i in range(len(Y.SignalPower)):
        print >> dfile, Y.periods[i],'|',Y.snrLC[i],'|',\
                        Y.snrFLAT[i],'|',normalizedPower[i]
    dfile.close()
    
LogFile.close()