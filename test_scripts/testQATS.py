
import uwpyKepler as kep
import sys
import pylab
import numpy as num
import cPickle as pickle

par = kep.pipelinepars.PipelinePars('ps1')
KIDList = kep.iofiles.readKIDsFromFile('KIDlist')

for trial in par.Trials(KIDList):
    tr = kep.pipelinepars.keptrial(trial)
    kw = kep.keplc.kw(**tr.kw)
    X = kep.keplc.keplc(tr.kid)
    X.runPipeline(kw)
    
    Y = kep.qatspars.qatslc(X.lcFinal,X.KID)
    Y.padLC()
    Y.addNoise()
    Y.runQATS3()
    #OutFile = open('ydata','wb')
    #pickle.dump(Y,OutFile)
    #OutFile.close()
    print len(Y.periods), len(Y.snr0)
    pylab.plot(Y.periods,Y.snr0,'b-')
    pylab.plot(Y.periods,Y.snr1,'r-')
    pylab.setp(pylab.gca().set_xscale('log'))
    pylab.show()