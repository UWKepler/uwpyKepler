
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
    
    Y = kep.qats.qatslc(X.lcFinal,X.KID)
    Y.padLC()
    Y.addNoise()
    Y.runQATS()
    #OutFile = open('ydata','wb')
    #pickle.dump(Y,OutFile)
    #OutFile.close()
    pylab.plot(Y.periods,Y.SignalPower,'b-')
    pylab.setp(pylab.gca().set_xscale('log'))
    pylab.show()