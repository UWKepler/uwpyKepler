import uwpyKepler as kep
import pylab
import numpy as num
import sys

skyNumber = sys.argv[1]

skyGroup = kep.dbinfo.returnSkyGroupKIDs(skyNumber)

FP = kep.dbinfo.returnKIDsInKEPFP(skyGroup)
PC = kep.dbinfo.returnKIDsInKEPPC(skyGroup)

FPCoords = kep.dbinfo.returnCoordKID(FP)
PCCoords = kep.dbinfo.returnCoordKID(PC)
skyCoords = kep.dbinfo.returnCoordKID(skyGroup)

pylab.plot(FPCoords[0], FPCoords[1], 'r.')
pylab.plot(PCCoords[0], PCCoords[1], 'r.')
pylab.plot(skyCoords[0], skyCoords[1], 'k.')
pylab.show()

#FPKids = kep.dbinfo.getKIDsFP()
#PCKids = kep.dbinfo.getKIDsPC()
#allKids = kep.dbinfo.getKIDsSource()

#FPCoords = kep.dbinfo.returnCoordKID(FPKids)
#PCCoords = kep.dbinfo.returnCoordKID(PCKids)
#allCoords = kep.dbinfo.returnCoordKID(allKids)

#pylab.plot(FPCoords[0], FPCoords[1], 'r.')
#pylab.plot(PCCoords[0], PCCoords[1], 'r.')
#pylab.plot(allCoords[0], allCoords[1], 'k.')
#pylab.show()

