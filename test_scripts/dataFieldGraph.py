import uwpyKepler as kep
import pylab
import numpy as num
import sys
import matplotlib.pyplot as plt

#skyNumber = sys.argv[1]
skyNumber = [1,2,3,4]
#for i in range(skyNumber):
	#skyGroup = kep.dbinfo.returnSkyGroupKIDs(skyNumber)

	#FP = kep.dbinfo.returnKIDsInKEPFP(skyGroup)
	#PC = kep.dbinfo.returnKIDsInKEPPC(skyGroup)
	
	#FPCoords = kep.dbinfo.returnCoordKID(FP)
	#PCCoords = kep.dbinfo.returnCoordKID(PC)
	#skyCoords = kep.dbinfo.returnCoordKID(skyGroup)
	
	#pylab.plot(skyCoords[0], skyCoords[1], 'k.')
	#pylab.plot(FPCoords[0], FPCoords[1], 'ro')
	#pylab.plot(PCCoords[0], PCCoords[1], 'co')
	#pylab.legend(('SG='+str(skyNumber),'NFP='+str(len(FP)),'NPC='+str(len(PC))))
lFP = []
lPC = []
lskyGroup = []
print type(lFP)
for i in skyNumber:
    skyGroup = kep.dbinfo.returnSkyGroupKIDs(i)
    
    FP = kep.dbinfo.returnKIDsInKEPFP(skyGroup)
    PC = kep.dbinfo.returnKIDsInKEPPC(skyGroup)
    print type(lFP)
    if len(FP) > 0:
        lFP.append(len(FP))
    else:
        lFP.append(0)
    if  len(PC) > 0:
        lPC.append(len(PC))
    else:
        lPC.append(0)
    
    lskyGroup.append(len(skyGroup))
    
   
    

N=len(skyNumber)
ind = num.arange(N)+1
width = .35

#p1 = plt.bar(ind, lskyGroup, width, color='b')
#p2 = plt.bar(ind, lFP,  width, color='r', bottom=lFP)
#p3 = plt.bar(ind, lPC, width, color='y')
plt.plot(ind,lskyGroup,'k.')
plt.plot(ind,lFP,'b.')
plt.plot(ind,lPC,'r.')
plt.setp(plt.gca().set_yscale('log'))
plt.ylabel('objects')
plt.title('PC and FP vs all')
#plt.xticks(ind+width/2., ('G1', 'G2', 'G3', 'G4', 'G5') )
#plt.yticks(np.arange(0,81,10))

plt.show()

    
    #####FPCoords = kep.dbinfo.returnCoordKID(FP)
    #####PCCoords = kep.dbinfo.returnCoordKID(PC)
    #####skyCoords = kep.dbinfo.returnCoordKID(skyGroup)
    
    #####pylab.plot(skyCoords[0], skyCoords[1], 'k.')
    #####pylab.plot(FPCoords[0], FPCoords[1], 'ro')
    #####pylab.plot(PCCoords[0], PCCoords[1], 'co')
    #####pylab.legend(('SG='+str(skyNumber),'NFP='+str(len(FP)),'NPC='+str(len(PC))))
    #####pylab.show()

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

