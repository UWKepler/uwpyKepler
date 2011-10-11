import uwpyKepler as kep
import pylab
import numpy as num
import sys
import matplotlib.pyplot as plt

def labels_above(rects):
    for rect in rects:
        height = rect.get_height()
        pylab.text(rect.get_x()+rect.get_width()/2.,1.10*height, '%d'%int(height))

def labels_beneath(rects):
    for rect in rects:
        height = rect.get_height()
        pylab.text(rect.get_x()+rect.get_width()/2.,0.70*height, '%d'%int(height)) 

#skyNumber = sys.argv[1]
skyNumber = [1,2,3]
pylab.figure(1)
for i in skyNumber:
	skyGroup = kep.dbinfo.returnSkyGroupKIDs(i)

	FP = kep.dbinfo.returnKIDsInKEPFP(skyGroup)
	PC = kep.dbinfo.returnKIDsInKEPPC(skyGroup)
	
	FPCoords = kep.dbinfo.returnCoordKID(FP)
	PCCoords = kep.dbinfo.returnCoordKID(PC)
	skyCoords = kep.dbinfo.returnCoordKID(skyGroup)
	
	pylab.plot(skyCoords[0], skyCoords[1], 'k.')
	pylab.plot(FPCoords[0], FPCoords[1], 'ro')
	pylab.plot(PCCoords[0], PCCoords[1], 'co')
	pylab.legend(('SG='+str(skyNumber),'NFP='+str(len(FP)),'NPC='+str(len(PC))))

pylab.figure(2)
lFP = []
lPC = []
lskyGroup = []
for i in skyNumber:
    skyGroup = kep.dbinfo.returnSkyGroupKIDs(i)
    FP = kep.dbinfo.returnKIDsInKEPFP(skyGroup)
    PC = kep.dbinfo.returnKIDsInKEPPC(skyGroup)
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
width = 1.0
indlabels = []
for i in skyNumber:
    indlabels.append('SG'+str(i))
    
rects1 = pylab.bar(ind, lskyGroup, width, color='b', bottom=0.1, label='LCs in SkyGroup')
rects2 = pylab.bar(ind, lPC, width, color='y', bottom=0.1, label='PCs in SkyGroup')
rects3 = pylab.bar(ind, lFP, width, color='r', bottom=0.1, label='FPs in SkyGroup')

pylab.setp(plt.gca().set_yscale('log'))
pylab.xticks(ind+width/2., indlabels, horizontalalignment='center')
pylab.ylabel('objects')
pylab.title('PC and FP vs all')
pylab.legend()

labels_above(rects1)
labels_above(rects2)
labels_beneath(rects3)
pylab.show()

    
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

