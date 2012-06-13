import uwpyKepler as kep
import numpy as num
import pdb

# returns KID from QATS .data file
# cursor must be at start of file
def getKID(dfile):
    lines = dfile.readline()
    KID = int(lines.split('|')[0].split('#')[1])
    return KID

def getBestPeriod(dfile):
    line = dfile.readline()
    return float(line.split('|')[1].split(' ')[1])

def getBestPeriodByKID(kid):
    sgStr = str(kep.dbinfo.getSkyGroup(kid)).zfill(3)
    dFileName = '/astro/store/student-scratch1/'\
                'johnm26/SPRING_BREAK_RUNS/SG' + sgStr + \
                '/signal.SG' + sgStr + \
                '.unflipped.' + kid + '.data'
    dfile = open(dFileName, 'r')
    return getBestPeriod(dfile)

def getQatsData(dfile):
    lines = dfile.readlines()
    if lines[0].startswith('#'):
        del lines[0]
    
    
    #Defining lists to be populated with data
    periods = []
    snrLC=[]
    snrFLAT=[]
    
    #Populating lists with data from file
    for line in lines:
        row = line.split(' | ')
        periods.append(float(row[0]))
        snrLC.append(float(row[1]))
        snrFLAT.append(float(row[2]))

    periods = num.array(periods)
    snrLC = num.array(snrLC)
    snrFLAT = num.array(snrFLAT)
    
    # get actual Signal to Noise Ratio
    snr = snrLC/snrFLAT
    
    return periods, snr, snrLC, snrFLAT