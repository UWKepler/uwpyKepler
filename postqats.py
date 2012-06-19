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
    kid = str(kid)
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
    
# used in flagQats fitting
def firstFareyVal(p0,p1,p2):
    m   = 1. # start on right side for when p0 = p1
    n   = 1.
    ln  = 0. # 'left numerator'
    ld  = 1. # 'left denominator'
    rn  = 1. # 'right numerator'
    rd  = 1. # 'right denominator'
    pRatio1 = p1/p0
    pRatio2 = p2/p0
    if pRatio2 > 1:
        temp = pRatio1
        pRatio1 = 1./pRatio2
        pRatio2 = 1./temp
    while m/n < pRatio1 or m/n > pRatio2:
        if m/n < pRatio1:
            ln = m
            ld = n
            m += rn
            n += rd
        else:
            rn = m
            rd = n
            m += ln
            n += ld
    return m, n
