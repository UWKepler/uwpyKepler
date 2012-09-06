import uwpyKepler as kep
import numpy as num

class KOI:
    def __init__(self, kid, period, t0, q):
        self.kid    = kid
        self.period = period
        self.t0     = t0
        self.q      = q
        
    def toString(self):
        return 'period=%s, t0=%s, q=%s' % \
            (self.period, self.t0, self.q)

# prompts user for which eData to use
def getPlanets(ed):
    count = 0
    koiTups = []
    print 'select an object to analyze:'
    for key in ed.keys():
        if len(ed[key]) > 0:
            print 'from ' + key
        for koi in ed[key]:
            count += 1
            print '\t' + str(count) + ': ' + koi.toString()
            koiTups.append((count, koi, key))
    prompt = 'input '
    choiceList = num.arange(count) + 1
    for choice in choiceList:
        if choice != choiceList[-1]:
            prompt += str(choice) + ' '
        else:
            prompt += 'or ' + str(choice)
    objChoice = int(raw_input(prompt + '\n> '))
    maskThese = []
    for koi in koiTups:
        if koi[0] == objChoice:
            choice = koi[1]
        elif koi[2] == 'database' \
        or koi[2] == 'transitFit file' \
        or koi[2] == 'binary file':
            maskThese.append(koi[1])

    return choice, maskThese

# retrieves all archived eData
def geteDataDict(lc):
    kid = lc.KID
    eds = {'database':[], 'eDataDiscovery file':[], \
        'transitFit file': [], 'binary file': []}
    eDataDiscovCount = 0
    transitFitFileCount = 0
    eData1 = lc.eData
    if eData1['eDataExists']:
        for key in eData1['KOI'].keys():
            if key.count('fake') == 0: # 'fake' aren't true KOIs
                p = eData1['KOI'][key]['Period']
                t = eData1['KOI'][key]['T0']
                q = eData1['KOI'][key]['Duration']
                koi = KOI(kid, p, t, q)
                eds['database'].append(koi)
    eData2 = kep.analysis.alleDataFromFile(lc.KID)
    for i in range(len(eData2[0])):
        koi = KOI(kid, eData2[0][i], eData2[1][i], eData2[2][i])
        eds['eDataDiscovery file'].append(koi)
    eData3 = kep.analysis.geteDataFromModelFile(lc.KID)
    for i in range(len(eData3[0])):
        koi = KOI(kid, eData3[0][i], eData3[1][i], eData3[2][i])
        eds['transitFit file'].append(koi)
    eData4 = kep.analysis.alleDataFromFile(lc.KID, binary=True)
    for i in range(len(eData4[0])):
        koi = KOI(kid, eData4[0][i], eData4[1][i], eData4[2][i])
        eds['binary file'].append(koi)

    return eds

def koisFromEDataDict(ed):
    kois = []
    for key in ed.keys():
        for koi in ed[key]:
            kois.append(koi)
    
    return kois

# mask transits that would otherwise disrupt plot
# by their depth
def removeOtherTransits(lc, maskThese):
    for i in range(len(maskThese)):
        lc.eData['eDataExists'] = True
        lc.eData['KOI']['fake' + str(i + 1)] = { \
            'Period':maskThese[i].period, \
            'T0':maskThese[i].t0, \
            'Duration':maskThese[i].q}
    kw = kep.quicklc.quickKW()
    lc.runPipeline(kw)
    kep.lcmod.ApplyMask(lc.lcFinal, 'eMask')
