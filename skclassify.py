from uwpyKepler.postqats import *
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
import os

TYPES = ('AverageObjects', 'VariableObjects', 'SinglePlanets', \
    'MultiPlanets', 'Binaries')
LISTS = ('AverageObjects.txt', 'VariableObjects.txt', \
    'SinglePlanets.txt', 'MultiPlanets.txt', 'EBList.txt')

def extraTreePredictions(kidCapsule):
    classifier = ExtraTreesClassifier(n_estimators=50, \
        max_depth=None, min_samples_split=1, random_state=0)
    classifier.fit(kidCapsule.train_M, kidCapsule.trainLabels)
    #import pdb; pdb.set_trace()
    return classifier.predict(kidCapsule.M)

def allClassifierPredictions(kidCapsule):
    decisionTree   = DecisionTreeClassifier(max_depth=None, min_samples_split=1, random_state=0)
    randomForest   = RandomForestClassifier(n_estimators=50, max_depth=None, min_samples_split=1, random_state=0)
    extraTrees     = ExtraTreesClassifier(n_estimators=50,   max_depth=None, min_samples_split=1, random_state=0)
    gradientBoost  = GradientBoostingClassifier(n_estimators=50, max_depth=1, learn_rate=1.0, random_state=0)
    decisionTree.compute_importances  = True
    randomForest.compute_importances  = True
    extraTrees.compute_importances    = True
    gradientBoost.compute_importances = True
    decisionTree.fit(kidCapsule.train_M, kidCapsule.trainLabels)
    randomForest.fit(kidCapsule.train_M, kidCapsule.trainLabels)
    extraTrees.fit(kidCapsule.train_M, kidCapsule.trainLabels)
    gradientBoost.fit(kidCapsule.train_M, kidCapsule.trainLabels)
    print decisionTree.feature_importances_
    print randomForest.feature_importances_
    print extraTrees.feature_importances_
    print gradientBoost.feature_importances_
    dt_pred = decisionTree.predict(kidCapsule.M)
    rf_pred = randomForest.predict(kidCapsule.M)
    et_pred = extraTrees.predict(kidCapsule.M)
    gb_pred = gradientBoost.predict(kidCapsule.M)
    #import pdb; pdb.set_trace()
    return dt_pred, rf_pred, et_pred, gb_pred

def mapPredictions(kidCapsule, predictions):
    Map = lambda index: TYPES[index]
    predictions_strings = map(Map, predictions.astype(int))
    return zip(kidCapsule.classifyKids, predictions_strings)

class QatsEnsemble(object):
    def __init__(self, name, infile, kids = None, matrix = None):
        self.name   = name
        self.infile = infile
        self.kids   = kids
        self.matrix = matrix

    def readKids(self):
        self.kids = num.loadtxt(self.infile, unpack = True)
        self.kids = self.kids.astype(int)

    def readMatrix(self):
        if self.kids == None:
            self.readKids()
        self.matrix, self.kidsFound = readInputs(self.kids)

class KidsLabelCapsule:
    #def __init__(self, balance = True, **kwargs):
        ## scikit recommends balancing dataset prior to learning
        #self.balance = balance
        #self.setTrainingData(**kwargs)
        #self.setClassifyData(classifyKids)
        
    def __str__(self):
        string = ''
        for classification, ensemble in zip(TYPES, self.ensembles):
            string += classification + ': ' + \
                str(len(ensemble.kidsFound)) + '\n'
        return string
    
    def setAllFields(self, classifyKids, balance=True, **kwargs):
        self.setTrainingData(balance=balance, **kwargs)
        self.setClassifyData(classifyKids)
        self.scaleMatrices()
    
    def setTrainingData(self, balance = True, **kwargs):
        normalKids      = None
        variableKids    = None
        onePlanetKids   = None
        multiPlanetKids = None
        binaryKids      = None
        for kw in kwargs:
            if kw == 'normalKids':
                normalKids = kwargs[kw]
            elif kw == 'variableKids':
                variableKids = kwargs[kw]
            elif kw == 'onePlanetKids':
                onePlanetKids = kwargs[kw]
            elif kw == 'multiPlanetKids':
                multiPlanetKids = kwargs[kw]
            elif kw == 'binaryKids':
                binaryKids = kwargs[kw]
        basedir       = "/astro/users/johnm26/kepArchive/allKidsOfTypeX/10_03_12Lists/kidListsFromCatalogs/testLearningLists_10_10_12"
        self.ensembles = []
        self.ensembles.append(QatsEnsemble("Normal", \
            os.path.join(basedir, "AverageObjects.txt"), \
            kids=normalKids))
        self.ensembles.append(QatsEnsemble("Variable", \
            os.path.join(basedir, "VariableObjects.txt"), \
            kids=variableKids))
        self.ensembles.append(QatsEnsemble("One Planet", \
            os.path.join(basedir, "SinglePlanets.txt"), \
            kids=onePlanetKids))
        self.ensembles.append(QatsEnsemble("Multi Planet", \
            os.path.join(basedir, "MultiPlanets.txt"), \
            kids=multiPlanetKids))
        self.ensembles.append(QatsEnsemble("EclBin", \
            os.path.join(basedir, "EBList.txt"), \
            kids=binaryKids))
        #self.normalClass.readMatrix()
        #self.variablClass.readMatrix()
        #self.onePlanetClass.readMatrix()
        #self.multiPlanetClass.readMatrix()
        #self.eclBinClass.readMatrix()
        for ensemble in self.ensembles:
            print ensemble.name
            ensemble.readMatrix()
        if balance:
            self.balanceEnsembles()

        self.extractEnsembleData()
    
    def balanceEnsembles(self):
        min_nObj = len(self.ensembles[0].kidsFound)
        for ensemble in self.ensembles[1:]:
            min_nObj = min(min_nObj, len(ensemble.kidsFound))
        for ensemble in self.ensembles:
            ensemble.matrix = ensemble.matrix[:min_nObj]
            ensemble.kidsFound = ensemble.kidsFound[:min_nObj]
    
    def extractEnsembleData(self):
        self.trainKids   = num.array([])
        self.trainLabels = num.array([])
        for label_num, ensemble in enumerate(self.ensembles):
            self.trainKids   = num.hstack( \
                (self.trainKids, ensemble.kidsFound) )
            self.trainLabels = num.hstack( \
                (self.trainLabels, label_num * num.ones(len(ensemble.kidsFound))) )
        self.train_M = self.ensembles[0].matrix.copy()
        for ensemble in self.ensembles[1:]:
            self.train_M = num.vstack( \
                (self.train_M, ensemble.matrix) )
    
    def setClassifyData(self, kids):
        self.M, self.classifyKids = readInputs(kids)

    # scale test and train matrices in terms of the
    # absolute largest values out of the two
    def scaleMatrices(self):
        # merge train and test matrices
        # use positive values to ensure negative numbers
        #   are between -1 and 0
        m = num.abs(num.vstack((self.train_M, self.M)))
        # the max values from each column of m
        maxes = num.max(num.rot90(m, k=3), axis=1)
        scaleMatrix(self.train_M, maxes)
        scaleMatrix(self.M, maxes)
    
    # scale test and train matrices in terms of the
    # largest values of the training matrix
    def scaleMatrices2(self):
        maxes = num.max(num.rot90(self.train_M, k=3), axis=1)
        scaleMatrix(self.train_M, maxes)
        scaleMatrix(self.M, maxes)
    
    # scale test and train matrices from -1 to 1
    def scaleMatrices3(self):
        m     = num.vstack((self.train_M, self.M))
        maxes = num.max(num.rot90(m, k=3), axis=1)
        mins  = num.min(num.rot90(m, k=3), axis=1)
        for i in range(len(m[0])):
            a = maxes[i]
            b = mins[i]
            # linearly scales values in each column from -1 to 1
            f = lambda x: (2*x - (a + b)) / (a - b)
            m[:, i] = f(m[:, i])
            #print m[:, i].max(), m[:, i].min()
            self.train_M = m[:self.train_M.shape[0]]
            self.M       = m[-self.M.shape[0]:]

#def scaleMatrix(m):
    #for i in range(len(m[0])):
        #m[:, i] /= m[:, i].max()

# scales columns of m by values in scale
# note m.shape[1] must = len(scale)
def scaleMatrix(m, scale):
    m /= scale

# helper method for getQatsFitFile
def getSkyGroupStr(kid):
    try:
        sg = kep.dbinfo.getSkyGroup(kid)
        return 'SG' + str(sg).zfill(3)
    except:
        #print "Wanring %d failed to find sky group" % kid
        return None

# helper method for readInputs
def getQatsFitPath(kid):
    fitPath = '/astro/store/student-scratch1/johnm26/fitQatsRuns/testPool'
    fitFile = '%dQatsFit.txt' % kid
    checkPath = os.path.join(fitPath, fitFile)
    if not os.path.isfile(checkPath):
        sgStr = getSkyGroupStr(kid)
        fitPath = '/astro/store/student-scratch1/johnm26/fitQatsRuns/%s' % sgStr
        checkPath = os.path.join(fitPath, fitFile)
        if not os.path.isfile(checkPath):
            return None
    return fitPath

def readInputs(kids, inputRoot = '/astro/store/student-scratch1/johnm26/SPRING_BREAK_RUNS'):
    kidsFound = []
    qatsFound = []
    for kid in kids:
        #try:
            #SG       = kep.dbinfo.getSkyGroup(kid)
        #except:
            #print "# WARNING", kid, "failed to find sky group"
            #continue
        #fitPath = \
            #'/astro/store/student-scratch1/johnm26/fitQatsRuns/testPool'
        #sgStr    = 'SG' + str(SG).zfill(3)
        ##filename = "signal.%s.unflipped.%d.data" % (sgStr, kid)
        ##infile   = os.path.join(inputRoot, sgStr, filename)
        #infile = getdFileName(kid)
        #if not os.path.isfile(infile):
            #continue
        #fitFile = '%dQatsFit.txt' % kid
        ## new path contains non-trainer KIDs
        #if not os.path.isfile(os.path.join(fitPath, fitFile)):
            ##print os.path.join(fitPath, fitFile)
            #fitPath = '/astro/store/student-scratch1/johnm26/fitQatsRuns/%s' % sgStr
            #if not os.path.isfile(os.path.join(fitPath, fitFile)):
        fitPath = getQatsFitPath(kid)
        if fitPath == None:
            print '# Warning %s failed to find qats fit file' % kid
            continue
        else:
            consecFails = 0
            qatsdFileName = getdFileName(kid)
            if not os.path.isfile(qatsdFileName):
                continue
            dfile = open(qatsdFileName, 'r')
            periods, snr, snrLC, snrFlat = getQatsData(dfile)
            qm = QatsFeaturesModel(kid, periods, snr)
            qm.fromFile(path=fitPath)
            dfile.close()
            
            kidsFound.append(kid)
            qatsFound.append(qm.returnFeatures(n_cmax=7, n_snrmax=5, n_dchi=5))

    inputMatrix = num.zeros((len(kidsFound), len(qatsFound[0])))
    for i in range(len(kidsFound)):
        inputMatrix[i] = qatsFound[i]

    return inputMatrix, num.array(kidsFound)

#def readInputsFast(kids, directory)

