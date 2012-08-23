import os
import pylab
import numpy as num
import sklearn
import uwpyKepler as kep
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.cross_validation import cross_val_score

class QatsEnsemble(object):
    def __init__(self, name, infile, kids = None, matrix = None):
        self.name   = name
        self.infile = infile
        self.kids   = kids
        self.matrix = matrix

    def readKids(self):
        self.kids = num.loadtxt(self.infile, unpack = True)

    def readMatrix(self):
        if self.kids == None:
            self.readKids()
        self.matrix = readInputs(self.kids)

def readInputs(kids, inputRoot = "/astro/store/student-scratch1/johnm26/SPRING_BREAK_RUNS"):
    kidsFound = []
    qatsFound = []
    for kid in kids:
        try:
            SG       = kep.dbinfo.getSkyGroup(kid)
        except:
            print "# WARNING", kid, "failed to find sky group"
            continue

        sgStr    = 'SG' + str(SG).zfill(3)
        filename = "signal.%s.unflipped.%d.data" % (sgStr, kid)
        infile   = os.path.join(inputRoot, sgStr, filename)
        if not os.path.isfile(infile):
            continue
        periods, snr, snrLC, snrFLAT = kep.postqats.getQatsData(open(infile))

        kidsFound.append(kid)
        qatsFound.append(snr)

    inputMatrix = num.zeros((len(kidsFound), len(qatsFound[0])))
    for i in range(len(kidsFound)):
        inputMatrix[i] = qatsFound[i]

    return inputMatrix


def plotConfusionMatrix(truth, predicted, labels, names):
    matrix = sklearn.metrics.confusion_matrix(truth, predicted, labels)
    fig = pylab.figure()
    ax  = fig.add_subplot(111)
    im  = ax.imshow(matrix, interpolation='nearest', cmap = pylab.cm.binary)
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            ax.text(x, y, "%d" % (matrix[y][x]))

    pylab.xticks(labels, names)
    pylab.yticks(labels, names)
    pylab.colorbar(im)
    pylab.show()

def testConfusionMatrix():
    truth      = num.random.randint(0, 3, 500)
    predicted  = num.random.randint(0, 3, 500)
    labels     = [0, 1, 2]
    names      = ["A", "B", "C"]
    plotConfusionMatrix(truth, predicted, labels, names)

if __name__ == '__main__':
    if 1:
        basedir       = "/astro/store/student-scratch1/martincj/condor/SpringBreakRuns/Scratch/TrainingWheels/"
        
        #Getting the lists and creating objects
        normalClass   = QatsEnsemble("Normal", os.path.join(basedir, "AverageObjects"))
        variableClass = QatsEnsemble("Variable", os.path.join(basedir, "VariableObjects"))
        testClass     = QatsEnsemble("Test", os.path.join(basedir, "AverageTestObjects"))

        #creating matrices
        normalClass.readMatrix()
        variableClass.readMatrix()
        testClass.readMatrix()

        # Merge matrices into a single matrix to classifier
        trainingSample = normalClass.matrix.copy()
        trainingSample = num.concatenate((trainingSample, variableClass.matrix))

        classLabels    = num.zeros(normalClass.matrix.shape[0])
        classLabels    = num.concatenate((classLabels, 1 * num.ones(variableClass.matrix.shape[0])))

        #Creating Classifiers
        decisionTree   = DecisionTreeClassifier(max_depth=None, min_samples_split=1, random_state=0)
        randomForest   = RandomForestClassifier(n_estimators=10, max_depth=None, min_samples_split=1, random_state=0)
        extraTrees     = ExtraTreesClassifier(n_estimators=10, max_depth=None, min_samples_split=1, random_state=0)

        #
        print "SCORE1", cross_val_score(decisionTree, trainingSample, classLabels)
        print "SCORE2", cross_val_score(randomForest, trainingSample, classLabels)
        print "SCORE3", cross_val_score(extraTrees, trainingSample, classLabels) #, cv = sklearn.cross_validation)???

        decisionTree   = decisionTree.fit(trainingSample, classLabels)
        randomForest   = randomForest.fit(trainingSample, classLabels)
        extraTrees     = extraTrees.fit(trainingSample, classLabels)

        print "GUESS1", decisionTree.predict(testClass.matrix)
        print decisionTree.predict_proba(testClass.matrix)

        print "GUESS2", randomForest.predict(testClass.matrix)
        print randomForest.predict_proba(testClass.matrix)

        print "GUESS3", extraTrees.predict(testClass.matrix)
        print extraTrees.predict_proba(testClass.matrix)

        import pdb; pdb.set_trace()

        kids      = num.loadtxt(infile, unpack = True)
        print kids
        readInputs(kids)
    if 0:
        testConfusionMatrix()
