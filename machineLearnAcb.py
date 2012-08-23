import os
import pylab
import numpy as num
import sklearn
import uwpyKepler as kep
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
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
        
        #import pylab
        #pylab.plot(periods, snr, label="snr")
        #pylab.plot(periods, snrLC, label="LC")
        #pylab.plot(periods, snrFLAT, label="FLAT")
        #pylab.legend()
        #pylab.show()

        kidsFound.append(kid)
        qatsFound.append(snr)

    inputMatrix = num.zeros((len(kidsFound), len(qatsFound[0])))
    for i in range(len(kidsFound)):
        inputMatrix[i] = qatsFound[i]

    return inputMatrix


def plotConfusionMatrix(truth, predicted, labels, names, title = "Confusion Matrix"):
    matrix = sklearn.metrics.confusion_matrix(truth, predicted, labels)
    fig = pylab.figure()
    ax  = fig.add_subplot(111)
    im  = ax.matshow(matrix, cmap = pylab.cm.OrRd)
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            ax.text(x, y, "%d" % (matrix[y][x]))

    pylab.xticks(labels, names)
    pylab.yticks(labels, names)
    pylab.ylabel("Truth")
    #pylab.colorbar(im)
    pylab.title(title)
    #pylab.show()

def testConfusionMatrix():
    truth      = num.random.randint(0, 3, 500)
    predicted  = num.random.randint(0, 3, 500)
    labels     = [0, 1, 2]
    names      = ["A", "B", "C"]
    plotConfusionMatrix(truth, predicted, labels, names)

if __name__ == '__main__':
    if 1:
        basedir       = "/astro/store/student-scratch1/martincj/condor/SpringBreakRuns/Scratch/TrainingWheels/"

        normalClass      = QatsEnsemble("Normal", os.path.join(basedir, "AverageObjects"))
        variableClass    = QatsEnsemble("Variable", os.path.join(basedir, "VariableObjects"))
        onePlanetClass   = QatsEnsemble("One Planet", os.path.join(basedir, "SinglePlanets.txt"))
        multiPlanetClass = QatsEnsemble("Multi Planet", os.path.join(basedir, "MultiPlanets.txt"))
        eclBinClass      = QatsEnsemble("EclBin", os.path.join(basedir, "EBList.txt"))

        normalClass.readMatrix()
        variableClass.readMatrix()
        onePlanetClass.readMatrix()
        multiPlanetClass.readMatrix()
        eclBinClass.readMatrix()

        ntovalidate = 10

        # Merge matrices into a single matrix to classifier
        trainingSample = normalClass.matrix[:-ntovalidate].copy()
        trainingSample = num.concatenate((trainingSample, variableClass.matrix[:-ntovalidate]))
        trainingSample = num.concatenate((trainingSample, onePlanetClass.matrix[:-ntovalidate]))
        trainingSample = num.concatenate((trainingSample, multiPlanetClass.matrix[:-ntovalidate]))
        trainingSample = num.concatenate((trainingSample, eclBinClass.matrix[:-ntovalidate]))

        classLabels    = num.zeros(normalClass.matrix[:-ntovalidate].shape[0])
        classLabels    = num.concatenate((classLabels, 1 * num.ones(variableClass.matrix[:-ntovalidate].shape[0])))
        classLabels    = num.concatenate((classLabels, 2 * num.ones(onePlanetClass.matrix[:-ntovalidate].shape[0])))
        classLabels    = num.concatenate((classLabels, 3 * num.ones(multiPlanetClass.matrix[:-ntovalidate].shape[0])))
        classLabels    = num.concatenate((classLabels, 4 * num.ones(eclBinClass.matrix[:-ntovalidate].shape[0])))

        # Do the machine learning!
        decisionTree   = DecisionTreeClassifier(max_depth=None, min_samples_split=1, random_state=0).fit(trainingSample, classLabels)
        randomForest   = RandomForestClassifier(n_estimators=50, max_depth=None, min_samples_split=1, random_state=0).fit(trainingSample, classLabels)
        extraTrees     = ExtraTreesClassifier(n_estimators=50,   max_depth=None, min_samples_split=1, random_state=0).fit(trainingSample, classLabels)
        gradientBoost  = GradientBoostingClassifier(n_estimators=50, max_depth=1, learn_rate=1.0, random_state=0).fit(trainingSample, classLabels)

        #ntovalidate    = 50
        truthLabels    = num.zeros(ntovalidate)
        truthLabels    = num.concatenate((truthLabels, 1 * num.ones(ntovalidate)))
        truthLabels    = num.concatenate((truthLabels, 2 * num.ones(ntovalidate)))
        truthLabels    = num.concatenate((truthLabels, 3 * num.ones(ntovalidate)))
        truthLabels    = num.concatenate((truthLabels, 4 * num.ones(ntovalidate)))

        names          = [x.name for x in (normalClass,variableClass,onePlanetClass,multiPlanetClass,eclBinClass)]

        for classifier, name in ((decisionTree, "Decision Tree"),
                                 (randomForest, "Random Forest"),
                                 (extraTrees, "Extra Trees"),
                                 (gradientBoost, "Gradient Boosting")):
            
            print "SCORE", name, cross_val_score(classifier, trainingSample, classLabels)

            predictions   = classifier.predict(normalClass.matrix[-ntovalidate:])
            predictions   = num.concatenate((predictions, classifier.predict(variableClass.matrix[-ntovalidate:])))
            predictions   = num.concatenate((predictions, classifier.predict(onePlanetClass.matrix[-ntovalidate:])))
            predictions   = num.concatenate((predictions, classifier.predict(multiPlanetClass.matrix[-ntovalidate:])))
            predictions   = num.concatenate((predictions, classifier.predict(eclBinClass.matrix[-ntovalidate:])))
        
            plotConfusionMatrix(truthLabels, predictions, num.arange(5), names, name)

        pylab.show()
        
        import pdb; pdb.set_trace()

        kids      = num.loadtxt(infile, unpack = True)
        print kids
        readInputs(kids)
    if 0:
        testConfusionMatrix()
