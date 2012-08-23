import os
import pylab
import numpy as num
import sklearn
import uwpyKepler as kep
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.cross_validation import cross_val_score, Bootstrap, KFold, LeaveOneOut
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import svm


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

def plotConfusionMatrix(matrix, title):
    names = ['Average', 'Variable', 'EB']
    fig = pylab.figure()
    ax  = fig.add_subplot(111)
    im  = ax.imshow(matrix, interpolation='nearest', cmap = pylab.cm.jet)
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            ax.text(x, y, "%d" % (matrix[y][x]))

    pylab.xticks(range(len(names)), names)
    pylab.yticks(range(len(names)), names)
    pylab.title(title)
    pylab.colorbar(im)
    pylab.show()


if __name__ == '__main__':
    if 1:
        basedir       = "/astro/store/student-scratch1/martincj/condor/SpringBreakRuns/Scratch/TrainingWheels/"

        #Getting the lists and creating objects
        averageObjectTraining   = QatsEnsemble("Average", os.path.join(basedir, "AverageObjects"))
        variableObjectTraining = QatsEnsemble("Variable", os.path.join(basedir, "VariableObjects"))
        EBTraining = QatsEnsemble("EB", os.path.join(basedir, "EBList.txt"))
        
        averageObjectTesting = QatsEnsemble("Average Test", os.path.join(basedir, "AverageTestObjects"))
        variableObjectTesting = QatsEnsemble("Average Test", os.path.join(basedir, "VariableTestingObjects"))
        EBTesting =  QatsEnsemble("EB Test", os.path.join(basedir, "EBTestList2"))


        #creating matrices
        averageObjectTraining.readMatrix()
        variableObjectTraining.readMatrix()
        EBTraining.readMatrix()
        
        averageObjectTesting.readMatrix()
        variableObjectTesting.readMatrix()
        EBTesting.readMatrix()

        # Merge matrices into a single large training set Matrix for classifiers
        trainingSet = num.concatenate((averageObjectTraining.matrix.copy(),\
                        variableObjectTraining.matrix,\
                        EBTraining.matrix))
        testingSet =  num.concatenate((averageObjectTesting.matrix.copy(),\
                        variableObjectTesting.matrix,\
                        EBTesting.matrix))
        trainingLabels = num.concatenate((num.zeros(averageObjectTraining.matrix.shape[0]),\
                        1 * num.ones(variableObjectTraining.matrix.shape[0]),\
                        2 * num.ones(EBTraining.matrix.shape[0])))
        testingLabels = num.concatenate((num.zeros(averageObjectTesting.matrix.shape[0]),\
                        1 * num.ones(variableObjectTesting.matrix.shape[0]),\
                        2 * num.ones(EBTesting.matrix.shape[0])))

        #Creating classifiers
        decisionTree   = DecisionTreeClassifier(max_depth=None, min_samples_split=1, random_state=0)
        randomForest   = RandomForestClassifier(n_estimators=10, max_depth=None, min_samples_split=1, random_state=0)
        extraTrees     = ExtraTreesClassifier(n_estimators=10, max_depth=None, min_samples_split=1, random_state=0)
        SVCrbf = svm.SVC(kernel='rbf', gamma=0.7)
        SVClinear = svm.SVC(kernel='linear')
        SVCpoly = svm.SVC(kernel='poly', degree=3)
        LinearSVC = svm.LinearSVC()
        NuSVC = svm.NuSVC()

        #Applying classifiers to data
        decisionTree   = decisionTree.fit(trainingSet, trainingLabels)
        randomForest   = randomForest.fit(trainingSet, trainingLabels)
        extraTrees     = extraTrees.fit(trainingSet, trainingLabels)
        SVCrbf = SVCrbf.fit(trainingSet, trainingLabels)
        SVClinear = SVClinear.fit(trainingSet, trainingLabels)
        SVCpoly = SVCpoly.fit(trainingSet, trainingLabels)
        LinearSVC = LinearSVC.fit(trainingSet, trainingLabels)
        NuSVC = NuSVC.fit(trainingSet, trainingLabels)

        #Finding predictions 
        decisionTreePredictions = decisionTree.predict(testingSet)
        randomForestPredictions = randomForest.predict(testingSet)
        extraTreesPredictions = extraTrees.predict(testingSet)
        SVCrbfPredictions = SVCrbf.predict(testingSet)
        SVClinearPredictions = SVClinear.predict(testingSet)
        SVCpolyPredictions = SVCpoly.predict(testingSet)
        LinearSVCPredictions = LinearSVC.predict(testingSet)
        NuSVCPredictions = NuSVC.predict(testingSet)

        #Finding probabilities of each prediction (confidence)
        decisionTreePredictionProbs = decisionTree.predict_proba(testingSet)
        randomForestPredictionProbs = randomForest.predict_proba(testingSet)
        extraTreesPredictionProbs = extraTrees.predict_proba(testingSet)



        ####################################################
        #Metrics evaluating how good our models are
        ####################################################
        #Cross validating
        decisionTreeCVS = cross_val_score(decisionTree, trainingSet, trainingLabels)
        randomForestCVS = cross_val_score(randomForest, trainingSet, trainingLabels)
        extraTreesCVS = cross_val_score(extraTrees, trainingSet, trainingLabels)
        SVCrbfCVS = cross_val_score(SVCrbf, trainingSet, trainingLabels)
        SVClinearCVS = cross_val_score(SVClinear, trainingSet, trainingLabels)
        SVCpolyCVS = cross_val_score(SVCpoly, trainingSet, trainingLabels)
        LinearSVCCVS = cross_val_score(LinearSVC, trainingSet, trainingLabels)
        NuSVCCVS = cross_val_score(NuSVC, trainingSet, trainingLabels)
        
        #decisionTreeKFCV = KFold(extraTrees, trainingSet, classLabels)
        #decisionTreeLOOCV = LeaveOneOut(extraTrees, trainingSet, classLabels)
        #decisionTreeBCV = Bootstrap(extraTrees, trainingSet, classLabels)
        
        #Classification reports
        decisionTreeClassReport = classification_report(testingLabels, decisionTreePredictions)
        randomForestClassReport = classification_report(testingLabels, randomForestPredictions)
        extraTreesClassReport = classification_report(testingLabels, extraTreesPredictions)
        SVCrbfClassReport = classification_report(testingLabels, SVCrbfPredictions)
        SVClinearClassReport = classification_report(testingLabels, SVClinearPredictions)
        SVCpolyClassReport = classification_report(testingLabels, SVCpolyPredictions)
        LinearSVCClassReport = classification_report(testingLabels, LinearSVCPredictions)
        NuSVCClassReport = classification_report(testingLabels, NuSVCPredictions)
        
        #Find Confusion Matrix
        decisionTreeConfusionMatrix = confusion_matrix(testingLabels, decisionTreePredictions)
        randomForestConfusionMatrix = confusion_matrix(testingLabels, randomForestPredictions)
        extraTreesConfusionMatrix = confusion_matrix(testingLabels, extraTreesPredictions)
        SVCrbfConfusionMatrix = confusion_matrix(testingLabels, SVCrbfPredictions)
        SVClinearConfusionMatrix = confusion_matrix(testingLabels, SVClinearPredictions)
        SVCpolyConfusionMatrix = confusion_matrix(testingLabels, SVCpolyPredictions)
        LinearSVCConfusionMatrix = confusion_matrix(testingLabels, LinearSVCPredictions)
        NuSVCConfusionMatrix = confusion_matrix(testingLabels, NuSVCPredictions)






        ########################################
        #Visualizing the learning
        ########################################
        plotConfusionMatrix(decisionTreeConfusionMatrix, 'decisionTreeConfusionMatrix')
        plotConfusionMatrix(randomForestConfusionMatrix, 'randomForestConfusionMatrix')
        plotConfusionMatrix(extraTreesConfusionMatrix, 'extraTreesConfusionMatrix')
        plotConfusionMatrix(SVCrbfConfusionMatrix, 'SVCrbfConfusionMatrix')
        plotConfusionMatrix(SVClinearConfusionMatrix, 'SVClinearConfusionMatrix')
        plotConfusionMatrix(SVCpolyConfusionMatrix, 'SVCpolyConfusionMatrix')
        plotConfusionMatrix(LinearSVCConfusionMatrix, 'LinearSVCConfusionMatrix')
        plotConfusionMatrix(NuSVCConfusionMatrix, 'NuSVCConfusionMatrix')
        
        #X = trainingSet
        #y = trainingLabels
        
        #X = X[y != 0, :2]
        #y = y[y != 0]
        
        #n_sample = len(X)
        
        #num.random.seed(0)
        #order = num.random.permutation(n_sample)
        #X = X[order]
        #y = y[order].astype(num.float)
        
        #X_train = X[:.9 * n_sample]
        #y_train = y[:.9 * n_sample]
        #X_test = X[.9 * n_sample:]
        #y_test = y[.9 * n_sample:]
        ## fit the model
        #for fig_num, kernel in enumerate(('linear', 'rbf', 'poly')):
            #clf = svm.SVC(kernel=kernel, gamma=10)

            #clf.fit(X_train, y_train)
        
            #pylab.figure(fig_num)
            #pylab.clf()
            #pylab.scatter(X[:, 0], X[:, 1], c=y, zorder=10, cmap=pylab.cm.Paired)
            #print 'eh'
            ## Circle out the test data
            #pylab.scatter(X_test[:, 0], X_test[:, 1],
                    #s=80, facecolors='none', zorder=10)
        
            #pylab.axis('tight')
            #x_min = X[:, 0].min()
            #x_max = X[:, 0].max()
            #y_min = X[:, 1].min()
            #y_max = X[:, 1].max()
            #print '2'
            #XX, YY = num.mgrid[x_min:x_max:200j, y_min:y_max:200j]
            #Z = clf.decision_function(num.c_[XX.ravel(), YY.ravel()])
        
            ## Put the result into a color plot
            #Z = Z.reshape(XX.shape)
            #pylab.pcolormesh(XX, YY, Z > 0, cmap=pylab.cm.Paired)
            #pylab.contour(XX, YY, Z, colors=['k', 'k', 'k'],
                    #linestyles=['--', '-', '--'],
                    #levels=[-.5, 0, .5])
        
            #pylab.title(kernel)
            #pylab.show()

    if 0:
        testConfusionMatrix()
