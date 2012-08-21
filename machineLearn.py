import pylab
import numpy as num
import sklearn
from sklearn.ensemble import RandomForestClassifier


def plotConfusionMatrix(truth, predicted, labels):
    matrix = sklearn.metrics.confusion_matrix(truth, predicted, labels)
    fig = pylab.figure()
    ax  = fig.add_subplot(111)
    ax.imshow(matrix)
    pylab.show()

def testConfusionMatrix():
    truth      = num.random.randint(0, 3, 500)
    predicted  = num.random.randint(0, 3, 500)
    labels     = [0, 1, 2]
    plotConfusionMatrix(truth, predicted, labels)

if __name__ == '__main__':
    testConfusionMatrix()
