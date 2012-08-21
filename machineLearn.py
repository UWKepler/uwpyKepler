import pylab
import numpy as num
import sklearn
from sklearn.ensemble import RandomForestClassifier




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
    testConfusionMatrix()
