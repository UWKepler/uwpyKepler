import numpy as num
import time

def linLeastSquares(x, y, funcs, nterms):
    # create an array of arrays; each equal in length to len(x)
    col = num.array([num.zeros(len(x))]*nterms)
    # apply functions to data to determine columns
    for i in range(len(col)):
        col[i] = funcs[i](x)
    # expand columns another dimension to attain rows
    row = num.array([col]*nterms)
    a = num.copy(row)
    b = num.zeros(nterms)
    for i in range(len(col)):
        a[i] = col[i] * row[i]
    for k in range(len(b)):
        b[k] = num.sum( (col[k] * y) )
    a = num.sum(a, axis=2)
    return num.linalg.solve(a, b), col[-1]