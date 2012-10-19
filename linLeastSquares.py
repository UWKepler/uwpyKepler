import numpy as num

# pre  : func(x) returns array of same length for all func in funcs
# pre  : elements of y correspond to elements of funcs(x)
# post : returns best fit coefficients for fitting x to y with funcs
def linLeastSquares(x, y, funcs, nterms, **kwargs):
    # create an array of arrays; 
    # each equal in length to len(func[0](x))
    firstEval = funcs[0](x)
    col = num.zeros( (nterms, len(firstEval)) )
    col[0] = firstEval
    # apply functions to data to determine columns
    for i in range(len(col))[1:]:
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
    # option to return values allows elimination of
    # time-consuming computations later
    #pdb.set_trace()
    for key in kwargs:
        if key == 'return_func_vals':
            returnFuncs = kwargs[key]
            return num.linalg.solve(a, b), col[returnFuncs]

    return num.linalg.solve(a, b)

# in some cases it is best for speed if
# the functions (see "funcs" argument above)
# are evaluated before this function is called
def linLeastSquares_funcsEvaluated(x, y, evaluatedFuncsOf_x, nterms):
    fx = evaluatedFuncsOf_x
    #col = num.array([num.zeros(len(x))]*nterms)
    col = num.zeros( (nterms, len(x)) )
    for i in range(len(col)):
        col[i] = fx[i]
    row = num.array([col]*nterms)
    a = num.copy(row)
    b = num.zeros(nterms)
    for i in range(len(col)):
        a[i] = col[i] * row[i]
    for k in range(len(b)):
        b[k] = num.sum( (col[k] * y) )
    a = num.sum(a, axis=2)

    return num.linalg.solve(a, b)