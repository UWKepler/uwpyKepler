import numpy as num
import scipy

# f_i_0 = 1
# f_i_1 = x
# f_i_2 = x**2
# f_i_3 = x**3
# f_i_4 = Q(x)
# 

def linLeastSquares(x, y, funcs, nterms):
    a = num.zeros(nterms, nterms)
    b = num.zeros(nterms)
    for k in range(len(b)):
        if isinstance(funcs[k], int):
            b[k] = (num.power(x, funcs[k]) * y).sum()
        else:
            b[k] = (funcs[k](x) * y).sum()
    for i in range(len(a)):
        for j in range(len(a[i])):
            if isinstance(funcs[i], int)\
               and isinstance(funcs[j], int):
                a[i][j] = (num.power(x, funcs[i])\
                         * num.power(x, funcs[j])).sum()
            elif isinstance(funcs[i], int)\
            and not isinstance(funcs[j], int):
                a[i][j] = (num.power(x, funcs[i])\
                         * funcs[j](x)).sum()
            elif not isinstance(funcs[i], int)\
            and isinstance(funcs[j], int):
                a[i][j] = (funcs[i](x)\
                         * num.power(x, funcs[j])).sum()
            else:
                a[i][j] = (funcs[i](x) * funcs[j](x)).sum()
    return scipy.linalg.solve(a, b)
    