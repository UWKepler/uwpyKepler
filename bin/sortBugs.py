#!/astro/apps/pkg/python64/bin//python

import sys
import uwpyKepler as kep

if __name__ == '__main__':
    
    logFile = sys.argv[1]
    kep.debug.getUniqueTraceback(logFile)