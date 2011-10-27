#!/astro/apps/pkg/python64/bin/python

import time, os
def readFunc():
    """ Lists all functions """
    
    HEADER = \
    " uwpyKepler used to access the local public Kepler lightcurves \n \n"+\
    " Authors - E. Agol, A. C. Becker, N. H. Thomas, C. J. Martin, \n"+\
    "           T. Steckler, J. Mehlhaff & P. Kundurthy \n \n"+\
    " Institution - University of Washington, Seattle \n"+\
    " Last Updated: "+\
    time.strftime("%Y-%b-%d %I:%M:%S %p", time.localtime())+"\n"
    
    OutReadMe = open('README','w')
    OutString = ''
    
    for file in os.listdir('.'):
        if file.endswith('.py') and not file.startswith('__init__'):
            OutString += "\n+ "+file+"\n"
            fileObj = open(file,'r')
            for line in fileObj.readlines():
                if line.startswith('def'):
                    OutString += "    "+line[3:]
                elif line.startswith('class'):
                    OutString += "     "+line
    print >> OutReadMe, HEADER
    print >> OutReadMe, OutString
    OutReadMe.close()

if __name__ == '__main__':

    readFunc()
