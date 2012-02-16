#!/astro/apps/pkg/python64/bin//python
#Copies tar files from the data archive online to
#'destination' (default Ben's scratch folder).
#Big files = run this in a scratch folder.
#you supply a list of quarters
#
import sys, os

root = 'public_'
exten = '.tgz'
path = 'http://archive.stsci.edu/pub/kepler/lightcurves/tarfiles/'
destination = '/astro/store/student-scratch1/bvegaff'
if __name__ == '__main__':
     
    
    arglen = len(sys.argv)
    for i in range(arglen-1):
        qnum = sys.argv[i+1] 
        tag = 'Q'+qnum
        url = path+root+tag+exten
        print 'wget %s' % url
        os.system('wget -P %s %s ' % (destination, url))
