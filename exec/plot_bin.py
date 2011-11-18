import sys
import uwpyKepler as kep
import numpy as num
import pylab
import optparse

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage=\
    "%prog KeplerID [mandatory arguments] [optional arguments]")
        
    #mandatory options (none mandatory but kid for this script)
    parser.add_option('-k','--kid',\
                        type=long,\
                        dest='kid',\
                        help='The Kepler ID')
    cChoice = ('LC','SC','')
    parser.add_option('-c','--ctype',\
                        choices=cChoice,\
                        type='choice',\
                        dest='ctype',\
                        default='LC',\
                        help='The expected cadence type'\
                        +', '.join(cChoice)\
                        +') [default: %default]')
    parser.add_option('-g','--gsize',\
                        type=int,\
                        dest='gsize',\
                        default=2,\
                        help='The gap-size used to split'\
                        +' split portions in index units')
    parser.add_option('-w','--owin',\
                        type=int,\
                        dest='owin',\
                        default=15,\
                        help='The window-size for'\
                        +' outlier rejection')
    parser.add_option('-t','--othresh',\
                        type=float,\
                        dest='othresh',\
                        default=5,\
                        help='The threshold in sigmas for'\
                        +' outlier rejection')
    parser.add_option('-d','--dwin',\
                        type=float,\
                        dest='dwin',\
                        default=150,\
                        help='The window-size for'\
                        +' detrending')
    parser.add_option('-p','--polyorder',\
                        type=int,\
                        dest='polyorder',\
                        default=6,\
                        help='The polynomial order for'\
                        +' detrending')

    #optional arguments
    parser.add_option('-a','--agap',\
                    type=int,\
                    dest='agap',\
                    default=1,\
                    help='The size of the artifical gap for'\
                    +' bad Kepler events')
    parser.add_option('-f','--durationfactor',\
                    type=float,\
                    dest='durfac',\
                    default=2,\
                    help='Correction multiplier for'\
                    +' bad Kepler transit duration') 
    
    #unique options for this binning script
    parser.add_option('-s','--sharex',\
                    type=int,\
                    dest='sharex',\
                    default=0,\
                    help='0 to display bins on same plot;'\
                    +' 1 to divide bins into subplots')
    parser.add_option('-o','--fileoptions',\
                    type=str,\
                    dest='fileoptions',\
                    default=None,\
                    help='Specify a parset file to read'\
                    +' in its parameters')

    (opts,args) = parser.parse_args()

    mandatory_options = \
    ['kid']
    Except = False
    for m in mandatory_options:
        if not opts.__dict__[m]:
            print "mandatory option %s is missing" % (m)
            Except = True
    if Except:
        parser.print_help()
        raise NameError('Mandatory options missing')

if opts.fileoptions != None:
    parlist = []
    file = open(opts.fileoptions, 'r')
    lines = file.readlines()
    for line in lines:
        line = str(line)
        line = line.split('|')
        parlist.append(line[1])
    opts.ctype     = parlist[2]
    opts.gsize     = int(parlist[3])
    opts.owin      = int(parlist[4])
    opts.othresh   = int(parlist[5])
    opts.dwin      = int(parlist[6])
    opts.polyorder = int(parlist[7])
    opts.agap      = int(parlist[8])
    opts.durfac    = int(parlist[9])

keplc = kep.keplc.keplc(opts.kid)
eData = keplc.eData
BJDREFI = keplc.BJDREFI
kw = kep.keplc.kw(ctype=opts.ctype,gsize=opts.gsize,owin=opts.owin,\
                      othresh=opts.othresh,dwin=opts.dwin,\
                      polyorder=opts.polyorder,\
                      agap=opts.agap,durfac=opts.durfac)
lcData = kep.keplc.lcData(opts.kid,eData,BJDREFI, kw).lcData
lst = [1, 5, 10, 20, 50, 100]
bind = kep.func.bin(lcData, lst)

colors = ['r','y','b','c','g','m','k']
if opts.sharex == 0:
    for binidx in range(len(bind.keys())):
        pylab.plot(bind[bind.keys()[binidx]]['x'], bind[bind.keys()[binidx]]['y'],colors[binidx]+'o')
elif opts.sharex == 1:
    nplots = len(lst)
    for binidx in range(len(bind.keys())):
        if binidx == 0:
            ax1 = pylab.subplot(nplots,1,binidx+1)
            pylab.plot(bind[bind.keys()[binidx]]['x'], bind[bind.keys()[binidx]]['y'],colors[binidx]+'o')
        else:
            pylab.subplot(nplots,1,binidx+1,sharex=ax1)
            pylab.plot(bind[bind.keys()[binidx]]['x'], bind[bind.keys()[binidx]]['y'],colors[binidx]+'o')
        pylab.title('Binsize: %s' % lst[binidx])
                
pylab.show()

