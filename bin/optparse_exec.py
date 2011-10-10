#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import optparse
import pylab

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage=\
    "%prog KeplerID [mandatory arguments] [optional arguments]")
        
    #mandatory options
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
                        help='The gap-size used to split'\
                        +' split portions in index units')
    parser.add_option('-w','--owin',\
                        type=int,\
                        dest='owin',\
                        help='The window-size for'\
                        +' outlier rejection')
    parser.add_option('-t','--othresh',\
                        type=float,\
                        dest='othresh',\
                        help='The threshold in sigmas for'\
                        +' outlier rejection')
    parser.add_option('-d','--dwin',\
                        type=float,\
                        dest='dwin',\
                        help='The window-size for'\
                        +' detrending')
    parser.add_option('-p','--polyorder',\
                        type=int,\
                        dest='polyorder',\
                        help='The polynomial order for'\
                        +' detrending')

    #optional arguments
    parser.add_option('-a','--agap',\
                    type=int,\
                    dest='agap',\
                    help='The size of the artifical gap for'\
                    +' bad Kepler events')
    parser.add_option('-f','--durationfactor',\
                    type=float,\
                    dest='dfac',\
                    help='Correction multiplier for'\
                    +' bad Kepler transit duration') 
    (opts,args) = parser.parse_args()

    mandatory_options = \
    ['kid','ctype','gsize','owin','othresh','dwin','polyorder']
    Except = False
    for m in mandatory_options:
        if not opts.__dict__[m]:
            print "mandatory option %s is missing" % (m)
            Except = True
    if Except:
        parser.print_help()
        exit(-1)
    
    x = kep.keplc.keplc(opts.kid) 
    kw = kep.keplc.kw(ctype=opts.ctype,gsize=opts.gsize,owin=opts.owin,\
                      othresh=opts.othresh,dwin=opts.dwin,\
                      polyorder=opts.polyorder)
    x.runPipeline(kw)
    #print x.lcFinal.keys()
    pylab.plot(x.lcFinal['x'],x.lcFinal['y'],'b.')
    pylab.plot(x.lcFinal['x'],x.lcFinal['correction'],'r-')
    pylab.show()
    pylab.plot(x.lcFinal['x'],x.lcFinal['ydt'],'b.')
    pylab.show()
    