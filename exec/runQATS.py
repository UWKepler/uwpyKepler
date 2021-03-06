#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import optparse

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
                        default=2,\
                        help='The gap-size used to split'\
                        +' split portions in index units')
    parser.add_option('-w','--owin',\
                        type=int,\
                        dest='owin',\
                        default=15,\    # note, this suits only LC
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
                        default=50,\   # note, this suits only LC
                        help='The window-size for'\
                        +' detrending')
    parser.add_option('-p','--polyorder',\
                        type=int,\
                        dest='polyorder',\
                        default=5,\
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
    
    # check agap size
    if opts.gize > opts.agap:
        opts.agap = opts.gize

    x = kep.keplc.keplc(opts.kid)
    kw = kep.keplc.kw(ctype=opts.ctype,gsize=opts.gsize,owin=opts.owin,\
                      othresh=opts.othresh,dwin=opts.dwin,\
                      polyorder=opts.polyorder,\
                      agap=opts.agap,durfac=opts.durfac)

    x.runPipeline(kw)