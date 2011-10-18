#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import optparse
import traceback

def testPipeline(KIDFile,ParFile, **kwargs):
    """
    Tests the pipeline on a list of KIDs
    and parameter options
    """
    
    par = kep.pipelinepars.PipelinePars(ParFile)
    KIDList = kep.iofiles.readKIDsFromFile(KIDFile)
    
    FileTag = par.parsetName+'.'+KIDList['File']
    PassLog = open(FileTag+'.Pass.log','w')
    FailLog = open(FileTag+'.Fail.log','w')
    BugLogFileName = FileTag+'.Bugs.log'
    BugLog = open(BugLogFileName,'w')
    SummaryLog = open(FileTag+'.summary.log','w')
    
    passCount = 0
    failCount = 0
    
    for trial in par.Trials(KIDList):
        tr = kep.pipelinepars.keptrial(trial)
        kw = kep.keplc.kw(**tr.kw)
        logString = tr.kid+' | '+kw.printString
        try:
            X = kep.keplc.keplc(tr.kid)
            X.runPipeline(kw)
            print >> PassLog, logString
            if kwargs['silent']: print 'Pass ',logString
            passCount += 1
        except:
            failCount += 1
            if kwargs['silent']: print 'Fail ', logString
            print >> FailLog, logString
            traceback.print_exc(file=open(BugLogFileName,"a"))
            print >> BugLog, '# '+logString

    print >> SummaryLog, '# passCount, failCount, totalTrials, NKIDs, Nkwopts'
    print >> SummaryLog, passCount, failCount, par.NTrials, par.NKIDs, par.Nkw
    SummaryLog.close()

if __name__ == '__main__':

    parser = optparse.OptionParser(usage=\
    "%prog -k [KIDListFile] -p [ParameterFile] (optionals->) -q")
        
    #mandatory options
    parser.add_option('-k','--kidfile',\
                        type=str,\
                        dest='kidfile',\
                        help='A file with list of KIDs')
    parser.add_option('-p','--paramfile',\
                        type=str,\
                        dest='pfile',\
                        help='The file with pipeline parameters')
                        
    parser.add_option('-q','--quiet',\
                        action="store_false",\
                        dest='verbose',\
                        help='suppress output',\
                        default=True)
    
    (opts,args) = parser.parse_args()

    mandatory_options = ['kidfile','pfile']
    Except = False
    for m in mandatory_options:
        if not opts.__dict__[m]:
            print "mandatory option %s is missing" % (m)
            Except = True
    if Except:
        parser.print_help()
        raise NameError('Mandatory options missing')

    testPipeline(opts.kidfile,opts.pfile, silent=opts.verbose)
