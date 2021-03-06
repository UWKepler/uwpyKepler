import sys
import MySQLdb
import numpy as num
#import pylab
#num.warning.warn(action = 'ignore')
            
def ApplyMask(lcdtData,mask):
    """ This function applies a given mask """
    
    # match data keys and apply mask    
    for key in lcdtData.keys():
        if key in 'xyerr':
            if mask != 'NoMask':
                lcdtData[key].mask = lcdtData['NoMask']
            lcdtData[key].mask = lcdtData[mask]
            
    return lcdtData

def ApplyMaskPortions(lcData,mask,portion):
    """ This function applies a given mask """
    
    # match data keys and apply mask    
    for key in lcData[portion].keys():
        if key in 'xyerr':
            if mask != 'NoMask':
                lcData[portion][key].mask = lcData[portion]['NoMask']
            lcData[portion][key].mask = lcData[portion][mask]
            
    return lcData

def returnData(lcdtData,typeTag):
    """
        typeTags
        all - All data (return all data without any masks)
        eonly - eclipse portions only (return only eclipse portions)
        elc - eclipse and flat portions (return data with outliers (o & k) removed)
        flat - only the flat part of the lightcurve (return data with o, k and eclipses removed)
        allflags - all flagged points, o, k and eclipse (opposite of flat)
        o - uwpyKepler outliers only (return outliers)
        k - Kepler flagged events (return Kepler flagged outliers)
        ok - uwpyKepler and Kepler outliers (return all outliers)
        noeclipse - return everything but the eclipse
    """

    lcdtData = ApplyMask(lcdtData,'NoMask')
    
    if typeTag.lower() == 'elc':
        idx = num.where( (lcdtData['OKMask'] == False))[0]
    elif typeTag.lower() == 'o':
        idx = num.where((lcdtData['OMask']))[0]
    elif typeTag.lower() == 'k':
        idx = num.where((lcdtData['KEMask']))[0]
    elif typeTag.lower() == 'ok':
        idx = num.where((lcdtData['OKMask']))[0]
    elif typeTag.lower() == 'eonly':
        idx = num.where((lcdtData['eMask']))[0]
    elif typeTag.lower() == 'flat' :
        idx = num.where((lcdtData['ALLMask'] == False))[0]
    elif typeTag.lower() == 'allflags' :
        idx = num.where((lcdtData['ALLMask']))[0]
    elif typeTag.lower() == 'all':
        idx = num.where((lcdtData['NoMask'] == False))[0]
    elif typeTag.lower() == 'noeclipse':
        idx0 = num.where((lcdtData['ALLMask'] == False))[0]
        idx1 = num.where((lcdtData['OMask']))[0]
        idx = list(set(num.hstack( (idx0,idx1) ) ) )
    else:
        'No match for typeTag found'
        idx = []
    
    x = lcdtData['x'][idx]
    y = lcdtData['ydt'][idx]
    yerr = lcdtData['yerrdt'][idx]
    cad = lcdtData['cadence'][idx]
    
    return x,y,yerr,cad