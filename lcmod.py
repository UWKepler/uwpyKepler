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
        all - All data (unmasked)
        eonly - eclipse portions only
        elc - eclipse and flat portions (outliers removed)
        flat - only the flat part of the lightcurve
        allflags - all flagged points, o, k and eclipse (opposite of flat)
        o - uwpyKepler outliers only
        k - Kepler flagged events (aka Kepler outliers)
        ok - uwpyKepler and Kepler outliers
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
    else:
        'No match for typeTag found'
        idx = []
    
    x = lcdtData['x'][idx]
    y = lcdtData['ydt'][idx]
    yerr = lcdtData['yerrdt'][idx]
    cad = lcdtData['cadence'][idx]
    
    return x,y,yerr,cad