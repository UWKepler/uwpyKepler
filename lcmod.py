import sys
import MySQLdb
import numpy as num
import pylab
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
        typeTag - all, eonly, elc, o, k, ok
    """

    lcdtData = ApplyMask(lcdtData,'NoMask')
    #print len(num.where(lcdtData['OKMask'])[0])
    #print len(num.where(lcdtData['OMask'])[0])
    #print len(num.where(lcdtData['eMask'])[0])
    #print len(num.where(lcdtData['KEMask'])[0])
    #print len(num.where(lcdtData['NoMask'])[0])
    
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
        idx = num.where( (lcdtData['OKMask'] == False) & (lcdtData['eMask'] == False) )[0]
    elif typeTag.lower() == 'all':
        idx = num.where((lcdtData['NoMask'] == False))[0]
    else:
        'No match for typeTag found'
        idx = []
    
    print len(idx), typeTag
    x = lcdtData['x'][idx]
    y = lcdtData['ydt'][idx]
    yerr = lcdtData['yerrdt'][idx]
    
    return x,y,yerr