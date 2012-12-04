"""EllipseFitTst*.py tries to fit an ellipse to a galaxy shape
   by determining points where the flux drops to a certain value
   as a function of distance from the weighted center.  I'll try
   Cts < 2 for starters.  Given b/a from the catalog, the rotational
   angle psi for a galaxy can be solved analytically, as can the semi-
   major axis a. Bob Abel, 7/2/12"""

# First import a bunch of modules, including pyfits for the HDUlist 
# (FITS header/data unit list),and the afw detection and image modules.
# matchTrimToDet probably matches the Trimfile coordinates to the 
# detection coordinates.  SIP (Simple Imaging Polynomial) convention 
# represents  non-linear geometric distortion as polynomials in FITS header keywords.
# I don't think gzip is used - how are the .gz files unpacked?
import os, sys, gzip, numpy, pyfits, math
import lsst.afw.detection as afwDetection
import lsst.afw.image as afwImage
import simplePlot, matchTrimToDet, xyToSkyTANSIP
from lsst.sims.catalogs.generation.db import queryDB
import lsst.afw.coord as afwCoord
import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
from matplotlib import pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.patches import Circle


# Load the SFD extinction maps
from lsst.sims.catalogs.measures.photometry import EBV
ebvMapNorth = EBV.EbvMap()
datadir = os.environ.get("CAT_SHARE_DATA")
ebvMapNorth.readMapFits(os.path.join(datadir, "data/Dust/SFD_dust_4096_ngp.fits"))
ebvMapSouth = EBV.EbvMap()
ebvMapSouth.readMapFits(os.path.join(datadir, "data/Dust/SFD_dust_4096_sgp.fits"))

def positionAngle(x,y):
   """Given (x,y) or (ra,dec), calculates the position angle
      in degrees, with respect to the positive y-axis, with
      the angle increases in the CCW direction."""
   pi = math.pi
   radeg = 180./pi
   xabs = abs(x)
   yabs = abs(y)

   if x <= 0:
      if y >= 0:
         psi = math.atan(xabs/yabs)*radeg
      else:
         psi = (math.atan(yabs/xabs) + pi/2.)*radeg
   elif x > 0:
      if y < 0:
         psi = (math.atan(xabs/yabs) + pi)*radeg
      else:
         psi = (math.atan(yabs/xabs) + 3.*pi/2.)*radeg
   else:
      print "positionAngle inputs invalid"
      psi = numpy.nan
   return psi


def xyPosAngle(x, y):
   """Given (x,y) coordinates, calculates the angle in radians
      with respect to the positive x-axis, with the angle increasing
      in the CCW direction."""
   pi = math.pi
   pi = 3.14159265
   xabs = abs(x)
   yabs = abs(y)
   if x >= 0:
      if y >= 0:
         Angle = math.atan(yabs/xabs)
      else:
         Angle = math.atan(xabs/yabs) + 3.*pi/2.
   elif x < 0:
      if y > 0:
         Angle = math.atan(xabs/yabs) + pi/2.
      else:
         Angle = math.atan(yabs/xabs) + pi
   else:
      Angle = numpy.nan
   return (Angle)

def ellipse(ra, rb, ang, x0, y0, Nb=50):
   """ellipseGenerator will produce an ideal ellipse given the 
      semi-major and semi-minor axes, the angle of inclination,
      the center (x,y) coordinates and the number of points on
      the ellipse.  from D.G. Long, BYU from Peter Blattner, 
      U.Neuchatel, Switz.
      ra, rb = major axis length, minor axis length,
      ang = rotation angle, (x0,y0) = center coordinates,
      Nb = # of points on the ellipse"""

   x = numpy.zeros(Nb)
   y = numpy.zeros(Nb)
   co, si = math.cos(ang), math.sin(ang)
   the = numpy.linspace(0,2*math.pi, Nb)
   for i in range(Nb):
      x[i] = ra * math.cos(the[i]) * co - si * rb * math.sin(the[i]) + x0
      y[i] = ra * math.cos(the[i]) * si + co * rb * math.sin(the[i]) + y0
   return x,y


def fitellipse( x, **kwargs ):  
   """ This alorithm for the ellipsoid method of Khachiyan for 2D point
       clouds.  This version in python will return: 1) the ellipse radius
       in x, 2) the ellipse radius in y, 3) the centre in x, 4) the centre
       in y, and 5) orientation of the ellipse, and the coordinates of the
       ellipse. x is a 2xN array [x,y] points.  Bob Abel, 7/9/12"""
   x = numpy.asarray( x )  
   ## Parse inputs  
   ## Default parameters  
   kwargs.setdefault( 'constraint', 'bookstein' )  
   kwargs.setdefault( 'maxits', 200 )  
   kwargs.setdefault( 'tol', 1e-5 )  
   #if numpy.shape(x[1]) == 2:  
   if x.shape[1] == 2:  
     x = x.T  
   centroid = numpy.mean(x, 1)  
   x    = x - centroid.reshape((2,1))  
   ## Obtain a linear estimate  
   if kwargs['constraint'] == 'bookstein':  
     ## Bookstein constraint : lambda_1^2 + lambda_2^2 = 1  
     z, a, b, alpha = fitbookstein(x)  
   ## Add the centroid back on  
   if z != None:
      z = z + centroid  
   return z, a, b, alpha  

def fitbookstein(x):  
   '''  
   function [z, a, b, alpha] = fitbookstein(x)  
   FITBOOKSTEIN  Linear ellipse fit using bookstein constraint  
     lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A  
   '''  
   ## Convenience variables  
   m = x.shape[1]  
   x1 = x[0, :].reshape((1,m)).T  
   x2 = x[1, :].reshape((1,m)).T  
   ## Define the coefficient matrix B, such that we solve the system  
   ## B *[v; w] = 0, with the constraint norm(w) == 1  
   B = numpy.hstack([ x1, x2, numpy.ones((m, 1)), numpy.power( x1, 2 ), 
                    numpy.multiply(math.sqrt(2) * x1, x2 ), numpy.power( x2, 2 ) ])  
   ## To enforce the constraint, we need to take the QR decomposition  
   Q, R = numpy.linalg.qr(B)  
   ## Decompose R into blocks  
   R11 = R[0:3, 0:3]  
   R12 = R[0:3, 3:6]  
   R22 = R[3:6, 3:6]  
   ## Solve R22 * w = 0 subject to norm(w) == 1  
   U, S, V = numpy.linalg.svd(R22)  
   V = V.T  
   w = V[:, 2]  
   ## Solve for the remaining variables  
   v = numpy.dot(numpy.linalg.solve( -R11, R12 ), w )  
   ## Fill in the quadratic form  
   A    = numpy.zeros((2,2))  
   A.ravel()[0]   = w.ravel()[0]  
   A.ravel()[1:3] = 1 / math.sqrt(2) * w.ravel()[1]  
   A.ravel()[3]   = w.ravel()[2]  
   bv    = v[0:2]  
   c    = v[2]  
   ## Find the parameters  
   z, a, b, alpha = conic2parametric(A, bv, c)
   print 'z, a, b, alpha = ', z, a, b, alpha
   return z, a, b, alpha  


def conic2parametric(A, bv, c):  
   '''  
   function [z, a, b, alpha] = conic2parametric(A, bv, c)  
   '''  
   ## Diagonalise A - find Q, D such at A = Q' * D * Q  
   D, Q = numpy.linalg.eig(A)  
   Q = Q.T  
   ## If the determinant < 0, it's not an ellipse  
   if numpy.prod(D) <= 0:  
     #raise RuntimeError, 'fitellipse:NotEllipse Linear fit did not produce an ellipse'
     z, a, b, alpha = None, None, None, None
     return z, a, b, alpha
   ## We have b_h' = 2 * t' * A + b'  
   t = -0.5 * numpy.linalg.solve(A, bv)  
   c_h = numpy.dot(numpy.dot( t.T, A ), t ) + numpy.dot( bv.T, t ) + c  
   z = t  
   a = math.sqrt(-c_h / D[0])  
   b = math.sqrt(-c_h / D[1])  
   alpha = math.atan2(Q[0,1], Q[0,0])  
   return z, a, b, alpha  

def uncMagCalc(redshift, tempGalwav, tempGalfla, multiFlux=None, finish=None):
   tempName = Sed(wavelen=tempGalwav, flambda=tempGalfla)
   if multiFlux is not None:
      tempName.multiplyFluxNorm(multiFlux)
   if finish is not None:
      tempName.redshiftSED(redshift, dimming=True)
      tempName.resampleSED(wavelen_min=300, wavelen_max=1200, wavelen_step=wavelen_step)
      tempName.flambdaTofnu()
   return tempName


def fileNames(fileLoop):
    """This function searches for LSST FITS files with associated trimfiles.  obss is a list 
       of the FITS subdirectories in the /astro/net/pogo3/.../obs/ directory.  They start with a v,
       then 9 numbers, then -fu, -fg, -fr, -fi, -fz, or -fy.  Each subdirectory has an E000 path and an 
       E001 path; they simulate consecutive exposures of the same area of sky.  Each of those has
       a subdirectory, R22 (the central raft in the camera) containing one zipped eimage (compressed
       FITS file).  There are 5287 v...-f* directories.  Load each directory in the obs directory 
       that starts with a 'v', name it and add '/E000/R22/' to it and call it the eDir.  For this run
       I only want images from the r filters.  Not all of the files have trim files associated with
       them, so I also skip the ones that don't."""

    #ePath = '/astro/net/pogo3/rgibson/testSim/testPM/obs/v858247512-fr/E000/R22/eimage_858247512_R22_S11_E000.fits.gz'
    ePath = '/astro/net/pogo3/rgibson/testSim/testPM/obs/v858903212-fr/E000/R22/eimage_858903212_R22_S11_E000.fits.gz'
    #eFile = 'eimage_858247512_R22_S11_E000.fits.gz'
    eFile = 'eimage_858903212_R22_S11_E000.fits.gz'
    trimPath = None; metaPath = None; eEPath = None
    trimDir = '/astro/net/pogo3/rgibson/testSim/testPM/trim/obsid'
    t0 = eFile.split('_')
    t1 = t0[1][0:-1]
    metaPath = '%s%s/metadata_%s.dat' % (trimDir, t1, t1)
    trimPath = '%s%s/pops/trim_%s_MSSTARS.dat.gz' % (trimDir, t1, t1)
    eEPath = '%s%s/pops/trim_%s_EASTEREGGS.dat.gz' % (trimDir, t1, t1)

    outDir = '/astro/net/pogo3/rgibson/testSim/testPM/out%s/' % eFile
    # [0:-3] removes the ".gz"
    newEFile = 'eImageTANSIPa' + eFile[0:-3]
    newEPath = outDir + newEFile
        
    fileLoop = 1
    # Get the FITS image and find the positive footprints therein.
    eIm = afwImage.ImageF(newEPath)
    metaData = afwImage.readMetadata(newEPath)
    wcs = afwImage.makeWcs(metaData)
    eHDUList = pyfits.open(newEPath)
    x0 = wcs.pixelToSky(0, 0).getPosition()[0]
    y0 = wcs.pixelToSky(0, 0).getPosition()[1]
    x1 = wcs.pixelToSky(0, 4096).getPosition()[0]
    y1 = wcs.pixelToSky(0, 4096).getPosition()[1]
    x2 = wcs.pixelToSky(4096, 4096).getPosition()[0]
    y2 = wcs.pixelToSky(4096, 4096).getPosition()[1]
    x3 = wcs.pixelToSky(4096, 0).getPosition()[0]
    y3 = wcs.pixelToSky(4096, 0).getPosition()[1]
    # You just want the minimum radius so you don't go off the edge of the image boundary.
    deltax = min(abs(x2-x1), abs(x2-x0), abs(x3-x1), abs(x3-x0))
    deltay = min(abs(y1-y0), abs(y2-y0), abs(y1-y3), abs(y2-y3))
    radius = min(deltax/2., deltay/2.)
    ra_center = wcs.pixelToSky(2048, 2048).getPosition()[0]
    dec_center = wcs.pixelToSky(2048, 2048).getPosition()[1]
    return ra_center, dec_center, radius, eIm, eHDUList, wcs, ePath, fileLoop




def findGalaxies(ra_center, dec_center, radius, rMagMax, baMin, baMax):

   # CONSTANTS
   rad2deg = 180./math.pi
   csize = 100000
   obj = "ASSEMBLEDGALAXY"
   expmjd = 0.
   gal_count = 0
   wavelen_step = 0.25
   Sixflts = ['float','float','float','float','float','float']


   # BANDPASSES 
   # Instantiate bandpasses and read into lsstbp.
   bpdir = os.getenv("LSST_THROUGHPUTS_BASELINE")
   filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
   lsstbp = {}
   for f in filterlist:
      lsstbp[f] = Bandpass(wavelen_min=300,wavelen_max=1200, wavelen_step=0.1)
      lsstbp[f].readThroughput(os.path.join(bpdir, 'total_' + f + '.dat'), wavelen_step=wavelen_step)

   # SEDS
   # Read in all of the galaxy seds in the root directory:
   galdir = "/astro/net/lsst1/shared/data/galaxySED/"
   gals = {}
   gallist = os.listdir(galdir)
   for i in range(len(gallist)):
      if gallist[i].endswith('.gz'):
         gallist[i] = gallist[i][:-3]
   for gal in gallist:
      gals[gal] = Sed()
      gals[gal].readSED_flambda(os.path.join(galdir, gal))

   # Check on resampling - want all galaxy seds to have the same wavelength range.
   # (although our LSST ranges are 300 - 1200 nm, the SED's are ~27 to 2290 nm).
   if ((gals[gallist[0]].wavelen.min() < 30) & (gals[gallist[0]].wavelen.max() > 2000)):
       # If true, then gals[gallist[0]] is okay to use as a template -- this ought to be true.
       wavelen_match = gals[gallist[0]].wavelen
   else:
       print "Had to use simple wavelength array for matching"
       wavelen_match = numpy.arange(30, 2200, 0.1, dtype='float')
   for gal in gallist:
       if gals[gal].needResample(wavelen_match = wavelen_match):
           gals[gal].resampleSED(wavelen_match = wavelen_match)

   # Create the galactic a, b values from CCM 89 for the source galaxy and ours.
   # adjust for redshift, add Milky Way dust (recalculate a and b),
   # normalize fluxes and calculate magnitudes:
   # First: calculate source galaxy a/b on wavelength range required for 
   # source galaxy (internal)from CCM 89.
   a_int, b_int = gals[gallist[0]].setupCCMab()

   # Second: calculate milky way a/b on wavelength range required for calculating 
   # magnitudes - i.e. 300 to 1200 nm.
   # Set up a Sed object that is the same for all galaxies, all chunks.
   # Start with a flat SED with F_AB = 3631 Jy 
   tmpgal = Sed()
   tmpgal.setFlatSED(wavelen_min=300, wavelen_max=1200, wavelen_step=wavelen_step)
   a_mw, b_mw = tmpgal.setupCCMab()  # so this is a/b on 300-1200 range. 

   # Set up phi, the wavelength-normalized system response for each filter, for each bandpass.
   # sb is the system response function (throughputs).  Also set up a bandpass list, for
   # manyMagCalc method and initiate mags w/dust.
   bplist = []
   for f in filterlist:
      lsstbp[f].sbTophi()
      bplist.append(lsstbp[f])
   phiarray, dlambda = tmpgal.setupPhiArray(bplist)

   objId = numpy.empty(0)
   ra = numpy.empty(0)
   dec = numpy.empty(0)
   diskFluxNorm = numpy.empty(0)
   bulgeFluxNorm = numpy.empty(0)
   diskSedFilename = numpy.empty(0)
   bulgeSedFilename = numpy.empty(0)
   a_disk = numpy.empty(0)
   b_disk = numpy.empty(0)
   a_bulge = numpy.empty(0)
   b_bulge = numpy.empty(0)
   pa_d = numpy.empty(0)
   rdshft = numpy.empty(0)
   Rmags = numpy.empty(0)

   #QUERY
   #myqdb = queryDB.queryDB(chunksize=csize,objtype=obj, filetypes=('REFERENCECATALOG',))
   myqdb = queryDB.queryDB(chunksize=csize,objtype=obj, filetypes=('TEST',))
   # Specify a circular field of view:
   ic = myqdb.getInstanceCatalogByCirc(ra_center, dec_center, radius, expmjd=expmjd)
   # Begin iteratively acquiring data
   gal_count = 0
   while ic is not None:
      gal_loop = len(ic.dataArray['raJ2000'])   
      gal_count += gal_loop
      objId = numpy.concatenate((objId, ic.dataArray['objId']), axis=0)
      ra = numpy.concatenate((ra, ic.dataArray['raJ2000']*rad2deg), axis=0)
      dec = numpy.concatenate((dec, ic.dataArray['decJ2000']*rad2deg), axis=0)
      diskFluxNorm = numpy.concatenate((diskFluxNorm, ic.dataArray['diskFluxNorm']), axis=0)
      bulgeFluxNorm = numpy.concatenate((bulgeFluxNorm, ic.dataArray['bulgeFluxNorm']), axis=0)
      diskSedFilename = numpy.concatenate((diskSedFilename, ic.dataArray['diskSedFilename']), axis=0)
      bulgeSedFilename = numpy.concatenate((bulgeSedFilename, ic.dataArray['bulgeSedFilename']), axis=0)
      a_disk = numpy.concatenate((a_disk, ic.dataArray['semiMajorDisk']), axis=0)
      b_disk = numpy.concatenate((b_disk, ic.dataArray['semiMinorDisk']), axis=0)
      a_bulge = numpy.concatenate((a_bulge, ic.dataArray['semiMajorBulge']), axis=0)
      b_bulge = numpy.concatenate((b_bulge, ic.dataArray['semiMinorBulge']), axis=0)
      pa_d = numpy.concatenate((pa_d, ic.dataArray['positionAngleDisk']), axis=0)
      rdshft = numpy.concatenate((rdshft, ic.dataArray['redshift']), axis=0)
      # Get next chunk, if it exists.
      ic = myqdb.getNextChunk()
   myqdb.closeSession()


   # Calculate galactic coordinates:
   gLon = []
   gLat = []
   for i in range(gal_count):
      gcoord = afwCoord.IcrsCoord(ra[i], dec[i]).toGalactic()
      gLon.append(gcoord.getL(afwCoord.DEGREES)*math.pi/180.)
      gLat.append(gcoord.getB(afwCoord.DEGREES)*math.pi/180.)
   ebv_mw = (EBV.calculateEbv(gLon, gLat, ebvMapNorth, ebvMapSouth, interp = True))
   del gLon
   del gLat
   print 'gLon, gLat calculated'

   # Now calculate magnitudes for each galaxy.  If you have a bulge, initiate its Sed
   # instance, multiply fnu by the bulge flux normalization and apply the bulge dust
   # model (currently zero).  Next initiate a Sed instance for the disk (if you have
   # one), multiply fnu by the disk flux normalization and apply the disk dust model.
   # If you have bulge and disk Sed's, add them together; if not you'll just use 
   # whichever one you have.  Correct for redshift, resample the Sed since now it's
   # shifted in wavelength, add the Milky Way dust and calculate magnitudes using the
   # manyMag method.  Uncorrected Magnitudes (no dust, reddening, redshifting added)
   # are calculated as well.

   uncmags = numpy.zeros(gal_count, dtype={'names':['u','g','r','i','z','y'], 'formats':Sixflts})
   raEdge = numpy.zeros(len(ra))
   decEdge = numpy.zeros(len(ra))
   option = numpy.zeros(len(ra)) + 5
   raEdge = ra - a_disk/3600.*math.sin((pa_d)*math.pi/180.)
   decEdge = dec + a_disk/3600.*math.cos((pa_d)*math.pi/180.)

   # Uncorrected Magnitudes (no dust, reddening, redshifting added)
   for i in range(gal_count):
      galdisk = diskSedFilename[i]
      galbulge = bulgeSedFilename[i]
       #raEdge = ra[i] - a_disk[i]/3600.*math.sin((pa_d[i])*math.pi/180.)
      #decEdge = dec[i] + a_disk[i]/3600.*math.cos((pa_d[i])*math.pi/180.)
      if (galbulge is not None) and (galdisk is not None):
         ba_disk = b_disk[i]/a_disk[i]
         ba_bulge = b_bulge[i]/a_bulge[i]
         baRatio = (diskFluxNorm[i]*ba_disk + bulgeFluxNorm[i]
                    *ba_bulge)/(diskFluxNorm[i] + bulgeFluxNorm[i])
         if baMin <= baRatio <= baMax:
            option[i] = 2
         else: continue
         tmpbulge = uncMagCalc(rdshft[i], gals[galbulge].wavelen, gals[galbulge].flambda, 
                               multiFlux=bulgeFluxNorm[i])
         tmpdisk = uncMagCalc(rdshft[i], gals[galdisk].wavelen, gals[galdisk].flambda,
                              multiFlux=diskFluxNorm[i])
         newgal = uncMagCalc(rdshft[i], tmpdisk.wavelen, (tmpdisk.flambda+tmpbulge.flambda), finish=1)
         tmpmags = newgal.manyMagCalc(phiarray, dlambda)
      elif galbulge is not None:
         baRatio = b_disk[i]/a_disk[i]
         if baMin <= baRatio <= baMax:
            option[i] = 0
         else: continue
         tmpbulge = uncMagCalc(rdshft[i], gals[galbulge].wavelen, gals[galbulge].flambda, 
                               multiFlux=bulgeFluxNorm[i], finish=1)
         tmpmags = tmpbulge.manyMagCalc(phiarray, dlambda)
      elif galdisk is not None:
         baRatio = b_disk[i]/a_disk[i]
         if baMin <= baRatio <= baMax:
            option[i] = 1
         else: continue
         tmpdisk = uncMagCalc(rdshft[i], gals[galdisk].wavelen, gals[galdisk].flambda,
                              multiFlux=diskFluxNorm[i], finish=1)
         tmpmags = tmpdisk.manyMagCalc(phiarray, dlambda)
      j = 0
      for f in filterlist:
         uncmags[f][i] = tmpmags[j]
         j += 1

   del diskSedFilename, bulgeSedFilename, tmpmags
   del Av_d, Rv_d, Av_b, Rv_b, rdshft, uncmags, raEdge, decEdge

   idx = numpy.where((tmpmags[2] <= rMagMax) & (option != 5))
   rmags = tmpmags[2][idx]
   objId = ojbId[idx]
   option.append(option)
   ra = ra[idx]
   dec = dec[idx]
   baRatio = baRatio[idx]
   a_disk = a_disk[idx]
   b_disk = b_disk[idx]
   a_bulge = a_bulge[idx]
   b_bulge = b_bulge[idx]
   diskFluxNorm = diskFluxNorm[idx]
   bulgeFluxNorm = bulgeFluxNorm[idx]
   pa_d = pa_d[idx]
   raEdge = raEdge[idx]
   decEdge = decEdge[idx]
   #print '%11.7f %11.7f %11.7f %11.7f' % (ra[i], dec[i], tmpmags[2], baRatio)


   print '# iterations: ', iteration + 1
   print 'Total # of galaxies: %d' % (gal_count)
   print 'Number of qualifying Galaxies within a radius of %s deg. is %s' % (radius, len(rmags))
   #print 'len(ra), len(rmags), len(baRatio), len(option): ', len(ra), len(rmags), len(baRatio), len(option)
   return ra, dec, rmags, baRatio, option, objId, a_disk, b_disk, a_bulge, b_bulge, diskFluxNorm, bulgeFluxNorm, pa_d, raEdge, decEdge

def findMatchingFootprints(eIm, eHDUList, wcs, ra, dec, rmags, baRatio, objId, option, a_disk, b_disk, a_bulge, b_bulge, diskFluxNorm, bulgeFluxNorm, pa_d, raEdge, decEdge, dxMin, dyMin, matchArcsec, plotno, ePath):
   # Get the set of footprints from the FITS image, setting the 
   # threshold by either count value or STDEV.  footprintsP is the list
   # of positive footprints for that particular FITS file, con-
   # sisting of a set of pixels.  getFootprints() returns the
   # footprints of detected objects.
   # Don't set the detection limit too bright, or sources will
   #  sometimes be missed. 
   detSetP = afwDetection.makeFootprintSet(
       eIm, afwDetection.createThreshold(2., 'value'))
       #eIm, afwDetection.createThreshold(28., 'stdev'))
   footprintsP = detSetP.getFootprints()
   nf = len(footprintsP)
   b1 = 1.999*1 - 0.327
   b4 = 1.999*4 - 0.327
   detxs = []; detys = []; nCts = []; estSizePix = []
   detMinx = []; detMaxx = []; detMiny = []; detMaxy = []
   detMaxRad2 = []
   n = 0; boxSize = 2
   minCtsBrightThresh = 1000
   #print 'Considering %i positive footprints.' % len(footprintsP)
   for fp in footprintsP:
      if n % 10000 == 0: print '%i of %i positive footprints analyzed.' % (n, nf)
      # Return "a list of BBoxs whose union contains exactly the pixels in
      # foot, neither more or less."
      # This looks at each boundary box in the list of boundary boxes (bboxes)
      # and gets the minimum and maximum x and y values of the list of coordinates,
      # and calculates the x and y lengths of the boundary box (dx, dy).
      # Lots of times they're only 1 pixel in length.  fAllFPP is a file
      # containing the center x,y coordinates of each footprint, as well as dx and dy.
      # minx, maxx, miny, and maxy will be the very minimum and maximum values of
      # all the footprints for all the boundary boxes in the list.  Create a
      # rectangular footprint:
      bboxes = afwDetection.footprintToBBoxList(fp)
      minx = 1e100; maxx = -1e100; miny = 1e100; maxy = -1e100
      for bbox in bboxes:
         x0, y0 = bbox.getMinX(), bbox.getMinY()
         x1, y1 = bbox.getMaxX(), bbox.getMaxY()
         dx = x1 - x0; dy = y1 - y0
         minx = min([minx, x0, x1]); maxx = max([maxx, x0, x1])
         miny = min([miny, y0, y1]); maxy = max([maxy, y0, y1])

      oldminx = minx; oldmaxx = maxx
      oldminy = miny; oldmaxy = maxy
      # The footprint regions are sometimes skewed, so center them
      #  as best we can by finding the maximum pixel and then weighting
      #  around that.
      # eIm.get(x,y) returns the number of counts (t0) in that pixel.
      # DS9 uses 1-based coords, so add 1
      highPixx = -1; highPixy = -1; highPixVal = -1.e100
      for x in range(minx, maxx+1):
         for y in range(miny, maxy+1):
            t0 = eIm.get(x, y)
            if t0 > highPixVal:
               highPixx = x; highPixy = y; highPixVal = t0
      # Now we know the max pixel.  Work in an aperture around it.  Use the
      #  footprints to set the maximum extent of the aperture.
      # Check to make sure the distance of the pixel in question is equal
      # to or less that the maximum radius of the aperture.
      # ctsWtSumx(or y) is the weighted center coordinate of the object 
      # in the  x (or y) direction.
      # detxs is a list of ctsWtSumx for each footprint.
      # Only keep the eImages that are above a minimum brightness threshold
      # and above a minimum x or y width in terms of pixels.
      tCts = 0; ctsWtSumx = 0; ctsWtSumy = 0; nPix = 0
      t0 = (maxx - highPixx)**2
      t1 = (minx - highPixx)**2
      maxx2 = max(t0, t1)
      t0 = (maxy - highPixy)**2
      t1 = (miny - highPixy)**2
      maxy2 = max(t0, t1)
      maxRad2 = maxx2 + maxy2
      dx = abs(maxx - minx)
      dy = abs(maxy - miny)
      if dx >= dxMin or dy >= dyMin:
         for x in range(minx, maxx+1):
            for y in range(miny, maxy+1):
               r2 = (x-highPixx)**2 + (y-highPixy)**2
               if r2 > maxRad2: continue
               t0 = eIm.get(x, y)
               tCts += t0
               ctsWtSumx += x * t0; ctsWtSumy += y * t0
               nPix += 1
         if tCts > minCtsBrightThresh:
            ctsWtSumx /= float(tCts); ctsWtSumy /= float(tCts)
            # Now iterate it again over a smaller region to improve
            # your chances of hitting the center.
            minx = int(ctsWtSumx - dx/5.); maxx = int(ctsWtSumx + dx/5.)
            miny = int(ctsWtSumy - dy/5.); maxy = int(ctsWtSumy + dy/5.)
            if maxx > 3999: maxx = 3999
            if minx < 1: minx = 0
            if maxy > 4071: maxy = 4071
            if miny < 1: miny = 0

            highPixx = -1; highPixy = -1; highPixVal = -1.e100
            for x in range(minx, maxx+1):
               for y in range(miny, maxy+1):
                  t0 = eIm.get(x, y)
                  if t0 > highPixVal:
                     highPixx = x; highPixy = y; highPixVal = t0

            tCts = 0; ctsWtSumx = 0; ctsWtSumy = 0; nPix = 0
            t0 = (maxx - highPixx)**2
            t1 = (minx - highPixx)**2
            maxx2 = max(t0, t1)
            t0 = (maxy - highPixy)**2
            t1 = (miny - highPixy)**2
            maxy2 = max(t0, t1)
            maxRad2 = maxx2 + maxy2
            for x in range(minx, maxx+1):
               for y in range(miny, maxy+1):
                  r2 = (x-highPixx)**2 + (y-highPixy)**2
                  if r2 > maxRad2: continue
                  t0 = eIm.get(x, y)
                  tCts += t0
                  ctsWtSumx += x * t0; ctsWtSumy += y * t0
                  nPix += 1
            ctsWtSumx /= float(tCts); ctsWtSumy /= float(tCts)
            minx = oldminx
            maxx = oldmaxx
            miny = oldminy
            maxy = oldmaxy
            t0 = (maxx - ctsWtSumx)**2
            t1 = (minx - ctsWtSumx)**2
            maxx2 = max(t0, t1)
            t0 = (maxy - ctsWtSumy)**2
            t1 = (miny - ctsWtSumy)**2
            maxy2 = max(t0, t1)
            maxRad2 = maxx2 + maxy2

            detxs.append(ctsWtSumx); detys.append(ctsWtSumy)
            nCts.append(tCts); estSizePix.append(nPix)
            detMinx.append(int(ctsWtSumx-dx/2.)); detMaxx.append(int(ctsWtSumx+dx/2.))
            detMiny.append(int(ctsWtSumy-dy/2.)); detMaxy.append(int(ctsWtSumy+dy/2.)) 
            detMaxRad2.append(maxRad2)
      n += 1

   # At this point you have data for bright enough, wide enough footprints.
   detxs = numpy.array(detxs); detys = numpy.array(detys)
   nCts = numpy.array(nCts); estSizePix = numpy.array(estSizePix)
   detMinx = numpy.array(detMinx); detMaxx = numpy.array(detMaxx)
   detMiny = numpy.array(detMiny); detMaxy = numpy.array(detMaxy)

   #brightRAs = numpy.zeros(len(detxs))
   #brightDecs = numpy.zeros(len(detxs))
   brightRAs = wcs.pixelToSky(detxs, detys).getPosition()[0]
   brightDecs = wcs.pixelToSky(detxs, detys).getPosition()[1]
        

   # for each value of your original coordinates of interest,
   # create a list (t0) of the square of the distance between each bright
   # coordinate and the galaxy coordinates.  t1 is the index of the minimum distance
   # squared, so t0[t1] is the minimum distance squared. The match to your
   # coordinates is brightRAs[t1] and BrightDecs[t1].
   #matchDistArcsec = numpy.zeros(len(ra))id
   psi = []
   psiFit = []
   psiCatRaDec = []
   psiCatXY = []
   psiImageRaDec = []
   psiImageXY = []
   pa_d = []
   aFit = []
   bFit = []
   baRatioFit = []
   idx = []
   tidx = []
   distArcsec = []
   gal_code = ['Lowba_bulge_Only', 'Lowba_disk_Only', 'Lowba_bulge_and_Disk']
   #gal_code = ['Bulge_Only', 'Disk_Only', 'Bulge_and_Disk']
   for i in range(len(ra)):
      # Assuming distance is small, keep it linear
      t0 = (brightRAs - ra[i])**2 + (brightDecs - dec[i])**2
      t1 = numpy.argmin(t0)
      diffArcsec = numpy.sqrt(t0[t1]) * 3600.
      if diffArcsec <= matchArcsec:
         idx.append(i)
         tidx.append(t1)
         distArcsec.append(diffArcsec)
           
   ra = ra[idx]
   dec = dec[idx]
   raEdge = raEdge[idx]
   decEdge = decEdge[idx]
   raCat = raEdge - ra
   decCat = decEdge - dec
   psiCatRaDec.append(positionAngle(raCat, decCat))

   xyCtrCat = wcs.skyToPixel(ra, dec)
   xCtrCat, yCtrCat = xyCtrCat
   xyPixPos = wcs.skyToPixel(raEdge, decEdge)
   xPixPos, yPixPos = xyPixPos

   xCat = xPixPos - xCtrCat
   yCat = yPixPos - yCtrCat
   psiCatXY.append(positionAngle(xCat, yCat))

   distArcsec = numpy.array(distArcsec)
   brightRAs = brightRAs[tidx]
   brightDecs = brightDecs[tidx]
   ra = ra[idx]
   dec = dec[idx]
   raEdge = raEdge[idx]
   decEdge = decEdge[idx]
   detMinx = detMinx[tidx]
   detMaxx = detMaxx[tidx]
   detMiny = detMiny[tidx]
   detMaxy = detMaxy[tidx]
   xwidth = abs(detMaxx - detMinx)
   ywidth = abs(detMaxy - detMiny)
   detxs = detxs[tidx]
   detys = detys[tidx]
   pa_d = pa_d[idx]
   Rad2Max = max((xwidth/2.)**2, (ywidth/2.)**2)
   nbins = int(max(xwidth/2., ywidth/2.))
   rmags = rmags[idx]
   baRatio = baRatio[idx]
   option = option[idx]
   objId = objId[idx]
   a_disk = a_disk[idx]
   b_disk = b_disk[idx]
   a_bulge = a_bulge[idx]          #for m in range(nbins[i]):
   b_bulge = b_bulge[idx]
   diskFluxNorm = diskFluxNorm[idx]
   bulgeFluxNorm = bulgeFluxNorm[idx]

   # Now you've got the object and its radius.  To calculate the flux as a
   # function of radius, I(r), break the radius squared into multiple segments
   # and bin the number of counts.  Then divide the total number of counts
   # in each bin by the area, if you want. 
   print 'len(detxs), len(objId), len(rmags), len(a_disk), len(detMinx), len(detMaxy): '
   print len(detxs), len(objId), len(rmags), len(a_disk), len(detMinx), len(detMaxy)
   for i in range(len(detxs)): 
      if i != -1:
         print 'Galaxy ID:  ', objId[i]
         print 'rmags:   ', rmags[i]
         print 'a_disk:  ', a_disk[i]
         print 'b_disk:  ', b_disk[i]
         print 'a_bulge: ', a_bulge[i]
         print 'b_bulge: ', b_bulge[i]
         print 'diskFluxNorm:  ', diskFluxNorm[i]
         print 'bulgeFluxNorm:  ', bulgeFluxNorm[i]
         print 'detMinx:  ', detMinx[i]
         print 'detMaxx:  ', detMaxx[i]
         print 'detMiny:  ', detMiny[i]
         print 'detMaxy:  ', detMaxy[i]
         print 'pa_d: ', pa_d[i]
         print 'detxs:  ', detxs[i]
         print 'detys:  ', detys[i]
         print 'ra (cat): ', ra[i]
         print 'dec (cat): ', dec[i]
         print 'raEdge, decEdge: ', raEdge[i], decEdge[i]
         print 'psiCatXY: xPixPos, xCtrCat, yPixPos, yCtrCat'
         print xPixPos[i], xCtrCat[i], yPixPos[i], yCtrCat[i]
         tCts = numpy.zeros(nbins[i])
         I_R = numpy.zeros(nbins[i])
         nPix = numpy.zeros(nbins[i])
         mew = numpy.zeros(nbins[i])
         rErr = numpy.zeros(nbins[i])
         Err_Up = numpy.zeros(nbins[i])
         Err_Dn = numpy.zeros(nbins[i])
         a2Bin = numpy.zeros(nbins[i])
         aBin = numpy.zeros(nbins[i])
         AveABin = numpy.zeros(nbins[i])
         xBdry = []
         yBdry = []
         ellipticity = b_disk[i]/a_disk[i]

         # Create an image plot for each footprint.
         C = numpy.zeros(((detMaxy[i]-detMiny[i]+1), (detMaxx[i]-detMinx[i]+1)), numpy.int)
         for y in range (detMiny[i], detMaxy[i], 1):
            for x in range(detMinx[i], detMaxx[i], 1):
               C[y-detMiny[i], x-detMinx[i]] = eIm.get(x,y)

         # Find the boundaries of an eimage using the pixel count threshold
         # as the discriminator.  Break the eimage into anglar bins (say, 
         # 50,so the distant outliers are less likely to be included)and
         # select the minimum value in each bin.  The original purpose was
         # to get data points for fitting ellipses to galactic disks.

         n = 50
         deltaAng = 2.*math.pi/float(n)
         xBdry = []
         yBdry = []          
         theta = []          
         radDist2 = []
         xVal = []
         yVal = []
         angleDeg = []
         for j in range(n):
            radDist2.append([])
            xVal.append([])
            yVal.append([])
            angleDeg.append([])

         # This section is devoted to finding psi, the angle of 
         # the ellipse's rotation CCW with respect to the positive
         # x-axis
         countMin = 12
         count = 0
         for x in range(detMinx[i] + 1, detMaxx[i]-1, 1):
            xp = float(x - detMinx[i])
            for y in range(detMiny[i], detMaxy[i], 1):
               yp = float(y - detMiny[i])
               if C[yp, xp] < countMin:
                  if C[yp, xp-1] < countMin or C[yp, xp+1] < countMin:
                     y0 = float(y) - detys[i]
                     x0 = float(x) - detxs[i]
                     Angle = xyPosAngle(x0, y0)
                     m = int(Angle/deltaAng)
                     radDist2[m].append(x0**2 + y0**2)
                     xVal[m].append(xp)
                     yVal[m].append(yp)
                     angleDeg[m].append(Angle*180./math.pi)

         for j in range(n):
            if len(radDist2[j]) >= 1.:
               t1 = numpy.argmin(radDist2[j])
               xBdry.append(xVal[j][t1])
               yBdry.append(yVal[j][t1])
               theta.append(angleDeg[j][t1])
                
         z, a, b, psiAve = fitellipse([xBdry, yBdry])
         aFit.append(a)
         bFit.append(b)
         if a >= b and a != None and b !=None:
            baRatioFit.append(b/a)
         elif b > a and a != None and b !=None:
            baRatioFit.append(a/b)
         else:
            baRatioFit.append(None)
         if psiAve != None:
            psiFit.append(float(psiAve)*180./math.pi)
         else:
            psiFit.append(None)
         print 'psiFit:  ', psiFit[i]
         if a and b and psiAve and z != None:
            xEllipse, yEllipse = ellipse(aFit[-1], bFit[-1], psiAve, detxs[i]-detMinx[i], detys[i]-detMiny[i])
            daMax2 = a*a/float(nbins[i])
         else:
            print 'a, b, psiAve, or z = None for i = ', i, a, b, psiAve, z
            psiImageXY.append(None)
            psiImageRaDec.append(None)
            print 'psiImageXY, psiImageRaDec: ', psiImageXY[i], psiImageRaDec[i]
            continue
         xMajor = aFit[-1]*math.cos(psiAve) + detxs[i]-detMinx[i]
         yMajor = aFit[-1]*math.sin(psiAve) + detys[i]-detMiny[i]
         print 'xMajor, yMajor = ', xMajor, yMajor

         x0 = xMajor - detxs[i]
         y0 = yMajor - detys[i]
         psiImageXY.append(positionAngle(x0, y0))
          
         ra3 = wcs.pixelToSky(xMajor, yMajor).getPosition()[0]
         dec3 = wcs.pixelToSky(xMajor, yMajor).getPosition()[1]
         print 'ra3, brightRAs[i], dec3, brightDecs[i] = ', ra3, brightRAs[i], brightDecs[i], dec3
         xyMajor3 = wcs.skyToPixel(ra3, dec3)
         xMajor3, yMajor3 = xyMajor3
         print "xMajor3, yMajor3 = ", xMajor3, yMajor3
         xyMajor2 = wcs.skyToPixel(raEdge[i], decEdge[i])
         xMajor2, yMajor2 = xyMajor2
         print 'brightRAs[i], brightDecs[i]:'
         print brightRAs[i], brightDecs[i]
         ra0 = ra3 - brightRAs[i]
         dec0 = dec3 - brightDecs[i]
         psiImageRaDec.append(positionAngle(ra0, dec0))
             
         t0 = 0.
         for m in range(nbins[i]):
            t0 += daMax2
            a2Bin[m] = t0
            aBin[m] = math.sqrt(t0)

         AveABin[0] = 10.**(math.log10(aBin[0])/2.)
         for m in range(nbins[i]-1):
            Ave = (math.log10(aBin[m+1]) + math.log10(aBin[m]))/2.
            AveABin[m+1] = 10.**Ave

         xbin = []
         ybin = []
         for x in range(detMinx[i], detMaxx[i]+1):
            for y in range(detMiny[i], detMaxy[i]+1):
               xprime = x - detxs[i]
               yprime = y - detys[i]
               aXY2 = (xprime*math.cos(psiAve) + yprime*math.sin(psiAve
                       ))**2 + (xprime*math.sin(psiAve) 
                       - yprime*math.cos(psiAve)/(b/a))**2
               if aXY2 > a*a: continue
               nval = min(int(aXY2/daMax2), nbins[i]-1)
               tCts[nval] += eIm.get(x,y)
               nPix[nval] += 1
               if nval <= 4:
                  xbin.append(x - detMinx[i])
                  ybin.append(y - detMiny[i])
         sumCts = 0
         sumPix = 0
         for m in range(nbins[i]):
            sumCts += tCts[m]
            sumPix += nPix[m]
            #mew, the mean, is in counts/pixel
            mew[m] = tCts[m]/nPix[m]
         for m in range(nbins[i]): 
            # I(R) in units of counts per square arcsec: 
            if mew[m] == 0:
               I_R[m] = numpy.nan
            else:
               I_R[m] = mew[m]/.04
               rErr[m] = math.sqrt(tCts[m])/nPix[m]/.04
         I_Reff = max(I_R)/math.e
         for m in range(nbins[i]):
            if I_R[m] < I_Reff:
               R_eff = 10**(math.log10(AveABin[m]) - (math.log10(I_R[m]) - math.log10(I_Reff))/(math.log10(I_R[m]) - math.log10(I_R[m-1]))*(math.log10(AveABin[m]) - math.log10(AveABin[m-1]))) 
               break 
         print 'xwidth, ywidth (in pixels) = ', abs(detMaxx[i] - detMinx[i]), abs(detMaxy[i] - detMiny[i])
         print 'Total Counts = %i, Total Pixels = %i' % (sumCts, sumPix)
         print 'log_Radius(") log_Intensity   n = 1 and 4      n = 1        n = 4    Pixels   Counts  Mean   Sigma'
         #Prepare Sersic n=1 and n = 4 model profiles for comparison
         In1 = numpy.zeros(nbins[i])
         In4 = numpy.zeros(nbins[i])
         InBoth = numpy.zeros(nbins[i])
         for m in range(nbins[i]):
           In1[m] = I_Reff*math.exp(-b1*((AveABin[m]/R_eff) - 1))
           In4[m] = I_Reff*math.exp(-b4*(((AveABin[m]/R_eff)**0.25) - 1))
           InBoth[m] = (diskFluxNorm[i]*In1[m] + bulgeFluxNorm[i]*In4[m])/(diskFluxNorm[i] + bulgeFluxNorm[i])
         for m in range(nbins[i]):
            # Convert pixels to arcsec (0.2 arcsec/pixel, roughly)
            AveABin[m] = math.log10(AveABin[m]*.2)
            if I_R[m] != None:
               Err_Up[m] = math.log10(I_R[m] + rErr[m]) - math.log10(I_R[m])
               if I_R[m] > rErr[m]:
                  Err_Dn[m] = math.log10(I_R[m]) - math.log10(I_R[m] - rErr[m])
               else:
                  Err_Dn[m] = 5
               I_R[m] = math.log10(I_R[m])
            if In1[m] != 0.: 
               In1[m] = math.log10(In1[m])
            else:
               In1[m] = numpy.nan
            if In4[m] != 0.: 
               In4[m] = math.log10(In4[m]) 
            else:
               In4[m] = numpy.nan
            if InBoth[m] != 0.: 
               InBoth[m] = math.log10(InBoth[m])
            else:
               InBoth[m] = numpy.nan
            if m > 8 and I_R[m] is 'nan' and I_R[m-1] is 'nan' and I_R[m-2] is 'nan' and I_R[m-3] is 'nan':
               nbins[i] = nbins[:(nbins[i]-m+1)]
         for m in range(nbins[i]):
            print '%11.2f %11.2f %11.2f %11.2f %11.2f  %11.2f  %11.2f %11.2f %11.2f' % (AveABin[m], I_R[m], InBoth[m], In1[m], In4[m], nPix[m], tCts[m], mew[m], rErr[m])
         fig = plt.figure()
         fig.suptitle('Galaxy %s, (ra,dec)=(%5.2f, %5.2f), %d pixels \n%d counts; Gal_ID %s R_eff = %5.2f arcsecs, PA = %5.2f deg.\n%s' %(gal_code[option[i]], brightRAs[i], brightDecs[i], sumPix, sumCts, objId[i], R_eff*.2, pa_d[i], ePath), fontsize=9)
         ax = fig.add_subplot(121)
         im = plt.imshow(C, aspect ='equal', origin='lower', cmap = plt.cm.gray, interpolation = 'nearest')
         cir = Circle((detxs[i]-detMinx[i],detys[i]-detMiny[i]), radius = .25, ec='r', fc = 'r')
         cir2 = Circle((xMajor, yMajor), radius = .25, ec='g', fc = 'g')
         print 'semi-major axis direction at (', xMajor, ', ', yMajor, ')'
         ax.add_patch(cir)
         ax.add_patch(cir2)
         line5, = ax.plot(xBdry, yBdry, 'r.')
         line6, = ax.plot(xEllipse, yEllipse, 'm.')
         line7, = ax.plot(xbin, ybin, 'y.')
         xcolname = 'log [Radius (arcsec)]'
         ycolname = 'Intensity (Log_10[Counts/Area])'
         ax2 = fig.add_subplot(122)
         plt.xlabel(xcolname, fontsize=10)
         plt.ylabel(ycolname, fontsize=10)
         line1 = ax2.errorbar(AveABin, I_R, yerr=[Err_Dn, Err_Up], fmt = 'bo', ms = 2)
         line2, = ax2.plot(AveABin, In1, 'm--')
         line3, = ax2.plot(AveABin, In4, 'g-.')
         line4, = ax2.plot(AveABin, InBoth, 'r-')
         props = font_manager.FontProperties(size=12)
         leg = ax2.legend([line1[0], line2, line3, line4], ['eImage', 'n=1', 'n=4', 'n=1+4'], loc='upper right', prop=props)
         plt.text(0.1, 0.1, 'r = %5.2f\nb/a = %5.3f' % (rmags[i], baRatio[i]), transform = ax2.transAxes)
         plt.savefig('%s/Sers70cM12%sh.png' % (gal_code[option[i]], objId[i]), format = 'png')
         plotno += 1
         plt.clf()
   print 'Length of the following variables:'
   print len(detxs), len(detys), len(xPixPos), len(yPixPos), len(psiImageXY), len(psiImageRaDec), len(psiCatXY), len(psiCatRaDec), len(psiFit), len(a_disk), len(b_disk), len(aFit), len(bFit)
   print 'objId    ra2-raEdge  dec2-decEdge   detxs   detys  xPixPos   yPixPos    a_disk b_disk   a     b b/a catalog  b/a fit'
   for i in range(len(detxs)):
      if aFit[i] and bFit[i] and a_disk[i] and b_disk[i] and psiFit[i] != None:
         print '%8i  %7.2e %7.2e %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f' % (objId[i], 
                ra2[i]-raEdge[i], dec2[i]-decEdge[i], detxs[i], 
                detys[i], xPixPos[i], yPixPos[i], a_disk[i], b_disk[i], aFit[i], bFit[i],
                b_disk[i]/a_disk[i], bFit[i]/aFit[i])
   print 'objId   pa_d psiCatRD psiImRD   psiCatXY  psiImXY  psiFit  b/a cat  b/a fit'
   for i in range(len(detxs)):
      if aFit[i] and bFit[i] and a_disk[i] and b_disk[i] and psiFit[i] != None:
         print '%8i  %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f  %7.2f' % (objId[i], pa_d[i], psiCatRaDec[i], psiImageRaDec[i], psiCatXY[i], psiImageXY[i], psiFit[i], b_disk[i]/a_disk[i], bFit[i]/aFit[i])
   return plotno

      
plotno = 1
rMagMax = 24
baMin = 0.0
baMax = 0.40
dxMin = 14
dyMin = 14
fileLoop = 0
theRange = 1
#theRange = 12
matchArcsec = 1.5
#matchArcsec = 2.
for j in range(theRange):
   print 'j = ',j, 'of ', theRange + 1, ' and fileLoop = ', fileLoop
   # Get information about an r filter FITS file.
   ra_center, dec_center, radius, eIm, eHDUList, wcs, ePath, fileLoop = fileNames(fileLoop)
   # Find galaxy catalog objects in that area that are bright and round.
   ra, dec, rmags, baRatio, option, objId, a_disk, b_disk, a_bulge, b_bulge, diskFluxNorm, bulgeFluxNorm, pa_d, raEdge, decEdge = findGalaxies(ra_center, dec_center, radius, rMagMax, baMin, baMax)
   # Find and analyze corresponding footprints if they're big enough and close enough to 
   # the galaxy catalog object's coordinates.
   plotno = findMatchingFootprints(eIm, eHDUList, wcs, ra, dec, rmags, baRatio, ojbId, option, a_disk, b_disk, a_bulge, b_bulge, diskFluxNorm, bulgeFluxNorm, pa_d, raEdge, decEdge, dxMin, dyMin, matchArcsec, plotno, ePath)
   



