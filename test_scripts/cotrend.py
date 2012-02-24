#This script does the cotrending for a specified quarter of data, assuming the 
#data and cotrending cbv files are contained in 'path' 
#(default path is ben's scratch folder)
#It spits out cotrended data in a new folder in the path directory

q = '6'		#the quarter to be updated
path = "/astro/store/student-scratch1/bvegaff/"	#location of the data

#Unfortunately, there are some preliminary steps:
#1. log into darkstar and do all the setup stuff 
#setenv iraf/local/iraf
#setenv LSST_DEVEL /astro/store/shared-scratch1/acbecker/LSST/lsst_devel
#source /local/acbecker/winter2012/loadLSST.csh
#setup numpy
#setup scipy
#setup matplotlib
#setup stsci 
#2. run this script using this command ->
#		$STSCI_DIR/bin/pyraf -s < /***wherever you've pulled this to***/cotrend.py
#
#
#this script:
#Goes thru the .fits files in the unzipped quarter data,
#sends them into the kepcotrend routine in the kepler lsst package,
#using the cbv files in path.
#-remember- need an extra line of white space at end for pyraf to read the script right
#The outfile has all of the original data plus the cotrended data in an extra column.

import os, re, sys

#the needed stuff for using the kepler package
reset kepler=/astro/users/acbecker/LSST/lsst_devel/PyKE/2.1.3/
task kepler.pkg = kepler$kepler.cl
kepler

print q
qpath = path + 'Q' + q + '_public'
broke = open(path + 'Q' + q + '_broke', 'w')
if q == '0':
	cbvfile = path + "cbv/kplr2009131105131-q00-d05_lcbv.fits"

if q == '1':
	cbvfile = path + "cbv/kplr2009166043257-q01-d05_lcbv.fits"

if q == '2':
	cbvfile = path + "cbv/kplr2009259160929-q02-d07_lcbv.fits"

if q == '3':
	cbvfile = path + "cbv/kplr2009350155506-q03-d04_lcbv.fits"

if q == '4':
	cbvfile = path + "cbv/kplr2010078095331-q04-d06_lcbv.fits"

if q == '5':
	cbvfile	= path + "cbv/kplr2010174085026-q05-d08_lcbv.fits"

if q == '6':
	cbvfile	= path + "cbv/kplr2010265121752-q06-d09_lcbv.fits"

if q == '7':
	cbvfile	= path + "cbv/kplr2010355172524-q07-d10_lcbv.fits"

if q == '8':
	cbvfile = path + "cbv/kplr2011073133259-q08-d11_lcbv.fits"

if q == '9':
	cbvfile	= path + "cbv/kplr2011177032512-q09-d12_lcbv.fits"

if q == '10':
	cbvfile	= path + "cbv/kplr2011271113734-q10-d13_lcbv.fits"


files = os.listdir(qpath)

for f in files:
	outf = path + 'Q' + q + '_cotrended/' + re.sub("_llc.fits", "_llc_cbv.fits", f)
	print f
	try:
		iraf.kepcotrend(infile = os.path.join(qpath, f),
                   outfile = outf,
                   cbvfile = cbvfile,
                   vectors = "1 2 3",
                   method = "llsq",
                   fitpower = 1.0,
                   iterate = "no",
                   sigmaclip = 3.0,
                   maskfile = "",
                   scinterp = "cubic",
                   plot = "no",
                   clobber = "yes")
	except:
		broke.write('failed: ' + f +'\n')

broke.close()

.exit
