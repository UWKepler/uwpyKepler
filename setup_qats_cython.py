from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy 
#export CFLAGS=-I/usr/lib/python2.7/site-packages/numpy/core/include/
#run using
#user]$python setup.py build_ext --inplace
incdir = '/astro/apps/pkg/python64/lib/python2.6/site-packages/numpy/core/include'
libdir = '/astro/apps/pkg/python64/lib/python2.6/site-packages/numpy/core/lib'
#incdir = ''
#libdir = ''
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("qats_cython", ["qats_cython.pyx"],library_dirs=[libdir],include_dirs=[incdir,numpy.get_include()],libraries=["m"])]
)
