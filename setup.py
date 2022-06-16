from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy 
from Cython.Distutils import build_ext

import subprocess as sbp
import os.path as osp

import os

GSL_INCLUDE = "/home/sugiymnn/cosmo/gsl/include"
GSL_LIB = "/home/sugiymnn/cosmo/gsl/lib"

FFTW_INCLUDE = "/home/sugiymnn/cosmo/fftw3/include"
FFTW_LIB     = "/home/sugiymnn/cosmo/fftw3/lib"

hitomipy_ext = Extension("hitomipy", sources=["hitomi.pyx"], 
                    extra_compile_args=['-std=c++11'],
                    include_dirs=[".", numpy.get_include(), GSL_INCLUDE, FFTW_INCLUDE],
                    library_dirs=[".", GSL_LIB, FFTW_LIB],
                    libraries=['m', 'gsl'],language="c++")

setup(name="hitomipy", cmdclass = {'build_ext': build_ext}, ext_modules=[hitomipy_ext])
