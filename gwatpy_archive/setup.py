import sys
import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
import multiprocessing

scriptdir = os.path.dirname(os.path.realpath(__file__))

compile_args = ['-fopenmp','-fPIC','-O3']

exts = [Extension("gwatpy.*",['src/*.pyx'],
                extra_compile_args=compile_args,
                include_dirs=[numpy.get_include(),scriptdir+"/../build/install/include/gwat","src"],
                extra_objects = [scriptdir+"/../build/install/lib/libgwat.so"],
                libraries = ['adolc','fftw3','gsl','gslcblas','gwat'],
                extra_link_args=['-fopenmp'],
                language='c++'
                )]

setup(
    name='gwatpy',
    packages=['gwatpy'],
    ext_modules=cythonize(exts,force=True,language_level=3),
    cmdclass={'build_ext':build_ext}
    )
