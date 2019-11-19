import sys
import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
import multiprocessing
scriptdir = os.path.dirname(__file__)
sys.path.append(os.path.abspath("."))
sys.path.append(os.path.abspath("./src"))
#sys.path.append(scriptdir+"/src")

#os.environ["CC"]="gcc-9"
#os.environ["CXX"]="g++-9"
compile_args = ['-fopenmp','-fPIC','-Wall','-O2']

exts = [Extension("gwatpy.*",['src/*.pyx'],
                extra_compile_args=compile_args,
                include_dirs=[numpy.get_include(),"../include/gwat","./src"],
                extra_objects = ["../lib/libgwat.so"],
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
