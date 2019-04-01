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

#os.environ["CC"]="gcc-7"
#os.environ["CXX"]="g++-7"
compile_args = ['-fPIC','-Wall','-O2']


exts = [Extension("gw_analysis_tools_py.*",['src/*.pyx'],
                extra_compile_args=compile_args,
                include_dirs=[numpy.get_include(),"../include","./src"],
                extra_objects = ["../lib/libgwanalysistools.a"],
                libraries = ['adolc','fftw3'],
                language='c++'
                )]

setup(
    name='gw_analysis_tools_py',
    packages=['gw_analysis_tools_py'],
    ext_modules=cythonize(exts,force=True,language_level=3),
    cmdclass={'build_ext':build_ext}
    )
