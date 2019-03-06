import sys
import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
import multiprocessing

#os.environ["CC"]="gcc-8"
#os.environ["CXX"]="g++-8"
compile_args = ['-fPIC','-Wall','-O2']


exts = [Extension("gw_analysis_tools_py.*",['src/*.pyx'],
                extra_compile_args=compile_args,
                include_dirs=[numpy.get_include(),"../include","/usr/include/",'/usr/local/include/'],
                extra_objects = ["../lib/libgwanalysistools.a"],
                libraries_dir=['/usr/local/lib/','/usr/lib64','/Users/sperkins/opt/adolc_base/lib64'],
                libraries = ['adolc','fftw3'],
                language='c++'
                )]

setup(
    name='gw_analysis_tools_py',
    packages=['gw_analysis_tools_py'],
    ext_modules=cythonize(exts,force=True),
    cmdclass={'build_ext':build_ext}
    )
