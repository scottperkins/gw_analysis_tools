from distutils.core import setup, Extension

gwat = Extension("testmodule",sources=["../build/install/lib/libgwat.so"])

setup(name="GWATPy",
    version='1.0',
    description="Description",
    ext_modules=[gwat])
