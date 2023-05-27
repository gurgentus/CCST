from distutils.core import setup, Extension
setup(name='omt', version='1.0',  \
ext_modules=[Extension('omt', ['omt.cpp'])])