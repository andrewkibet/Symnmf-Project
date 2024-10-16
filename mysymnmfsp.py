from setuptools import setup, Extension
import numpy

module = Extension('mysymnmfsp',
                   sources=['mysymnmfsp.c', 'symnmf.c'],  # Your C files
                   include_dirs=[numpy.get_include(), r'C:\Program Files\Python312\include'],
                   library_dirs=[r'C:\Program Files\Python312\libs'],
                   libraries=['python312'])

setup(name='SymNMF',
      version='1.0',
      description='SymNMF Python wrapper',
      ext_modules=[module])
