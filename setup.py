# run "python3 setup.py build_ext --inplace" to installing the module
from setuptools import setup, Extension

module = Extension('mysymnmfsp',
                   sources=['mysymnmfsp.c'],
                   include_dirs=[r'C:\Program Files\Python312\include'],   # Adjust this path to your Python include directory
                   library_dirs=[r'C:\Program Files\Python312\libs'],      # Adjust this path to your Python library directory
                   libraries=['python312'])                 # Adjust to your Python version (e.g., python38 for Python 3.8)

setup(name='SymNMF',
      version='1.0',
      description='SymNMF Python wrapper',
      ext_modules=[module])
