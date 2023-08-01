from distutils.core import setup, Extension
from Cython.Build import cythonize

package = Extension('search_funcs', ['search_funcs.pyx'])
setup(ext_modules=cythonize([package]))