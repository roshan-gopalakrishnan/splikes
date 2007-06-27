from distutils.core import setup
from distutils.extension import Extension
from Pyrex.Distutils import build_ext

import numpy

setup(
  name = 'main',
  ext_modules=[ 
    Extension("simutils",["simutils.pyx"],
    include_dirs = [numpy.get_include()]),
    ],
  cmdclass = {'build_ext': build_ext}
)







