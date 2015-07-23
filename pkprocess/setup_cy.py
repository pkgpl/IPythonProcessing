from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
import numpy as np

setup(name = 'velanc',
      # This includes the NumPy headers when compiling.
      cmdclass={'build_ext':build_ext},
      include_dirs = [np.get_include()],
      ext_modules = [Extension("velan_cy",["velan_cy.pyx"])]
      )

setup(name = 'corrc',
      # This includes the NumPy headers when compiling.
      cmdclass={'build_ext':build_ext},
      include_dirs = [np.get_include()],
      ext_modules = [Extension("corr_cy",["corr_cy.pyx"])]
      )

setup(name = 'traveltime',
      # This includes the NumPy headers when compiling.
      cmdclass={'build_ext':build_ext},
      include_dirs = [np.get_include()],
      ext_modules = [Extension("traveltime_cy",["traveltime_cy.pyx"])]
      )

setup(name = 'kirchhoff1',
      # This includes the NumPy headers when compiling.
      cmdclass={'build_ext':build_ext},
      include_dirs = [np.get_include()],
      ext_modules = [Extension("kirchhoff_cy",["kirchhoff_cy.pyx"])]
      )
