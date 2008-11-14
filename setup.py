from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages
from distutils.core import Extension
import os
import numpy

""" Setup script for the roche python extension"""

library_dirs = []
include_dirs = []

# need to direct to where includes and  libraries are
if os.environ.has_key('TRM_SOFTWARE'):
    library_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'lib'))
    include_dirs.append(os.path.join(os.environ['TRM_SOFTWARE'], 'include'))
else:
    print >>sys.stderr, "Environment variable TRM_SOFTWARE pointing to location of shareable libraries and includes not defined!"

include_dirs.append(numpy.get_include())

roche = Extension('trm.roche._roche',
                  define_macros   = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '1')],
                  undef_macros    = ['USE_NUMARRAY'],
                  include_dirs    = include_dirs,
                  library_dirs    = library_dirs,
                  runtime_library_dirs = library_dirs,
                  libraries       = ['roche'],
                  sources         = [os.path.join('trm', 'roche', 'roche.cc')])

setup(name='trm.roche',
      namespace_packages = ['trm'],
      version='0.1',
      package_dir = {'trm.roche' : os.path.join('trm', 'roche')},
      packages =find_packages(),
      ext_modules=[roche],

      description='Python interface to roche geometry routines.',
      author='Tom Marsh',
      author_email='t.r.marsh@warwick.ac.uk',
      url='http://www.astro.warwick.ac.uk/',
      long_description="""
trm.roche provides an interface to a set of C++ routines for computing Roche lobe
shapes, gas streams and the like.
""",

      )

