#!/usr/bin/env python

"""Setup script for the Object3DCTools module distribution."""

import os, sys, glob
try:
    import numpy
except ImportError:
    text  = "You must have numpy installed.\n"
    text += "See http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103\n"
    raise ImportError, text

from distutils.core import setup
from distutils.extension import Extension
if sys.platform == 'win32':
    libraries = ['opengl32', 'glu32']
    define_macros = [('WIN32', None)]
else:
    libraries = ['GL', 'GLU']
    define_macros = []

sources = glob.glob('*.c')
#sources = glob.glob('*.cpp')
setup (
        name         = "Object3DCTools",
        version      = "1.0",
        description  = "Object3D helper module",
        author       = "V.A. Sole - ESRF",
        author_email = "sole@esrf.fr",
        url          = "http://www.esrf.fr/computing/bliss/",

        # Description of the modules and packages in the distribution
        ext_modules  = [
                       Extension(
                            name          = 'Object3DCTools',
                            sources       = sources,
                            define_macros = define_macros,
                            libraries  = libraries, 
                            include_dirs  = [numpy.get_include()]
                       ),
       ],
)
