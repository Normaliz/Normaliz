#!/usr/bin/env python

from distutils.core import setup, Extension
import sys

if sys.version_info < (3,5):
    macro_list = [ ( "PYTHON_VERSION_OLDER_THREE_FIVE", "1" ) ]
else:
    macro_list = [ ]

setup(
    name = 'PyNormaliz',
    version = '1.16',
    description = 'An interface to Normaliz',
    author = 'Sebastian Gutsche, Richard Sieg',
    author_email = 'sebastian.gutsche@gmail.com',
    url = 'https://github.com/Normaliz/PyNormaliz',
    py_modules = [ "PyNormaliz" ],
    ext_modules = [ Extension( "PyNormaliz_cpp",
                              [ "NormalizModule.cpp" ],
                              extra_link_args=['-lnormaliz', '-lgmp' ],
                              define_macros = macro_list ) ],
    
    package_data = {'': [ "COPYING", "GPLv2", "Readme.md" ] },
)
