#! /usr/bin/env python

import glob
from distutils.core import setup

setup(name = 'anoisetools',
      version = '0.1',
      author = 'Connor McCoy',
      author_email = 'cmccoy@fhcrc.org',
      packages = ['anoisetools', 'anoisetools.test', 'anoisetools.scripts'],
      requires = ['Python (>= 2.7)'],
      scripts = glob.glob('scripts/*'),
      )
