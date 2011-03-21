#! /usr/bin/env python

import glob
from distutils.core import setup

setup(name = 'ampiclonnoise',
      version = '0.1',
      author = 'Connor McCoy',
      author_email = 'cmccoy@fhcrc.org',
      packages = ['ampiclonnoise', 'ampiclonnoise.test'],
      requires = ['Python (>= 2.7)'],
      scripts = ['scripts/anoise_clean', 'scripts/anoise_split'],
      )
