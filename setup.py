#! /usr/bin/env python

import glob
from distutils.core import setup

setup(name = 'ampliconnoise',
      version = '0.1',
      author = 'Connor McCoy',
      author_email = 'cmccoy@fhcrc.org',
      packages = ['ampliconnoise', 'ampliconnoise.test'],
      requires = ['Python (>= 2.7)'],
      scripts = ['scripts/anoise-clean', 'scripts/anoise-split'],
      )
