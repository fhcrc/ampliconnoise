``anoisetools``
===============

Python package for prepping 454 data for use with `AmpiclonNoise`_::

    raw.sff -> flower -> anoisetools -> AmpiclonNoise


Installation
------------

1. Make sure your system meets the requirements (currently, just Python version 2.7)
2. Run::

    git clone git://github.com/fhcrc/ampiclonnoise.git
    cd ampiclonnoise && python2.7 setup.py

Scripts
-------

Running setup.py installs:

``anoise_split``
  Split sequencing results based on sample barcodes

``anoise_clean``
  Clean bad reads, trim long reads, and require a minimum length


.. _AmpiclonNoise: http://code.google.com/p/ampliconnoise/
