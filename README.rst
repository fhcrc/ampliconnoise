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

Running setup.py installs the ``ampiclonnoise`` package, most accessible from 
the ``anoise`` script.

``anoise`` needs a subcommand to be useful:

* ``clean``: Trim bad flowgrams, enforce a minimum number of flows, remove any
  flows beyond a defined maximum
* ``split``: Split sequences from a flower-decoded SFF file by tag sequence, 
  writing the results to a file for each tag
* ``truncate``: Trim tags from FASTA sequences, truncate length


.. _AmpiclonNoise: http://code.google.com/p/ampliconnoise/
