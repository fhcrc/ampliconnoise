``anoisetools``
=================

Python package for prepping 454 data for use with `AmpliconNoise`_
(`Quince et al BMC Bioinformatics 2011`_, `Quince et al Nature Methods 2009`_)::

    raw.sff -> anoisetools -> Processed Data

The source for `AmpliconNoise`_ is also included.

For flowgram data, we target the original ``.sff`` files.


Installation
------------

1. Make sure your system meets the requirements (currently, just Python
   version 2.7 and `Biopython`_)
2. Run::

    git clone git://github.com/fhcrc/ampliconnoise.git
    cd ampliconnoise && python2.7 setup.py

Running setup.py installs the ``anoisetools`` package, mostly accessible from
the ``anoise`` script.

Overview
--------

``anoise`` is called with a subcommand::

    anoise [subcommand]

In our analyses, we might follow a process along the lines of::

  #!/bin/sh
  MPIARGS="-np 12"

  # Run PyroNoise
  anoise pyronoise \
    --mpi-args "$MPIARGS" \
    --temp-dir $TMP_DIR \
    sample.sff

  anoise truncate "{barcode}" 400 < sample-pnoise_cd.fa > sample-pnoise_trunc.fa

  # Run SeqNoise
  anoise seqnoise \
    --mpi-args "$MPIARGS" \
    --stub sample \
    --temp-dir $TMP_DIR \
    sample-pnoise_trunc.fa \
    sample-pnoise.mapping


.. _AmpliconNoise: http://code.google.com/p/ampliconnoise/
.. _Quince et al BMC Bioinformatics 2011: http://dx.doi.org/10.1186/1471-2105-12-38
.. _Quince et al Nature Methods 2009: http://dx.doi.org/10.1038/nmeth.1361
.. _Biopython: http://biopython.org/wiki/Main_Page
