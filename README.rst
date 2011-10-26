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

Help can be accessed via ``anoise -h`` or ``anoise subcommand -h``.

For our analyses, initial preprocessing is two steps:

* Split the original ``.sff`` file into one ``.sff`` per barcoded sample
* Process each sample using wrappers for PyroNoise and SeqNoise

Splitting Sequences
^^^^^^^^^^^^^^^^^^^

To split an ``.sff``, use ``anoise split``, providing a file with comma-delimited
``base_path_for_output,barcode,primer`` records, e.g.::

    sample1/sample1,ATAG,TAAATGGCAGTCTAGCAGAARAAG

will fill ``./sample1/sample1.sff``, with all sequences starting with
``ATAG``, followed by ``TAAATGGCAGTCTAGCAGAARAAG``.

Degenerate primers should be specified as such.

If the barcode map is names ``barcodes.csv``, and the full SFF is ``G0YK51K01.sff``,
one would call [1]_::

  anoise split barcodes.csv G0YK51K01.sff

Running PyroNoise and SeqNoise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each sample in our analyses, we follow a process along the lines of::

  #!/bin/sh
  MPIARGS="-np 12"
  TMP_DIR="."

  # Run PyroNoise
  # This cleans flowgrams prior
  anoise pyronoise \
    --mpi-args "$MPIARGS" \
    --temp-dir $TMP_DIR \
    sample1.sff

  anoise truncate "{barcode}" 400 < sample1-pnoise_cd.fa > sample1-pnoise_trunc.fa

  # Run SeqNoise
  anoise seqnoise \
    --mpi-args "$MPIARGS" \
    --stub sample1 \
    --temp-dir $TMP_DIR \
    sample-pnoise_trunc.fa \
    sample-pnoise.mapping

Both ``pyronoise`` and ``seqnoise`` create a temporary direcory for processing.
If running MPI jobs spanning multiple nodes, be sure to set ``--temp-dir`` to a
location accessible from all.

.. _AmpliconNoise: http://code.google.com/p/ampliconnoise/
.. _Quince et al BMC Bioinformatics 2011: http://dx.doi.org/10.1186/1471-2105-12-38
.. _Quince et al Nature Methods 2009: http://dx.doi.org/10.1038/nmeth.1361
.. _Biopython: http://biopython.org/wiki/Main_Page

.. [1] Note: the split step creates a child process for each sample.
