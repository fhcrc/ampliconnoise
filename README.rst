``anoisetools``
=================

Python package for prepping 454 data for use with `AmpliconNoise`_
(`Quince et al BMC Bioinformatics 2011`_, `Quince et al Nature Methods 2009`_)::

    raw.sff -> anoisetools -> Processed Data

The source for `AmpliconNoise`_ is also included.

For flowgram data, we target the original ``.sff`` files.

Installation
------------

0. Ensure that your computer meets the minimum requirements. Currently, just
   Python 2.7, plus the requirements of AmpliconNoise.
1. Install `BioPython`_ if you don't have it. Note that numpy is *not* required
   for this project. If you're not planning to use BioPython, you can answer
   "No" when the BioPython installer prompts you about numpy.
2. Download and install::

    curl -L https://github.com/fhcrc/ampliconnoise/tarball/master | tar xjf -
    cd fhcrc-ampliconnoise-*
    python2.7 setup.py install  # may require sudo

   See `Installing Python Modules`_ for more information and options.
3. Build the AmpliconNoise binaries, and ensure they're present in your
   ``path``

Running setup.py installs the ``anoisetools`` package, mostly accessible from
the ``anoise`` script.

Overview
--------

``anoise`` is called with a subcommand::

    anoise [subcommand]

Help can be accessed via ``anoise -h`` or ``anoise <subcommand> -h``.

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

If the barcode map is named ``barcodes.csv``, and the full SFF is ``G0YK51K01.sff``,
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
    sample1-pnoise_trunc.fa \
    sample1-pnoise.mapping

Both ``pyronoise`` and ``seqnoise`` create a temporary direcory for processing.
If running MPI jobs spanning multiple nodes, be sure to set ``--temp-dir`` to a
location accessible from all.

.. _AmpliconNoise: http://code.google.com/p/ampliconnoise/
.. _Quince et al BMC Bioinformatics 2011: http://dx.doi.org/10.1186/1471-2105-12-38
.. _Quince et al Nature Methods 2009: http://dx.doi.org/10.1038/nmeth.1361
.. _Biopython: http://biopython.org/wiki/Main_Page
.. _Installing Python Modules: http://docs.python.org/install/index.html

.. [1] Note: the split step creates a child process for each sample.
