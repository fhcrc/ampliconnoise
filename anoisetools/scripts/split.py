
# Class
#  - Read barcode file
#  - Generate list of file descriptors
#  - Read flower input
#  - Find a destination based on barcode and primer
#  - Write results to appropriate output file

import collections
import csv
import re

def _load_barcodes(fp):
    """
    Read barcodes from fp

    Raises a ValueError on duplicate key *or* value.

    fp should contain (at least) sequence_identifier,barcode_sequence
    pairs.
    """
    reader = csv.reader(fp)
    values = set()
    d = {}
    for line in reader:
        k, v = line[1], line[0]
        if k in d:
            raise ValueError("Duplicate barcode: {0}".format(k))
        if v in values:
            raise ValueError("Duplicate value: {0}".format(v))
        d[k] = v
        values.add(v)
    return d

def _close_all(files):
    """
    Closes an iterable of files, aggregating any exceptions
    """
    exceptions = []
    for f in files:
        try:
            f.close()
        except Exception, e:
            exceptions.append(e)

    if exceptions:
        raise IOError("Could not close {0} files: {1}".format(
            len(exceptions), exceptions))


class SFFFlowerSplitter(object):

    def __init__(self, barcode_map, primer, dest_dir, unmatched_dest):
        """
        Initialize a new instance

        :param map barcode_map: map from sequence barcode to name
        :param str primer: Primer sequence
        :param str dest_dir: destination directory for output
        :param str unmatched_dest: name of file for all input
          sequences not matching any barcodes.
        """
        self.barcode_map = barcode_map
        self.primer = primer
        self._barcode_re = re.compile(r'TCAG(\w+){0}'.format(primer),
                                      re.IGNORECASE)
        self.dest_dir = dest_dir
        self.unmatched_dest = unmatched_dest
        self._handles = None

    def _open(self):
        """
        Opens files for each barcode, stores in handles
        """
        if self._handles:
            raise IOError("Already initialized with {0} handles".format(
                     len(self._handles)))

        for barcode, name in self.barcode_map.items():
            outpath = os.path.join(self.dest_dir, name)
            self._handles[barcode] = open(outpath, 'w')

        # Create a default handle
        default_outpath = os.path.join(self.dest_dir, self.unmatched_dest)
        self._handles[None] = open(default_outpath, 'w')

    def close(self):
        """
        Closes all open file handles
        """
        try:
            _close_all(self._handles.values())
            self._handles = None
        finally:
            self._default_handle.close()
            self._default_handle = None

    def _handle_record(record):
        """
        Identifies the barcode in the record, writes the record
        to the appropriate outfile, and returns the barcode
        """
        bases = record.bases_from_flows()
        m = self._barcode_re.match(bases)
        barcode = m.group(1) if m else None
        fp = self._handles[barcode]
        print >> fp, record
        return barcode

    def split(iterable):
        """
        Takes an iterable generating :ref:`anoisetools.sff.SFFRead`s,
        writes them to a set of output files
        returns a dictionary of mapping barcode -> # of reads
        """
        counts = collections.defaultdict(int)
        for record in iterable:
            barcode = self._handle_record(record)
            counts[barcode] += 1

        return counts
