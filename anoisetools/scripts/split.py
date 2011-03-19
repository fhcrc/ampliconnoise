
# Class
#  - Read barcode file
#  - Generate list of file descriptors
#  - Read flower input
#  - Find a destination based on barcode and primer
#  - Write results to appropriate output file

import csv

def load_barcodes(fp):
    reader = csv.reader(fp)
    d = {}
    for line in reader:
        d[line[1]] = line[0]
    return d


class SFFFlowerSplitter(object)u:

    def __init__(self, barcode_map, primer, dest_dir):
        """
        Initialize a new instance

        :param map barcode_map: map from sequence barcode to name
        :param str primer: Primer sequence
        :param str dest_dir: destination directory for output
        """
        self.barcode_map = barcode_map
        self.primer = primer
        self.dest_dir = dest_dir
        self._handles = None

    def _open(self):
        """
        Opens files for each barcode, stores in handles
        """
        if self._handles:
            raise IOError("Already initialized with {0} handles".format(
                     len(self._handles)))

        for barcode, name in self.barcode_map.items():
            outpath = os.path.join(dest_dir, name)
            self._handles[barcode] = open(outpath, 'w')

    def close(self):
        """
        Closes all open file handles
        """
        for fh in self._handles.values():
            fh.close()

    def handle_record(record):
        raise NotImplementedError("TODO: handle records")
