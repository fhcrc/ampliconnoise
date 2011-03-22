"""
Tools for reading Ampiclonnoise outputs
"""
import shutil
import tempfile

from ampiclonnoise import sff

class AnoiseRawReader(object):
    """
    Reader for Ampiclonnoise .raw files,
    which consist of a row with:
    record_count key_barcode file_name
    >id
    flow_count flow1 [flow2 [flow3...]]

    This function skips the first record, only returning the flow.

    Record length is truncated to flow_count
    """
    def __init__(self, fp):
        self._fp = (i.rstrip('\n') for i in fp)
        header = next(self._fp)
        count, barcode, fname = header.split(None, 2)
        self.count = int(count)
        self.barcode = barcode
        self.fname = fname

    def __iter__(self):
        while True:
            header = next(self._fp)[1:]
            split_line = next(self._fp).split()
            flow_length = int(split_line[0])
            flows = split_line[1:flow_length + 1]
            flows = map(float, flows)
            record = sff.SFFRead(header)
            record.flows = flows
            yield record

def _record_to_anoise_raw(record):
    """
    Generates a string suitable for using as input to Ampiclonnoise,
    consisting of the identifier, a newline, the integer length of the
    flow, a space, and the float flow readings.
    """
    return '>{identifier}\n{length} {flow}'.format(
            identifier=record.identifier, length=record.right_clip,
            flow=record._flow_to_string())

class AnoiseRawWriter(object):
    """
    Writer for AmpiclonNoise .raw files

    .raw files require a row at the top with a count of records.
    To get around that, this class writes to a temporary file, then
    copies the results to fp when closed.
    """

    def __init__(self, fp, identifier):
        self.identifier = identifier
        self._fp = fp
        self.file_name = getattr(fp, 'name', "Unknown File")
        self._temp = tempfile.TemporaryFile()

        self.count = 0

    def _header(self):
        return '{0.count} {0.identifier} {0.file_name}'.format(self)

    def write(self, record):
        print >>self._temp, _record_to_anoise_raw(record)
        self.count += 1

    def _copy_to_output(self):
        """
        Add a header to _fp, writes all the records in _temp
        """
        # Rewind
        self._temp.seek(0)
        # Write header
        print >>self._fp, self._header()
        # Copy contents
        shutil.copyfileobj(self._temp, self._fp)

    def close(self):
        """
        Close the writer

        Copies all records into the destination file, closes the temp handle
        """
        # No action if already closed
        if self._temp.closed and self._fp.closed:
            return
        try:
            self._copy_to_output()
        finally:
            self._temp.close()
            self._fp.close()
