"""
Tools for reading Ampiclonnoise outputs
"""

from anoisetools import flower

class SplitKeysReader(object):
    """
    Reader for the output of ampiclonnoise SplitKeys*.pl scripts,
    which consist of a row with:
    0 key_barcode file name
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
            record = flower.FlowerRecord(header)
            record.flows = flows
            yield record

