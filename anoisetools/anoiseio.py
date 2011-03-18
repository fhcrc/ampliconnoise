"""
Tools for reading Ampiclonnoise outputs
"""

from anoisetools import flower

def splitkeys_reader(fp):
    """
    Reader for the output of ampiclonnoise SplitKeys*.pl scripts,
    which consist of a row with:
    0 key_barcode file name
    >id
    flow_count flow1 [flow2 [flow3...]]

    This function skips the first record, only returning the flow.

    Record length is truncated to flow_count
    """
    fp = (i.rstrip('\n') for i in fp)
    header = next(fp)
    while True:
        header = next(fp)[1:]
        split_line = next(fp).split()
        flow_length = int(split_line[0])
        flows = split_line[1:flow_length + 1]
        flows = map(float, flows)
        record = flower.FlowerRecord(header)
        record.flows = flows
        yield record
