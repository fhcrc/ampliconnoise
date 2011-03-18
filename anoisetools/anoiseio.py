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

    This function skips the first record, only returning the flow
    """
    fp = (i.rstrip('\n') for i in fp)
    header = next(fp)
    while True:
        header = next(fp)[1:]
        flows = next(fp).split()[1:]
        flows = map(float, flows)
        record = flower.FlowerRecord(header)
        record.flows = flows
        yield record
