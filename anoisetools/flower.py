"""
Tools for working with raw flower output.
"""

import itertools
import re

# Sequence in which nucleotides are flowed:
FLOW_ORDER = 'TACG'


def flow_to_seq(flowgram, flow_order=FLOW_ORDER):
    """
    Takes a list of float flow results and transforms them to a sequence

    Taken from Chris Quince's Ampiclonnoise perl script
    """

    bases = itertools.cycle(flow_order)

    result = []
    for read_value, base in zip(flowgram, bases):
        # Round to the nearest value
        signal = int(round(read_value, 0))
        result.extend([base] * signal)

    return ''.join(result)


class FlowerRecord(object):
    def __init__(self, identifier):
        self.identifier = identifier
        self.info = None
        self.clip = None
        self.flows = None
        self.index = None
        self.bases = None
        self.quals = None


#>FTWCYXX01BTPDQ
#  Info:   2009-04-08 11:31:35 R1 (631,1068)
#  Clip:   5 282
#  Flows:  1.04 0.01 1.02 0.07 0.05 0.97 0.06 2.00 0.96 1.09 0.98 0.08 0.02 1.18   Index:  1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 23 23 23 25 27 27 27 3  Bases:  tcagGTACAGAAAAATTCCCCTCCCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT  Quals:  37 37 37 37 37 37 37 37 38 34 20 20 20 20 15 34 34 37 37 37 37 38 34 3>FTWCYXX01B0TTO
#  Info:   2009-04-08 11:31:35 R1 (712,1642)
#  Clip:   5 278
#  Flows:  1.04 0.04 1.01 0.09 0.06 0.99 0.09 1.97 1.02 1.06 0.98 0.07 0.01 1.20
#  Index:  1 3 6 8 8 9 10 11 14 16 18 18 18 18 18 21 21 23 25 27 27 29 31 33 35 3  Bases:  tcagGTACAGAAAAATTCTCCTCTCAATTAAAACTATGTGATGTGATTTCTATGTCCCCTCCTGAGGGTT  Quals:  37 37 37 37 37 37 37 37 38 38 28 28 28 28 28 38 38 37 37 37 37 37 40
#

_HEADER_REGEXP = re.compile(r'^\s*>(.*)')
def _is_header(line):
    return line and _HEADER_REGEXP.match(line)

_LINE_REGEXP = re.compile('^\s*(\w+):\s+(.*)$')

def read_flower(iterable):
    header = next(iterable)
    if not _is_header(header):
        raise ValueError("Invalid record identifier: {}".format(header))

    record = FlowerRecord(header[1:])

    while True:
        try:
            line = next(iterable)
        except StopIteration:
            yield record
            raise StopIteration()
        if _is_header(line):
            yield record
            record = FlowerRecord(_HEADER_REGEXP.match(line).group(1))
        m = _LINE_REGEXP.match(line)
        if not m:
            # Ignore blank lines for now
            continue
        key, value = m.groups()

        if key in ('Info', 'Bases'):
            setattr(record, key.lower(), value)
        elif key in ('Clip', 'Index', 'Quals'):
            setattr(record, key.lower(), map(int, value.split()))
        elif key == 'Flows':
            record.flows = map(float, value.split())
        else:
            raise ValueError("Unknown key: {}".format(key))
