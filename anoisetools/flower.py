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
    """
    Flower record
    """

    __slots__ = ['identifier', 'info', 'clip',
                 'flows', 'index', 'bases', 'quals']

    def __init__(self, identifier):
        self.identifier = identifier
        self.info = None
        self.clip = None
        self.flows = None
        self.index = None
        self.bases = None
        self.quals = None

    def __str__(self):
        """
        Return the record as a string
        """
        # TODO: make this a little less repetitive
        lines = []
        lines.append('>{}'.format(self.identifier))
        lines.append('  Info: \t{}'.format(self.info))
        lines.append('  Clip: \t{}'.format(' '.join(map(str, self.clip))))
        lines.append('  Flows:\t{}'.format(' '.join(['{:.2f}'.format(i)
                                                   for i in self.flows])))
        lines.append('  Index:\t{}'.format(' '.join(map(str, self.index))))
        lines.append('  Bases:\t{}'.format(self.bases))
        lines.append('  Quals:\t{}'.format(' '.join(map(str, self.quals))))
        return '\n'.join(lines)


# Matches read headers in flower output
_HEADER_REGEXP = re.compile(r'^\s*>(.*)')
# Matches data lines in flower output
_LINE_REGEXP = re.compile('^\s*(\w+):\s+(.*)$')

def _is_header(line):
    """
    True if the record is a FlowerRecord header, matching
    _HEADER_REGEXP.
    """
    return line and _HEADER_REGEXP.match(line)


def read_flower(iterable):
    """
    Reads an input containing lines from a flower sff.txt output,
    returning a generator yielding FlowerRecords
    """
    header = next(iterable)
    if not _is_header(header):
        raise ValueError("Invalid record identifier: {}".format(header))

    # Populate key handlers - each is a function taking a value,
    # Returning the value to set the lower case attribute on the record
    key_handlers = {}
    key_handlers['Info'] = lambda value: value
    key_handlers['Bases'] = key_handlers['Info']
    key_handlers['Clip'] = lambda value: map(int, value.split())
    key_handlers['Index'] = key_handlers['Clip']
    key_handlers['Quals'] = key_handlers['Clip']
    key_handlers['Flows'] = lambda value: map(float, value.split())

    record = FlowerRecord(_HEADER_REGEXP.match(header).group(1))

    while True:
        try:
            line = next(iterable)
        except StopIteration:
            yield record
            break
        if _is_header(line):
            yield record
            record = FlowerRecord(_HEADER_REGEXP.match(line).group(1))
        m = _LINE_REGEXP.match(line)
        if not m:
            # Ignore blank lines for now
            continue
        key, value = m.groups()

        try:
            if getattr(record, key.lower(), None) is not None:
                raise ValueError("{0} already set: {1}",
                                 getattr(record, key.lower()))
            else:
                setattr(record, key.lower(), key_handlers[key](value))
        except KeyError:
            raise ValueError("Unknown key: {}".format(key))
