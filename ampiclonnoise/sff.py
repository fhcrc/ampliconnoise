"""
Tools for working with raw flower output.
"""

import itertools
import math
import re


# Sequence in which nucleotides are flowed:
FLOW_ORDER = 'TACG'

# Sequence at the beginning of reads
READ_BEGIN = 'TCAG'


def flow_to_seq(flowgram, flow_order=FLOW_ORDER):
    """
    Takes a list of float flow results and transforms them to a sequence

    Taken from Chris Quince's Ampiclonnoise perl script
    """

    bases = itertools.cycle(flow_order)

    result = []
    for read_value, base in zip(flowgram, bases):
        # Round to the nearest value
        signal = int(math.floor(read_value + 0.5))
        result.extend([base] * signal)

    return ''.join(result)


class SFFRead(object):
    """
    SFF Reading
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
        lines.append('>{0}'.format(self.identifier))
        lines.append('  Info: \t{0}'.format(self.info))
        lines.append('  Clip: \t{0}'.format(' '.join(map(str, self.clip))))
        lines.append('  Flows:\t{0}'.format(self._flow_to_string()))
        lines.append('  Index:\t{0}'.format(' '.join(map(str, self.index))))
        lines.append('  Bases:\t{0}'.format(self.bases))
        lines.append('  Quals:\t{0}'.format(' '.join(map(str, self.quals))))
        return '\n'.join(lines)

    @property
    def right_clip(self):
        return self.clip[1]

    @property
    def left_clip(self):
        return self.clip[0]

    def bases_from_flows(self):
        """
        Return the bases associated with the flow records, trimmed
        to the right clip value
        """
        rtrim = self.right_clip
        return flow_to_seq(self.flows[:rtrim])

    def validate(self):
        """
        Makes sure all attributes are set
        """
        failures = []
        for i in SFFRead.__slots__:
            if getattr(self, i, None) is None:
                failures.append(i)
        if failures:
            raise ValueError("Missing attribute(s): {0}".format(
                ', '.join(failures)))

    @property
    def string_flows(self):
        """
        Returns the float flows as a list of properly formatted strings
        (2 decimal places)
        """
        return ['{0:.2f}'.format(i) for i in self.flows]

    def _flow_to_string(self):
        """
        Converts the flow floats to a string
        """
        return ' '.join(self.string_flows)

# Matches read headers in flower output
_HEADER_REGEXP = re.compile(r'^\s*>(.*)')

# Matches data lines in flower output
_LINE_REGEXP = re.compile('^\s*(\w+):\s+(.*)$')


def _is_header(line):
    """
    True if the record is a SFFRead header, matching
    _HEADER_REGEXP.
    """
    return line and _HEADER_REGEXP.match(line)


def parse_flower(iterable):
    """
    Reads an input containing lines from a flower sff.txt output,
    returning a generator of FlowerRecords contained in the input.
    """
    # Possible TODO: track line numbers for reporting errors?
    header = next(iterable)
    if not _is_header(header):
        raise ValueError("Invalid record identifier: {0}".format(header))

    # Populate key handlers - each is a function taking a value,
    # Returning the value to set the lower case attribute on the record
    # e.g. key_handlers['Info'] sets record.info when the Info: key is
    # is found
    key_handlers = {}
    key_handlers['Info'] = lambda value: value
    key_handlers['Bases'] = key_handlers['Info']
    key_handlers['Clip'] = lambda value: map(int, value.split())
    key_handlers['Index'] = key_handlers['Clip']
    key_handlers['Quals'] = key_handlers['Clip']
    key_handlers['Flows'] = lambda value: map(float, value.split())

    record = SFFRead(_HEADER_REGEXP.match(header).group(1))

    # Iterate over the input, generating records
    while True:
        try:
            line = next(iterable)
        except StopIteration:
            record.validate()
            yield record
            break
        if _is_header(line):
            # new record in input - validate and yield the current
            record.validate()
            yield record
            record = SFFRead(_HEADER_REGEXP.match(line).group(1))
        m = _LINE_REGEXP.match(line)
        if not m:
            # Ignore blank lines for now
            continue
        key, value = m.groups()

        try:
            if getattr(record, key.lower(), None) is not None:
                raise ValueError("{0} already set: {1}".format(key,
                        getattr(record, key.lower())))
            else:
                setattr(record, key.lower(), key_handlers[key](value))
        except KeyError:
            raise ValueError("Unknown key: {0}".format(key))
