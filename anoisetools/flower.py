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
        signal = int(read_value + 0.5)
        #signal = int(round(read_value, 0))
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
        lines.append('>{0}'.format(self.identifier))
        lines.append('  Info: \t{0}'.format(self.info))
        lines.append('  Clip: \t{0}'.format(' '.join(map(str, self.clip))))
        lines.append('  Flows:\t{0}'.format(self._flow_to_string()))
        lines.append('  Index:\t{0}'.format(' '.join(map(str, self.index))))
        lines.append('  Bases:\t{0}'.format(self.bases))
        lines.append('  Quals:\t{0}'.format(' '.join(map(str, self.quals))))
        return '\n'.join(lines)

    def _flow_to_string(self):
        """
        Converts the flow floats to a string
        """
        return ' '.join(('{:.2f}'.format(i) for i in self.flows))

    @property
    def right_trim(self):
        return self.clip[1]

    def to_anoise_raw(self):
        """
        Generates a string suitable for using as input to Ampiclonnoise,
        consisting of the identifier, a newline, the integer length of the
        flow, a space, and the float flow readings.
        """
        return '>{identifier}\n{length} {flow}'.format(
                identifier=self.identifier, length=len(self.flows),
                flow=self._flow_to_string())

    def validate(self):
        """
        Makes sure all attributes are set
        """
        failures = []
        for i in FlowerRecord.__slots__:
            if getattr(self, i, None) is None:
                failures.append(i)
        if failures:
            raise ValueError("Missing attribute(s): {0}".format(', '.join(failures)))


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

    record = FlowerRecord(_HEADER_REGEXP.match(header).group(1))

    # Iterate over the input, building records
    while True:
        try:
            line = next(iterable)
        except StopIteration:
            record.validate()
            yield record
            break
        if _is_header(line):
            record.validate()
            yield record
            record = FlowerRecord(_HEADER_REGEXP.match(line).group(1))
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
