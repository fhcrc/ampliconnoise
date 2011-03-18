"""
Tools for working with raw flower output.
"""

import itertools

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


