"""
Tools for working with raw SFF files
"""

import itertools


# Complement of sequence in which nucleotides are flowed.
# THe standard bases read would be ATGC
FLOW_ORDER = 'TACG'

# Sequence at the beginning of reads
READ_BEGIN = 'TCAG'

def cumsum(iterable):
    """
    Generates a cumulative sum from iterable
    """
    s = 0
    for i in iterable:
        s += i
        yield s

def find_clip(seqrecord):
    """
    Determine flowgram clipping locations based on base clipping locations
    """
    a = seqrecord.annotations
    flow_index = list(cumsum(a['flow_index']))
    left_flow = flow_index[a['clip_qual_left']]
    right_flow = flow_index[a['clip_qual_right'] - 1]
    return left_flow, right_flow

def flow_to_seq(flowgram, flow_order=FLOW_ORDER):
    """
    Takes a list of float flow results and transforms them to a sequence

    Taken from Chris Quince's AmpliconNoise perl scripts
    """
    # Bases are repeatedly flowed in flow_order
    bases = itertools.cycle(flow_order)

    result = []

    # Build the sequence
    for read_value, base in zip(flowgram, bases):
        # Get base counts: round to the nearest integer value
        if isinstance(read_value, float):
            signal = int(round(read_value, 0))
        else:
            # Integer float values from BioPython -
            # float value * 100
            signal = int(round(read_value, -2) / 100)
        result.extend([base] * signal)

    return ''.join(result)


def bases_from_flows(seq_record):
    """
    Return the bases associated with the flow records, trimmed
    to the right clip value
    """
    right_clip = find_clip(seq_record)[1]
    return flow_to_seq(seq_record.annotations['flow_values'][:right_clip])


def flows_to_string(seq_record):
    """
    Return the flows
    """
    flows = seq_record.annotations['flow_values']
    flows = [float(i) / 100 for i in flows]
    string_flows = ['{0:.2f}'.format(i) for i in flows]
    return ' '.join(string_flows)
