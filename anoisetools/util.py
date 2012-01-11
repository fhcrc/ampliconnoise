"""
Utility functions
"""

import re

# Map from Ambiguous Base to regex
_AMBIGUOUS_MAP = {
       'R': '[GA]',
       'Y': '[TC]',
       'K': '[GT]',
       'M': '[AC]',
       'S': '[GC]',
       'W': '[AT]',
       'B': '[GTC]',
       'D': '[GAT]',
       'H': '[ACT]',
       'V': '[GCA]',
       'N': '[AGCT]',
}

def ambiguous_pattern(sequence_str):
    """
    Convert an IUPAC ambiguous string to a regular expression pattern string
    """
    return ''.join(_AMBIGUOUS_MAP.get(c, c) for c in sequence_str)

def ambiguous_regex(sequence_str):
    """
    Convert an IUPAC ambiguous string to a regular expression
    """
    return re.compile(ambiguous_pattern(sequence_str), re.IGNORECASE)
