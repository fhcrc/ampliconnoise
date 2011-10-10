"""
Clean sff.txt or anoise raw data file:

    * Truncate each sequence to the first occurrence of bad data, defined by:
        * 4 flows in a row with no reads
        * Reading above some high signal cutoff
        * Positive, but excessively low readings
    * Trim long sequences to a defined maximum length
    * Only accept sequences with a defined minimum length

Creating an output .fasta and a .dat file, suitable for feeding into the
AmpliconNoise C code.

Per Quince et al, 2011,
* For GS FLX runs, require min_flows of 360, max_flows of 360
* For Titanium runs, require min_flows of 360, max_flows of 720
"""

import argparse
import re
import shutil
import sys
import tempfile

from anoisetools import sff, anoiseio


# Default minimum clean flows to keep a read
DEFAULT_MIN_FLOWS = 360
# Default maximum flows to keep
DEFAULT_MAX_FLOWS = 720

FLOWGRAM_SIZE = 4

def build_parser(subparsers):
    """
    Given an ArgumentParser, adds options taken by ``clean``
    """
    parser = subparsers.add_parser("clean",
            help='Clean flowgrams prior to running AmpliconNoise')
    parser.description = """Clean flowgrams prior to
processing with AmpliconNoise: enforce minimum length, remove invalid data,
trim."""
    parser.add_argument('--min-flows', type=int,
            help="Minimum length to accept sequence"
            " default: %(default)s", default=DEFAULT_MIN_FLOWS)
    parser.add_argument('--max-flows', type=int,
            help="Maximum length to trim sequences to "
            "default: %(default)s", default=DEFAULT_MAX_FLOWS)
    parser.add_argument('primer', metavar='PRIMER',
            help="Regexp to identify primer sequence")
    parser.add_argument('--max-empty',
            type=int, default=0, help="""Maximum empty flowgrams (default: %(default)s)""")
    parser.add_argument('outname', metavar="OUTNAME",
            help="base name for output files - OUTNAME.fa and OUTNAME.dat")
    parser.add_argument('--input', metavar='INPUT', default=sys.stdin,
            type=argparse.FileType('r'),
            help='Input data file (default: stdin)')
    return parser


class FlowgramFilter(object):
    """
    Instrumented flowgram filter_record
    """

    def __init__(self, primer_re=None, high_signal_cutoff=9.49,
            low_signal_cutoff=0.7, signal_start=0.5, max_empty=0, min_flows =
            DEFAULT_MIN_FLOWS, max_flows=DEFAULT_MAX_FLOWS):
        self.high_signal_cutoff = high_signal_cutoff
        self.low_signal_cutoff = low_signal_cutoff
        self.signal_start = signal_start

        self.primer_re = primer_re or re.compile('.')
        self.min_flows = min_flows
        self.max_flows = max_flows
        self.max_empty = max_empty

        self.empty_flow = 0
        self.high_signal = 0
        self.ambiguous_flow = 0
        self.passed = 0
        self.failed = 0

    def is_empty(self, flowgram):
        """
        Returns true if the flowgram has no signal
        """
        r = all(i < self.signal_start for i in flowgram)
        if r:
            self.empty_flow += 1
        return r

    def is_high(self, flowgram):
        """
        Returns true if the flowgram contains any excessively high values
        """
        r = any(i > self.high_signal_cutoff for i in flowgram)
        if r:
            self.high_signal += 1
        return r

    def is_ambiguous(self, flowgram):
        """
        Returns true if the flowgram contains any ambiguous values
        """
        r =  any(i > self.signal_start and i < self.low_signal_cutoff
                   for i in flowgram)
        if r:
            self.ambiguous_flow += 1
        return r

    def is_valid(self, flowgram):
        if not len(flowgram) <= FLOWGRAM_SIZE:
            raise ValueError("Unexpected flowgram length: {0}".format(
                len(flowgram)))

        return not (self.is_empty(flowgram) or self.is_high(flowgram) or
                    self.is_ambiguous(flowgram))

    def _clean_flow_iter(self, flows):
        """
        Generator of cleaned flowgrams
        """
        empty_count = 0

        flowgrams = (flows[i:i + FLOWGRAM_SIZE]
                     for i in xrange(0, len(flows), FLOWGRAM_SIZE))
        for flowgram in flowgrams:
            if self.is_empty(flowgram):
                self.empty_flow += 1
                empty_count += 1

                # If maximum empty flowgram count exceeded, stop
                if empty_count > self.max_empty:
                    break
                else:
                    continue

            # Handle high and ambiguous flowgrams
            # For now, just truncate when seen
            if not self.is_high(flowgram) and not self.is_ambiguous(flowgram):
                yield flowgram
            else:
                break

    def filter_record(self, flows):
        return [flow for flowgram in self._clean_flow_iter(flows)
                for flow in flowgram]

    def filter_records(self, records):
        for record in records:
            flows = record.annotations['flow_values'][:record.annotations['clip_qual_right']]
            flows = [float(i) / 100 for i in flows]
            flows = self.filter_record(flows)
            trimmed_reading = sff.flow_to_seq(flows)

            m = self.primer_re.match(trimmed_reading)
            if (self.min_flows is None or len(flows) >= self.min_flows) and m:
                sequence = m.group(1)

                # Truncate the flow result to a maximum length
                # Note that the FASTA result is unchanged
                flow_result = flows[:self.max_flows]

                self.passed += 1
                yield (record.id, sequence, flow_result)
            else:
                self.failed += 1


def invoke(reader, fa_handle, dat_path, primer, min_flows, max_flows,
        max_empty):
    """
    Main routine. Loop over the records in reader, write good results to
    fa_handle and dat_handle respectively.
    """
    primer_re = re.compile(r'^{0}.*?({1}.*)'.format(sff.READ_BEGIN, primer))
    flowgram_filter = FlowgramFilter(primer_re, min_flows=min_flows,
            max_flows=max_flows, max_empty=max_empty)

    cleaned_iter = flowgram_filter.filter_records(reader)

    # On the first pass, write to a temporary file: we'll need to include the
    # record count and maximum length at the beginning of the file later.
    with tempfile.TemporaryFile() as dat_handle:
        for identifier, sequence, flows in cleaned_iter:
            # Write FASTA
            print >> fa_handle, '>{0}\n{1}'.format(identifier,
                                                   sequence)
            # Write data
            print >> dat_handle, identifier, len(flows),
            print >> dat_handle, ' '.join(map('{:.2f}'.format, flows))

        print dat_path, '{0} passed; {1} failed'.format(flowgram_filter.passed,
                flowgram_filter.failed)

        # Second pass - add a line with number clean and maximum length
        dat_handle.seek(0)
        with open(dat_path, 'w') as dat_fp:
            print >> dat_fp, flowgram_filter.passed, max_flows
            # Copy from temp file to final file
            shutil.copyfileobj(dat_handle, dat_fp)

    return flowgram_filter.passed, flowgram_filter.failed


def main(parsed):
    """
    Check arguments, call invoke(...)
    """
    try:
        reader = anoiseio.AnoiseRawReader(parsed.input)
        with open(parsed.outname + '.fa', 'w') as fasta_handle:
            invoke(reader, fasta_handle, parsed.outname + '.dat',
                   parsed.primer, parsed.min_flows, parsed.max_flows,
                   parsed.max_empty)
    finally:
        parsed.input.close()
