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

from Bio import SeqIO

from anoisetools import sff, anoiseio


# Default minimum clean flows to keep a read
DEFAULT_MIN_FLOWS = 360
# Default maximum flows to keep
DEFAULT_MAX_FLOWS = 720


class RejectedSequenceError(ValueError):
    pass


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
    parser.add_argument('outname', metavar="OUTNAME",
            help="base name for output files - OUTNAME.fa and OUTNAME.dat")
    parser.add_argument('--input', metavar='INPUT', default=sys.stdin,
            type=argparse.FileType('r'),
            help='Input data file (default: stdin)')
    parser.add_argument('--save-rejected', metavar='FASTQ',
            dest="rej_handle",
            help='Save rejected files in FASTQ format',
            type=argparse.FileType('w'))
    return parser


# Possible TODO: Make the thresholds configurable
def is_flowgram_valid(flowgram, high_signal_cutoff=9.49,
                      low_signal_cutoff=0.7, signal_start=0.5):
    """
    Checks for a valid flowgram, defined as:
    * no reads above high_signal_cutoff
    * no reads in (signal_start, low_signal_cutoff)
    * At least one bp read
    """

    # Check arguments
    if len(flowgram) > 4:
        raise ValueError("Unexpected flowgram length: {0}".format(
            len(flowgram)))

    # Number of noisy reads
    noisy = 0
    # Number of items in flowgram with signal
    signal = 0

    for read in flowgram:
        if read > signal_start:
            # Positive reading (1+ nucleotides present)
            signal += 1

            # Check for noise
            if read < low_signal_cutoff or read > high_signal_cutoff:
                noisy += 1

    return noisy == 0 and signal > 0


def trim_noise(flows):
    """
    Trims a list of flows to the first noisy flowgram found
    """
    flowgram_size = 4

    # Iterate through the reads in groups of 4,
    # Making sure that there's
    # A) 1+ with signal above threshold
    # B) 0 with excessively high or low signal
    for i in xrange(0, len(flows), flowgram_size):
        flowgram = flows[i:i + flowgram_size]

        if not is_flowgram_valid(flowgram):
            # Return all data up to this point
            return flows[:i]

    # No invalid flowgrams found: return the whole sequence
    return flows


def handle_record(flows, primer_re, min_flows, max_flows):
    """
    Processes a record, checks for validity and returns the result

    primer_re must return a single group containing the primer
    and remainder of the sequence

    If the record has good data and meets the min_flows requirement:
    Returns a length-2 tuple containing the sequence to be
    written to the FASTA file, and a set of flows trimmed to
    max_flows length.

    Raises RejectedSequenceError on invalid sequence.
    """
    trimmed_flows = trim_noise(flows)
    trimmed = len(trimmed_flows) != len(flows)
    if min_flows is not None and len(trimmed_flows) < min_flows:
        raise RejectedSequenceError("Ins. Flows. Bad Flowgram: {0}".format(trimmed))

    trimmed_reading = sff.flow_to_seq(trimmed_flows)

    m = primer_re.match(trimmed_reading)

    if (min_flows is None or len(trimmed_flows) >= min_flows) and m:
        sequence = m.group(1)

        # Truncate the flow result to a maximum length
        # Note that the FASTA result is unchanged
        flow_result = trimmed_flows[:max_flows]

        return (sequence, flow_result)
    else:
        raise RejectedSequenceError("No primer match")


def invoke(reader, fa_handle, rej_handle, dat_path, , primer, min_flows,
           max_flows):
    """
    Main routine. Loop over the records in reader, write good results to
    fa_handle and dat_handle respectively.
    """
    good = 0
    bad = 0
    primer_re = re.compile(r'^{0}.*?({1}.*)'.format(sff.READ_BEGIN, primer))

    # On the first pass, write to a temporary file: we'll need to include the
    # record count and maximum length at the beginning of the file later.
    with tempfile.TemporaryFile() as dat_handle:
        for record in reader:
            try:
                sequence, flows = handle_record(record.flows, primer_re,
                                                min_flows, max_flows)
            except RejectedSequenceError:
                if rej_handle:
                    SeqIO.write([sequence], rej_handle, 'fastq')
            if sequence and flows:
                # Write FASTA
                print >> fa_handle, '>{0}\n{1}'.format(record.identifier,
                                                       sequence)
                # Write data
                print >> dat_handle, record.identifier, len(flows),
                print >> dat_handle, ' '.join(map('{:.2f}'.format, flows))
                good += 1
            else:
                bad += 1

        print dat_path, '{0} good; {1} bad '.format(good, bad)

        # Second pass - add a line with number clean and maximum length
        dat_handle.seek(0)
        with open(dat_path, 'w') as dat_fp:
            print >> dat_fp, good, max_flows
            # Copy from temp file to final file
            shutil.copyfileobj(dat_handle, dat_fp)


def main(parsed):
    """
    Check arguments, call invoke(...)
    """
    try:
        reader = anoiseio.AnoiseRawReader(parsed.input)
        with open(parsed.outname + '.fa', 'w') as fasta_handle:
            invoke(reader, fasta_handle, parsed.rej_handle, parsed.outname + '.dat',
                   parsed.primer, parsed.min_flows, parsed.max_flows)
    finally:
        parsed.input.close()
