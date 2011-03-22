"""
Clean sff.txt or anoise raw data file:

    * Truncate each sequence to the first occurrence of bad data, defined by:
        * 4 flows in a row with no reads
        * Reading above some high signal cutoff
        * Positive, but excessively low readings
    * Trim long sequences to a defined maximum length
    * Only accept sequences with a defined minimum length

Creating an output .fasta and a .dat file, suitable for feeding into the
AmpiclonNoise C code.
"""

import argparse
import re
import shutil
import sys
import tempfile

from ampiclonnoise import sff, anoiseio

DEFAULT_MIN_FLOWS = None
DEFAULT_MAX_FLOWS = 360

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

    Otherwise, returns ``(None, None)``
    """
    flow_length = len(flows)
    flowgram_size = 4

    trimmed_flows = trim_noise(flows)
    trimmed_reading = sff.flow_to_seq(trimmed_flows)

    m = primer_re.match(trimmed_reading)
    if (min_flows is None or len(trimmed_flows) >= min_flows) and m:
        sequence = m.group(1)

        # Truncate the flow result to a maximum length
        # Note that the FASTA result is unchanged
        flow_result = trimmed_flows[:max_flows]

        return (sequence, flow_result)
    else:
        return (None, None)


def invoke(reader, fa_handle, dat_path, primer, min_flows, max_flows):
    """
    Main routine. Loop over the records in reader, write good results to
    fa_handle and dat_handle respectively.
    """
    good = 0
    bad = 0
    primer_re = re.compile(r'^TCAG.*({0}.*)'.format(primer))

    # On the first pass, write to a temporary file: we'll need to include the
    # record count and maximum length at the beginning of the file later.
    with tempfile.TemporaryFile() as dat_handle:
        for record in reader:
            sequence, flows = handle_record(record.flows, primer_re,
                                            min_flows, max_flows)
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

        print '{0} good records; {1} bad records'.format(good, bad)

        # Second pass - add a line with number clean and maximum length
        dat_handle.seek(0)
        with open(dat_path, 'w') as dat_fp:
            print >> dat_fp, good, max_flows
            # Copy from temp file to final file
            shutil.copyfileobj(dat_handle, dat_fp)


def main(args=sys.argv[1:]):
    """
    Check arguments, call invoke(...)
    """
    parser = argparse.ArgumentParser(description="""Clean flowgrams prior to
processing with AmpiclonNoise: enforce minimum length, remove invalid data,
trim.""")
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
    parser.add_argument('--flower-input', action='store_true',
            help='Input is in flower-decoded .sff.txt format '
                 '(default: %(default)s)',
            default=False)
    parsed = parser.parse_args(args)

    try:
        if parsed.flower_input:
            reader = sff.parse_flower(parsed.input)
        else:
            reader = anoiseio.SplitKeysReader(parsed.input)
        with open(parsed.outname + '.fa', 'w') as fasta_handle:
            invoke(reader, fasta_handle, parsed.outname + '.dat',
                   parsed.primer, parsed.min_flows, parsed.max_flows)
    finally:
        parsed.input.close()
