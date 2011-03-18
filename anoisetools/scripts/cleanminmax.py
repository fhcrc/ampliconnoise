"""

"""

import argparse
import re
import sys

from anoisetools import flower

DEFAULT_MIN_LENGTH = 400
DEFAULT_MAX_LENGTH = 600
MAX_ANOISE_LENGTH = 360

def is_flowgram_valid(flowgram, high_signal_cutoff=9.49,
                      low_signal_cutoff=0.7, signal_start=0.5):
    """
    Checks for a valid flowgram, defined as:
    * no reads above high_signal_cutoff
    * no reads in (signal_start, low_signal_cutoff)
    * At least one reading
    """

    # Check arguments
    if len(flowgram) > 4:
        raise ValueError("Unexpected flowgram length: {0}".format(len(flowgram)))

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

def handle_record(flows, primer_re, min_length, max_length):
    """
    Processes a record, checks for validity and returns the result

    primer_re must return a single group containing the primer
    and remainder of the sequence

    If the record has good data and meets the min_length requirement:
    Returns a length-2 tuple containing the sequence to be
    written to the FASTA file, and a set of flows trimmed to
    MAX_ANOISE_LENGTH.

    Otherwise, returns ``(None, None)``
    """
    flow_length = len(flows)
    flowgram_size = 4

    # Iterate through the reads in groups of 4,
    # Making sure that there's
    # A) 1+ with signal above threshold
    # B) 0 with excessively high or low signal
    trimmed_length = flow_length
    for i in xrange(0, flow_length, flowgram_size):
        flowgram = flows[i:i + flowgram_size]
        if not is_flowgram_valid(flowgram):
            trimmed_length = i

    trimmed_flows = flows[:trimmed_length]
    trimmed_reading = flower.flow_to_seq(trimmed_flows)

    m = primer_re.match(trimmed_reading)
    if trimmed_length >= min_length and m:
        sequence = m.group(1)

        flow_result = trimmed_flows[:MAX_ANOISE_LENGTH]

        return (sequence, flow_result)
    else:
        print 'failed', m, trimmed_length
        print trimmed_reading,
        return (None, None)


def invoke(reader, fa_handle, dat_handle, primer, min_length, max_length):
    """
    Main routine. Loop over the records in reader, write good results to
    fa_handle and dat_handle respectively.
    """
    good = 0
    bad = 0
    primer_re = re.compile(r'^TCAG.*({0}.*)'.format(primer))
    for record in reader:
        sequence, flows = handle_record(record.flows, primer_re, min_length,
                                        max_length)
        if sequence and flows:
            # Write FASTA
            print >> fa_handle, '>{0}\n{1}'.format(record.identifier, sequence)
            # Write data
            print >> dat_handle, record.identifier, len(flows),
            print >> dat_handle, ' '.join(map(str, flows))
            good += 1
        else:
            bad += 1

    print '{0} good records, {1} bad records'.format(good, bad)


def main(args=sys.argv[1:]):
    """
    Check arguments, call invoke(...)
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--min-length', type=int,
            help="Minimum length to accept sequence"
            " default: %(default)d", default=DEFAULT_MIN_LENGTH)
    parser.add_argument('--max-length', type=int,
            help="Maximum length to trim sequences to "
            "default: %(default)d", default=DEFAULT_MAX_LENGTH)
    parser.add_argument('primer', metavar='PRIMER',
            help="Regexp to identify primer sequence")
    parser.add_argument('outname', metavar="OUTNAME",
            help="base name for output files - OUTNAME.fa and OUTNAME.dat")
    parser.add_argument('input', metavar='INPUT', default=sys.stdin,
            type=argparse.FileType('r'),
            help='Input flower-processed data file (default: stdin)')
    parsed = parser.parse_args(args)

    try:
        reader = flower.read_flower(parsed.input)
        with open(parsed.outname + '.fa', 'w') as fasta_handle:
            with open(parsed.outname + '.dat', 'w') as dat_handle:
                invoke(parsed.input, fasta_handle, dat_handle,
                       parsed.primer, parsed.min_length, parsed.max_length)
    finally:
        parsed.input.close()

