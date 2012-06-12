"""
Tools for reading Ampliconnoise outputs
"""
import itertools
import os.path
import re
import shutil
import tempfile

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from anoisetools import sff


class AnoiseRawReader(object):
    """
    Reader for AmpliconNoise .raw files,
    which consist of a row with:
    record_count key_barcode file_name
    >id
    flow_count flow1 [flow2 [flow3...]]

    This function skips the first record, only returning the flow.

    Yields BioPython SeqRecords
    """
    def __init__(self, fp):
        self._fp = (i.rstrip('\n') for i in fp)
        header = next(self._fp)
        count, barcode, fname = header.split(None, 2)
        self.count = int(count)
        self.barcode = barcode
        self.fname = fname

    def __iter__(self):
        while True:
            header = next(self._fp)[1:]
            split_line = next(self._fp).split()
            flow_length = int(split_line[0])
            flows = split_line[1:]
            flows = map(float, flows)
            bases = sff.flow_to_seq(flows)
            sequence = Seq(bases, generic_dna)
            record = SeqRecord(sequence, id=header)
            # Add flows and clip as annotations
            record.annotations['flow_values'] = [int(i * 100) for i in flows]
            record.annotations['clip_flow_right'] = flow_length
            yield record


def _record_to_anoise_raw(seq_record):
    """
    Generates a string suitable for using as input to AmpliconNoise,
    consisting of the identifier, a newline, the integer length of the
    flow, a space, and the float flow readings.
    """
    l, r = sff.find_clip(seq_record)
    return '>{identifier}\n{length} {flow}'.format(
            identifier=seq_record.id,
            length=r,
            flow=sff.flows_to_string(seq_record))


class AnoiseRawWriter(object):
    """
    Writer for AmpliconNoise .raw files

    .raw files require a row at the top with a count of records.
    To get around that, this class writes to a temporary file, then
    copies the results to fp when closed.
    """

    def __init__(self, fp, identifier):
        self.identifier = identifier
        self._fp = fp
        self.file_name = getattr(fp, 'name', "Unknown File")
        self._temp = tempfile.TemporaryFile()

        self.count = 0

    def _header(self):
        return '{0.count} {0.identifier} {0.file_name}'.format(self)

    def write(self, record):
        print >> self._temp, _record_to_anoise_raw(record)
        self.count += 1

    def write_records(self, records):
        for record in records:
            self.write(record)

    def _copy_to_output(self):
        """
        Add a header to _fp, writes all the records in _temp
        """
        # Rewind
        self._temp.seek(0)
        # Write header
        print >> self._fp, self._header()
        # Copy contents
        shutil.copyfileobj(self._temp, self._fp)

    def close(self):
        """
        Close the writer

        Copies all records into the destination file, closes the temp handle
        """
        # No action if already closed
        if self._temp.closed and self._fp.closed:
            return
        try:
            self._copy_to_output()
        finally:
            self._temp.close()
            self._fp.close()

def _read_count(read_name):
    m = re.search(r'noise_\d+_(\d+)$', read_name)
    if m:
        return int(m.group(1))
    return 1

def read_mapping(fp):
    """
    Read an ampliconnoise .mapping file
    """
    # Base name
    bp = os.path.splitext(os.path.basename(fp.name))[0]

    indexes = itertools.count()

    lines = (i.strip() for i in fp)
    for line in lines:
        used_sequence, seqs = line.split(None, 1)
        seqs = seqs.split(',')
        total_weight = sum(_read_count(i) for i in seqs)
        read_name = '{0}_{1}_{2}'.format(bp, next(indexes), total_weight)
        yield read_name, seqs

def merge_mapping(snoise_map, pnoise_map):
    """
    Merge two mapping files, from SeqNoise and PyroNoise, ending with a mapping
    from seqnoise_id -> list of original sequence ids
    """
    pnoise_map = dict(pnoise_map)
    for i, seqs in snoise_map:
        orig_seqs = [o for s in seqs for o in pnoise_map[s]]
        yield i, orig_seqs
