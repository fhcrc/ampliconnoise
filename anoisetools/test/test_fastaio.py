from cStringIO import StringIO
import unittest
from anoisetools import fastaio


FASTA_SEQS = """>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY
>hgi|55211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLV
EWIWGGFSVDKASLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLG
LLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVIL
GLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGX
IENY"""

class ParseFastaTestCase(unittest.TestCase):

    def setUp(self):
       self.input_fasta = iter(FASTA_SEQS.splitlines())

    def tearDown(self):
        pass

    def test_read_type(self):
        reader = fastaio.parse_fasta(self.input_fasta)
        first = next(reader)
        self.assertEquals(fastaio.Sequence, type(first))

    def test_read_id(self):
        reader = fastaio.parse_fasta(self.input_fasta)
        first = next(reader)
        self.assertEquals("gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]", first.id)

    def test_read_seq(self):
        reader = fastaio.parse_fasta(self.input_fasta)
        first = next(reader)
        self.assertEquals("""LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGTNLVEWIWGGFSVDKATLNRFFAFHFILPFTMVALAGVHLTFLHETGSNNPLGLTSDSDKIPFHPYYTIKDFLGLLILILLLLLLALLSPDMLGDPDNHMPADPLNTPLHIKPEWYFLFAYAILRSVPNKLGGVLALFLSIVILGLMPFLHTSKHRSMMLRPLSQALFWTLTMDLLTLTWIGSQPVEYPYTIIGQMASILYFSIILAFLPIAGXIENY""", first.seq)

    def test_read_multiple(self):
        reader = fastaio.parse_fasta(self.input_fasta)
        seqs = [i for i in reader]
        self.assertEquals(2, len(seqs))
        second = seqs[1]
        self.assertEquals("hgi|55211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]", second.id)


class WriteFastaTestCase(unittest.TestCase):

    def setUp(self):
        self.fastas = [fastaio.Sequence('test sequence 1', 'ACTG'),
                       fastaio.Sequence('test sequence 2', 'TTGCAGT')]

        self.output = StringIO()

    def tearDown(self):
        pass

    def test_write_basic(self):
        fastaio.write_fasta(self.fastas, self.output, wrap=None)

        actual = self.output.getvalue()
        self.assertEquals('>test sequence 1\nACTG\n>test sequence 2\nTTGCAGT\n',
                actual)

    def test_wrap_true(self):
        """wrap=True should wrap to 80 char/line"""
        fastas = [fastaio.Sequence('test1', 'ACTG'*40)]
        fastaio.write_fasta(fastas, self.output, wrap=True)
        expected = """>test1
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
"""
        self.assertEquals(expected, self.output.getvalue())
