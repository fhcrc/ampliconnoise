from cStringIO import StringIO
import unittest
from ampliconnoise import fastaio

class ParseFastaTestCase(unittest.TestCase):

    def setUp(self):
       fasta = """>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
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
       self.input_fasta = iter(fasta.splitlines())


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


