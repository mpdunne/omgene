import os
from unittest import TestCase

from utils.sequence_io import read_seqs, write_seqs, make_seq
from tempfile import mktemp


class TestSequenceIO(TestCase):

    def test_read_seqs(self):
        seqs = read_seqs("test/resources/aa.fasta")
        self.assertEqual(2, len(seqs))

    def test_write_seqs(self):
        tmp = mktemp(suffix=".fasta")
        write_seqs([make_seq('gene_1', 'AAAAAA')], tmp)
        reread = read_seqs(tmp)
        self.assertEqual(1, len(reread))
        os.remove(tmp)

    def test_make_seq(self):
        seq = make_seq('gene_1', 'AAAAAA')
        self.assertEqual('gene_1', seq.id)
        self.assertEqual('AAAAAA', seq.seq)
