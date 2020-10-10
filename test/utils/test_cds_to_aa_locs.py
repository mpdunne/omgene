from src.utils_todo.alignment import cds_to_aa_locs

from unittest import TestCase


class TestBlosumMatrices(TestCase):

    def test_alignment_score(self):
        pass

    def test_cds_to_aa_locs(self):
        check1 = cds_to_aa_locs([0, 12, 18, 24])
        self.assertListEqual([0, 4, 6, 8], check1)

        check2 = cds_to_aa_locs([0, 11, 18, 24])
        self.assertListEqual([0, 3, 4, 6, 8], check2)