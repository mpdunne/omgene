import itertools

from unittest import TestCase
from utils.blosum_matrices import blosum_dict, blosum_matrix


class TestBlosumMatrices(TestCase):

    def test_blosum_dict_has_all_aas(self):
        blos = blosum_dict()
        for a, b in itertools.product('ACBEDGFIHKMLNQPSRTWVYXZ*-', repeat=2):
            self.assertTrue((a, b) in blos)

    def test_blosum_dict_is_symmetric(self):
        blos = blosum_dict()
        for a, b in itertools.combinations('ACBEDGFIHKMLNQPSRTWVYXZ*-', 2):
            self.assertEqual(blos[(a, b)], blos[(b, a)])

    def test_blosum_matrix_is_nonzero(self):
        blos = blosum_dict()
        blos_mat = blosum_matrix(blos)
        self.assertTrue(not (blos_mat == 0).all())
