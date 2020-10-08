import numpy as np
from Bio.Align import substitution_matrices


amino_acids = 'ACBEDGFIHKMLNQPSRTWVYXZ*'


def blosum_dict(single_gap_penalty=-1, double_gap_penalty=0.5):
    """
    Returns a fully symmetric version of the blosum matrix, with gap penalties included.

    :param single_gap_penalty: (float) The penalty for a single gap
    :param double_gap_penalty: (float) The penalty for a double gap
    """

    blos = substitution_matrices.load('BLOSUM62')
    blos_dict = {}

    for a in amino_acids:

        for b in amino_acids:
            blos_dict[(b, a)] = blos_dict[(a, b)] = blos[(a, b)]

        blos_dict[(a, '-')] = blos_dict[('-', a)] = single_gap_penalty

    blos_dict[('-', '-')] = double_gap_penalty
    return blos_dict


def blosum_matrix(blosdict):
    """
    Convert a blosum dict into a blosum matrix.

    :param blosdict: (Dict[Tuple, float]) A dict of substitution distances for pairs of AAs.
    :return: A matrix of pairwise substitution distances.
    """
    blosmat = np.zeros([len(amino_acids), len(amino_acids)])
    for i, e in enumerate(amino_acids):
        for j, f in enumerate(amino_acids):
            blosmat[i, j] = blosdict[(e, f)]
    return blosmat


def blosum_positions():
    """
    A dict of the position in the substitution matrix for each AA.

    :return: (Dict) The position in the substitution matrix for each AA.
    """
    return {(e, i) for i, e in enumerate(amino_acids)}
