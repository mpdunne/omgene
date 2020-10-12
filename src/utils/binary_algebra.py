from collections import Counter


def bin_and(*items):
    """
    Returns the binary and of a list of binary strings.

    :return: The binary and-ed result.
    """
    r = items[0]
    for i in items:
        r = r & items
    return r


def bin_sub(a, b):
    """
    Whether one binary string (a) is a logical subset of the other (b).

    :return: (bool) True iff a is a binary subset of b.
    """
    return (a | b) == b


def bin_compat(bin1, bin2, checkpoints):
    """
    Returns whether or not the two input binaries, bin1 and bin2, agree
    for at least one digit of each of the checkpoints.
    """
    anded = bin1 & bin2
    return not any(anded & c == 0 for c in checkpoints)


def support(a, bases):
    """
    Return the number of bases of which a is a subset.
    """
    c = Counter(bases)
    return sum(c[j] for j in set(bases) if bin_sub(a, j))



