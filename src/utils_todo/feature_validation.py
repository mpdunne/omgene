# These get set on entry.
static_do = ["gc", "gt"]
static_ac = ["ag"]
static_sc = ["atg"]


def is_start_codon(c):
    return is_feature(c, 3, static_sc)


def is_stop_codon(c):
    return is_feature(c, 3, ["tga", "taa", "tag"])


def is_donor(c):
    return is_feature(c, 2, static_do)


def is_acceptor(c):
    return is_feature(c, 2, static_ac)


def is_donor_site(p, cds):
    return is_donor(cds[p:p + 2])


def is_acceptor_site(p, cds):
    return is_acceptor(cds[p - 3:p - 1])


def is_feature(c, length, choices):
    if len(c) != length:
        return False
    return c.lower() in choices
