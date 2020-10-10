global_donor_sites = []
global_acceptor_sites = []
global_start_codons = []


def initialise_features(donor_sites=("gc", "gt"), acceptor_sites=("ag"), start_codons=("atg")):
    """
    Initialise the global set of codon / splice site feature patterns.
    """
    global global_donor_sites, global_acceptor_sites, global_start_codons

    global_donor_sites = donor_sites
    global_acceptor_sites = acceptor_sites
    global_start_codons = start_codons


def is_start_codon(c):
    """
    Is it a start codon?
    """
    return is_feature(c, global_start_codons)


def is_stop_codon(c):
    """
    Is it a stop codon?
    """
    return is_feature(c, ["tga", "taa", "tag"])


def is_donor(c):
    """
    Is it a donor codon?
    """
    return is_feature(c, global_donor_sites)


def is_acceptor(c):
    """
    Is it an acceptor codon?
    """
    return is_feature(c, global_acceptor_sites)


def is_donor_site(p, cds):
    """
    Is there a donor in the CDS at the given position?
    """
    return is_donor(cds[p:p + 2])


def is_acceptor_site(p, cds):
    """
    Is there an acceptor in the CDS at the given position?
    """
    return is_acceptor(cds[p - 3:p - 1])


def is_feature(c, choices):
    """
    Is the provided string the appropriate length and does it appear in the list?
    """
    return c.lower() in choices
