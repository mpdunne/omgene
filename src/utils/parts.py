from utils.feature_validation import *
from utils_todo.gtf_tools import  *


def is_donor_part(gtfline, cds):
    """
    Whether this gtf line could possibly contain a donor in the correct place.
    """
    if is_terminal_part(gtfline, cds):
        return False

    if is_donor(cds[gtfline.end: gtfline.end + 2]):
        return True


def is_acceptor_part(gtfline, cds):
    """
    Whether this gtf line could possibly contain an acceptor in the correct place.
    """
    return is_acceptor(cds[gtfline.start - 3: gtfline.start - 1])


def is_initial_part(gtfline, cds):
    """
    Whether this gtf line starts with a start codon and is in the right frame.
    """
    if gtfline.frame != 0:
        return False

    part = get_cds_live([gtfline], cds)

    # If it's a micromicropart, can't tell at this stage if it's a start. return true.
    if len(part) < 3:
        return True

    return is_start_codon(part[0:3])


def is_terminal_part(gtfline, cds):
    """
    Whether this gtf line ends with a stop codon.
    """
    # Must be zero-ended
    endframe = get_end_frame(gtfline)
    if endframe != 0:
        return False

    part = get_cds_live([gtfline], cds)

    # If it's a micromicropart, can't tell at this stage if it's a stop. return true.
    if len(part) < 3:
        return True

    return is_stop_codon(part[-3:])


def seq_from_parts_list(parts_list, cds):
    """
    Translate a list of parts and return the corresponding CDS.
    """
    return translate_gtf_live([part["gtfline"] for part in parts_list], cds)


def get_part_string(prefix, parts):
    """
    Get a string representation of the part.
    """
    return prefix + "".join([f".b{part['gtfline'].start}e{part['gtfline'].end}" for part in parts])
