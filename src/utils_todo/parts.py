from utils.feature_validation import *


def is_donor_part(gtfline, cds):
    if is_terminal_part(gtfline, cds):
        return False
    if is_donor(cds[int(gtfline[4]):
        int(gtfline[4]) + 2]): return True


def is_acceptor_part(gtfline, cds):
    if is_acceptor(cds[int(gtfline[3]) - 3: int(gtfline[3]) - 1]):
        return True


def is_initial_part(gtfline, cds):
    # Must be zero-framed
    frame = int(gtfline[7])
    if frame != 0:
        return False
    # If it's a micromicropart, can't tell at this stage if it's a start. return true.
    part = get_cds_live([gtfline], cds)
    if len(part) < 3:
        return True
    return is_start_codon(part[0:3])


def is_terminal_part(gtfline, cds):
    # Must be zero-ended
    endframe = get_end_frame(gtfline)
    if endframe != 0:
        return False
    # If it's a micromicropart, can't tell at this stage if it's a stop. return true.
    part = get_cds_live([gtfline], cds)
    if len(part) < 3:
        return True
    return is_stop_codon(part[-3:])


def seq_from_parts_list(plist, cds):
    return translate_gtf_live([p["gtfline"] for p in plist], cds)


def get_part_string(prefix, parts):
    return prefix + "".join([".b" + str(part["gtfline"][3]) + "e" + str(part["gtfline"][4]) for part in parts])
