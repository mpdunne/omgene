from Bio.Seq import Seq
import copy
import re
import tempfile

from utils.misc import call_function
from utils.sequence_io import read_seqs, write_seqs


def align(f_in, f_out):
    """
    Run MAFFT L-INSI on a fasta file.

    """
    call_function("linsi --quiet " + f_in + " > " + f_out)






def cds_to_aa_locs(locs):
    return list(set([a / 3 for a in locs] + [(a + ((3 - a) % 3)) / 3 for a in locs]))


def align_ref(f_in, f_out, p_ref, p_ref_out):
    call_function("sed -r \"s/>/>dummy./g\" " + p_ref + " > " + p_ref_out)
    call_function("linsi --quiet --add " + f_in + " " + p_ref_out + " > " + f_out)
    refseqs = [a for a in read_seqs(f_out) if "dummy" in a.id]
    write_seqs(refseqs, p_ref_out)
    clean_dummies(f_out)


def align_seeded(list_f_in, f_out, prealigned=False, safety=False):
    # Use each of the input fasta files as seeds. Need to double-check
    # that they each have more than one entry.
    function_string = "linsi --quiet "
    f_all = tempfile.mktemp()
    for f_in in list_f_in:
        call_function("cat " + f_in + " >> " + f_all)
        seqs = read_seqs(f_in)
        if len(seqs) == 0:
            continue
        if len(seqs) == 1:
            dummy = copy.deepcopy(seqs[0])
            dummy.id = seqs[0].id + ".dummy"
            seqs += [dummy]
        write_seqs(seqs, f_in)
        if prealigned:
            function_string += " --seed " + f_in
        else:
            f_in_aln = re.sub(r"\fa", r"", f_in) + ".aln"
            align(f_in, f_in_aln)
            function_string += " --seed " + f_in_aln
    function_string += " /dev/null > " + f_out
    call_function(function_string)
    # Remove the _seed_ prefix and remove any dummy entries
    clean_dummies(f_out)


def clean_dummies(f_out):
    toclean = read_seqs(f_out)
    cleaned = []
    dummies = []
    for c in toclean:
        if "dummy" in c.id:
            continue
        s = copy.deepcopy(c)
        s.id = re.sub(r"_seed_", r"", c.id)
        cleaned.append(s)
    write_seqs(cleaned, f_out)
    write_seqs(dummies, f_out + ".dummies")


def alignment_score(sequences, omit_empty=False, scaled=False):
    sequences = [s for s in sequences if not omit_empty or list(set(str(s.seq))) != ["-"]]
    if not sequences:
        return 0
    slen = len(sequences[0])
    if not slen:
        return 0
    dist = dict((i, [s.seq[i] for s in sequences]) for i in range(0, slen))
    scores = dict((c, col_score(dist[c])) for c in dist)
    score = sum(scores.values())
    return score if not scaled else score / (1.0 * slen)


def col_score(column):
    counts = Counter(column)
    l = len(amino_acids)
    cl = len(column)
    countsmatrix = np.zeros((l, l))
    for k in set(counts.values()):
        if not k in static_binoms:
            static_binoms[k] = bn(k, 2)
    if not cl in static_binoms:
        static_binoms[cl] = bn(cl, 2)
    # Else-else
    for i, j in itertools.combinations(counts.keys(), 2):
        ipos = static_blospos[i.upper()]
        jpos = static_blospos[j.upper()]
        countsmatrix[ipos, jpos] = counts[i] * counts[j]
    # Self-self
    for i in counts:
        ipos = static_blospos[i.upper()]
        countsmatrix[ipos, ipos] = static_binoms[counts[i]]
    # Don't count things twice.
    scoresmatrix = static_blosmat * countsmatrix
    score = np.sum(scoresmatrix) / static_binoms[cl]
    return score


def flatten_alignment(alignedseqs, gtfs, path_out="", pre_sorted=False):
    # Flatten_alignment accepts a list of aligned Seqs, and a dict of loose leaf gtfs.
    # It returns a list of flattened Seqs, and a set of aa ands cds coordinates for the gtfs, and the sorted set of gtfs
    aa = ["A", "G"]
    flatseqs = []
    good_gtfs = {}
    locations_aa = {}
    locations_cds = {}
    for s in alignedseqs:
        gtf = gtfs[s.id] if pre_sorted else safe_gtf(gtfs[s.id])
        good_gtfs[s.id] = gtf
        if not gtf:
            continue
        # The stop codon won't be included in the alignment.
        seq_expanded = "".join([a * 3 for a in str(s.seq)])
        seq_expanded += "@@@"
        seq_expanded = re.sub(r"(-*)@@@$", r"***\1", seq_expanded)
        # Furthermore, partial codons are possible, so tag some crap on the end for safety
        seq_expanded += "###"
        res = {}
        pos = 0
        flat = ""
        strand = gtf[0][6]
        for i, line in enumerate(gtf):
            cdslen = int(line[4]) - int(line[3]) + 1
            localpos = 0
            localres = []
            while localpos < cdslen:
                if len(seq_expanded) <= pos or seq_expanded[pos] != "-":
                    localres.append(pos)
                    flat += aa[i % 2]
                    localpos += 1
                else:
                    flat += "-"
                pos += 1
            res[i] = sorted(localres)
        outstr = ""
        for qpos in range(len(s)):
            pos = 3 * qpos
            triplet = flat[pos:pos + 3]
            ender = "A" if triplet == "AAA" else ("G" if triplet == "GGG" else ("-" if triplet == "---" else "Q"))
            outstr += ender
        outstr = re.sub(r"Q([Q-]*$)", r"", outstr)
        outstr = outstr + (len(s) - len(outstr)) * "-"
        outstr = re.sub(r"(-*)E", r"E\1", outstr + "E")
        outstr = re.sub(r"^(-*)[AG]", r"\1S", outstr)
        t = copy.deepcopy(s)
        t.seq = Seq(outstr)
        flatseqs.append(t)
        locations_cds[s.id] = res
        locations_aa[s.id] = dict((i, [a / 3 for a in res[i]]) for i in res)
    if path_out:
        write_seqs(flatseqs, path_out)
    return locations_aa, locations_cds, good_gtfs





def chop_alignment(alnseqs, chopper, negative=False):
    """chop an alignment based on a list of coordinates
    """
    if not alnseqs:
        return []
    if negative:
        chopper = [a for a in range(len(alnseqs[0])) if not a in chopper]
    return list(chop(alnseqs, chopper))


def chop(alnseqs, chopper):
    res = []
    for a in alnseqs:
        s = copy.deepcopy(a)
        s.seq = Seq("".join([a[i] for i in chopper if 0 <= i < len(a.seq)]))
        res.append(s)
    return res
