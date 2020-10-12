from Bio.Seq import Seq
from collections import Counter
import copy
import itertools
from scipy.special import binom
import re
import tempfile

from utils.blosum_matrices import *
from utils.misc import call_function, sed_file
from utils.sequence_io import read_seqs, write_seqs


def align(f_in, f_out):
    """
    Run MAFFT L-INSI on a fasta file.

    """
    call_function(f"linsi --quiet {f_in} > {f_out}")


def cds_to_aa_locs(locs):
    """
    Convert coordinates of a CDS alignment to corresponding coordinates
    of an amino acid alignment. If a CDS coordinate is not divisible by three,
    return both the lower and upper values for its corresponding AA location.

    :param locs: (list of ints) The locations
    """
    lower_locs = [a // 3 for a in locs]
    upper_locs = [(a + ((3 - a) % 3)) // 3 for a in locs]
    return sorted(list(set(lower_locs + upper_locs)))


def alignment_score(sequences, omit_empty=False, scaled=False):
    """
    Return an score for the alignment, calculated by summing column scores
    based on the Blosum matrix.

    :param sequences: (List of SeqRecords) An iterable of SeqRecord items
    :param omit_empty: (bool) Ignore any empty sequences
    :param scaled: (bool) Whether to scale the alignment score by alignment length.
    """
    if omit_empty:
        sequences = [s for s in sequences if not re.match(r'^-$', s.seq)]

    if not sequences:
        return 0

    slen = len(sequences[0])
    if not slen:
        return 0

    columns = {i: [s.seq[i] for s in sequences] for i in range(0, slen)}
    scores = {i: col_score(columns[i]) for i in columns}
    score = sum(scores.values())

    return score if not scaled else score / (1.0 * slen)


def col_score(column):
    """
    Get the alignment score for an individual column.

    :param column: A column of amino acid values.
    """
    counts = Counter(column)
    column_length = len(column)
    countsmatrix = np.zeros((len(amino_acids), column_length))

    # Else-else
    for i, j in itertools.combinations(counts.keys(), 2):
        ipos = amino_acids[i.upper()]
        jpos = amino_acids[j.upper()]
        countsmatrix[ipos, jpos] = counts[i] * counts[j]

    # Self-self
    for i in counts:
        ipos = global_blosum_positions[i.upper()]
        countsmatrix[ipos, ipos] = binom(counts[i], 2)

    # Don't count things twice.
    scoresmatrix = global_blosum_matrix * countsmatrix
    score = np.sum(scoresmatrix) / binom(column_length, 2)
    return score


def chop_alignment(alnseqs, chopper, negative=False):
    """
    Return only the elements of the alignned sequences indexed by the values in chopper

    :param alnseqs: A list of aligned sequences
    :param chopper: A list of integer indices
    :param negative: whether or not to invert the chopper list.
    """
    if not alnseqs:
        return []

    if negative:
        chopper = [a for a in range(len(alnseqs[0])) if not a in chopper]

    res = []
    for a in alnseqs:
        s = copy.deepcopy(a)
        s.seq = Seq("".join([a[i] for i in chopper if 0 <= i < len(a.seq)]))
        res.append(s)
    return res


def align_ref(fasta_new_seqs, fasta_new_seqs_out, fasta_ref_seqs, fasta_ref_seqs_out):
    """
    Add the sequences in fasta_new_seqs to the sequences in fasta_ref_seqs, performs an MSA (but
    keeping the reference alignments unchanged except for possible gaps).
    The resulting alignments for the new sequences and for the reference sequences are outputted to
    fasta_new_seqs_out and fasta_ref_seqs_out respectively.
    """
    # Perform the alignment using linsi --add.
    sed_file(fasta_ref_seqs, fasta_ref_seqs_out, r">", r">dummy.")
    call_function(f"linsi --quiet --add {fasta_new_seqs} {fasta_ref_seqs_out} > {fasta_new_seqs_out}")

    # Want to remove the ref seqs from fasta_new_seqs_out, and remove new seqs from fasta_ref_seqs_out.
    refseqs = [a for a in read_seqs(fasta_new_seqs_out) if re.match(r"^>dummy", a.id)]
    newseqs = [a for a in read_seqs(fasta_new_seqs_out) if not re.match(r"^>dummy", a.id)]

    write_seqs(newseqs, fasta_new_seqs_out)
    write_seqs(refseqs, fasta_ref_seqs_out)






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



def align_seeded(list_f_in, f_out, prealigned=False):
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
            function_string += f" --seed {f_in}"
        else:
            f_in_aln = re.sub(r"\fa", r"", f_in) + ".aln"
            align(f_in, f_in_aln)
            function_string += " --seed " + f_in_aln
    function_string += " /dev/null > " + f_out
    call_function(function_string)
    # Remove the _seed_ prefix and remove any dummy entries
    clean_dummies(f_out)


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






