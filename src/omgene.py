#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2017 Michael Dunne
#
# OMGene version 0.1.0
#
# This program (OMGene) is distributed under the terms of the GNU General Public License v3
#
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# For any enquiries send an email to Michael Dunne
# michael.dunne@worc.ox.ac.uk

__author__ = "Michael Dunne"
__credits__ = "Michael Dunne, Steve Kelly"

import csv
import re
import os
import sys
import itertools
import copy
import subprocess
import multiprocessing
import commands
import datetime
import tempfile
import pickle
import string
import pyfaidx
import argparse
import math
import random

import networkx as nx
from shutil import copyfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.special import binom as bn

from utils.blosum_matrices import *

static_blosdict = blosum_dict(-1, 0)
static_blosmat = blosum_matrix(static_blosdict)
static_blospos = blosum_positions()

# Save time later by storing calculated binomial values
static_binoms = {}


##########################################
# Check programs
##########################################

def can_run_command(command, q_allow_stderr=False):
    capture = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = [x for x in capture.stdout]
    stderr = [x for x in capture.stderr]
    return check(len(stdout) > 0 and (q_allow_stderr or len(stderr) == 0))


def can_run_specific(line, line_formatted):
    return True if can_run_command(line) else output_bool(False, "ERROR: Cannot run " + line_formatted,
                                                          "    Please check " + line_formatted + " is installed and that the executables are in the system path\n")


def can_run_minus_h(package, package_formatted):
    return can_run_specific(package + " -h", package_formatted)


def can_run_awk():
    return returns_expected("awk '{print 4 + 5}'", "awk '{print 4 + 5}'", "a", "9\n")


def can_run_man(package, package_formatted):
    return can_run_specific("man " + package, package_formatted)


def check(bool_val, true_msg=" - ok", false_msg=" - failed", true_extra="", false_extra=""):
    return output_bool(True, true_msg, true_extra) if bool_val else output_bool(False, false_msg, false_extra)


def returns_expected(message, command, contents_in, expected_out):
    # sys.stdout.write("Test can run \""+message+"\"\n")
    path_tf1 = tempfile.mktemp()
    path_tf2 = tempfile.mktemp()
    write(contents_in, path_tf1)
    call_function(command + " " + path_tf1 + " > " + path_tf2)
    res = read(path_tf2)
    return res == expected_out


def output_bool(bool_val, msg, msg_extra):
    if msg_extra:
        print(msg_extra)
    return bool_val


def check_shell():
    """ Run quick checks to ensure the relevant packages are installed, and that they do
        the things we need them to do (applies in particular to bedtools and to R).
        - awk           - exonerate   - bedtools
        - echo          - sed
        - cut           - grep
        - sort          -tac/cat
    """
    sprint("Checking installed programs...")
    checks = []
    for program in ["sed", "echo", "cut", "sort", "grep", "uniq", "tac", "cat", "mktemp"]:
        checks.append(can_run_man(program, program))
    for program in ["exonerate", "bedtools"]:
        checks.append(can_run_minus_h(program, program))
    checks += [can_run_awk()]
    if not all(checks):
        print(
            "\n_some programs required to run omgene are not installed or are not callable.\n_please ensure all of the above programs are installed and in the system path.")
        sys.exit()


########################################
# Generic Utilities
########################################

def lprint(l):
    """Print a list piece by piece
    """
    for i in l:
        print(i)


def dprint(d):
    """Print a dict piece by piece
    """
    for k in d:
        print(k, d[k])


def print_seqs(seqs):
    for s in seqs:    sprint(str(s.seq))


def sprint(string):
    print(string)


def read(path_file, tag="r"):
    with open(path_file, tag) as f:
        return f.read()


def read_csv(path_file, ignore_blank=True, ignore_hash=False):
    with open(path_file, "r") as f:
        data = [line for line in csv.reader(f, delimiter="\t") if \
                (not ignore_blank or ''.join(line).strip()) and (not ignore_hash or not line[0].startswith('#'))]
    return data


def write_csv(rows, path_file):
    with open(path_file, "w") as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, quotechar='')
        writer.writerows(rows)


def write(msg, path_dest):
    with open(path_dest, "w") as f:
        f.write(str(msg))


def make_if_absent(path_dir):
    if not os.path.exists(path_dir): 
        os.makedirs(path_dir)
    return path_dir


def concat_files(fileslist, dest):
    with open(dest, "w") as o:
        for file in fileslist:
            with(open(file, "r")) as b: o.write(b.read())


def call_function(str_function):
    """Call a function in the shell
    """
    subprocess.call([str_function], shell=True)


def grab_lines(command):
    return commands.getstatusoutput(command)[1].split("\n")


def interpolate(l1, l2=[]):
    return range(min(l1 + l2), max(l1 + l2))


def delete_if_present(path):
    """Delete a folder which may or may not exist
    """
    try:
        os.remove(path)
    except OSError:
        pass


def pasync(pool, function, args):
    """Run asynchronously
    """
    pool.apply_async(function, args=args)


def pause():
    call_function("sleep 10000")


def check_file_exists(path_file):
    return path_file if os.path.exists(path_file) else sys.exit("File does not exist: " + path_file)


def idict(keys, defaultvalue):
    return dict((k, copy.deepcopy(defaultvalue)) for k in keys)


def del_indices(l, indices):
    # Only delete indices that exist, and do them in the right order!
    k = l[:]
    indices_sorted = sorted(set(i % len(l) for i in indices if -len(l) <= i < len(l)), key=lambda x: -x)
    for i in indices_sorted: del k[i]
    return k


def flatten_lol(lol):
    """Flatten a list of lists
    """
    return [y for x in lol for y in x]


def de_dup(dlist):
    nlist = []
    for i in dlist:
        if i not in nlist:
            nlist.append(i)
    return nlist


def merge_dicts(dicts):
    # Note - identical keys are flattened.
    res = {}
    for d in dicts:
        for k in d:
            res[k] = d[k]
    return res


def wrap_string(string, n):
    return [string[start:start + n] for start in range(0, len(string), n)]


def blank_file(path_file):
    dirname = os.path.dirname(path_file)
    make_if_absent(dirname)
    call_function("echo -n \"\" >" + path_file)
    return path_file


def anyperm(listin, checker):
    return any([(list(i) in checker) for i in itertools.permutations(listin)])


########################################
# Feature validation
########################################

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


########################################
# Binary algebra
########################################

from utils.binary_algebra import *


#########################################
# Sequence I/O
#########################################

def read_seqs(path_in, tolist=True):
    seqs = SeqIO.parse(path_in, "fasta")
    if tolist:
        return list(seqs)
    else:
        return seqs


def write_seqs(seqs, path_out):
    SeqIO.write(seqs, path_out, "fasta")


def make_seq(label, aa):
    return SeqRecord(id=label, name="", description="", seq=Seq(aa))


#########################################
# Live gtf operations
#########################################

def abs_frame(gtfline):
    return (int(gtfline[3]) + int(gtfline[7])) % 3


def gtf_compatible_frame(gtfline1, gtfline2):
    return abs_frame(gtfline1) == abs_frame(gtfline2)


def get_protoexon(gtf):
    return [int(gtf[3]), int(gtf[4]), int(gtf[7])]


def gtf_length(line):
    return int(line[4]) - int(line[3]) + 1


def gtf_lines_overlap(i, j):
    return (int(i[3]) <= int(j[3]) and int(j[3]) <= int(i[4])) or (int(j[3]) <= int(i[3]) and int(i[3]) <= int(j[4]))


def same_frame(i, j):
    return (int(i[3]) + int(i[7])) % 3 == (int(j[3]) + int(j[7])) % 3


def overlap_in_frame(i, j):
    return gtf_lines_overlap(i, j) and same_frame(i, j)


def gtf_line_equals(gtf1, gtf2):
    return int(gtf1[3]) == int(gtf2[3]) and int(gtf1[4]) == int(gtf2[4]) and gtf1[0] == gtf2[0]


def safe_gtf(gtflines):
    return sorted(de_dup([a for a in gtflines if a[2].lower() == "cds"]),
                  key=lambda x: (math.pow(-1, x[6] == "-")) * int(x[3]))


def clean_gtf_live(gtf):
    return [line for line in gtf if not is_bad_gtf([line])]


def is_bad_gtf(gtf):
    return any(int(line[3]) >= int(line[4]) for line in gtf)


def get_cds_live(gtf, cds):
    return "".join([cds[int(line[3]) - 1: int(line[4])] for line in sorted(gtf, key=lambda x: int(x[3]))])


def translate_gtf_live(gtf, cds):
    return translate_part(get_cds_live(gtf, cds))


def translate_part(cds_part):
    return str(Seq(cds_part[0:len(cds_part) - (len(cds_part) % 3)]).translate())


def get_end_frame(gtf):
    return (int(gtf[4]) - int(gtf[3]) - int(gtf[7]) + 1) % 3


def get_overlaps(list_gtf):
    overlaps = []
    safe = list(list_gtf)
    for i in list_gtf:
        for j in list_gtf:
            if i == j: 
                continue
            if gtf_lines_overlap(i, j):
                if not i in overlaps: 
                    overlaps.append(i)
                if i in safe: 
                    safe.remove(i)
    return safe, overlaps


def get_non_overlapping_subsets(overlaps):
    # For a set of possibly overlapping gtflines, returns
    # all maximally non-overlapping subsets
    d_o = {}
    for i, e in enumerate(overlaps): d_o[i] = e
    G = nx.Graph()
    keys = d_o.keys()
    for i in keys: G.add_node(i)
    for i in itertools.product(keys, keys):
        if i[0] != i[1] and gtf_lines_overlap(d_o[i[0]], d_o[i[1]]):
            G.add_edge(i[0], i[1])
    H = nx.complement(G)
    C = nx.find_cliques(H)
    return [[d_o[i] for i in k] for k in C]


########################################
# Part group operations
########################################

def group_parts(dict_parts, ranges=True):
    part_graph = nx.Graph()
    for generegion in dict_parts:
        for key in dict_parts[generegion]:
            tag = generegion + "." + str(key)
            if not tag in part_graph.nodes(): 
                part_graph.add_node(tag)
            part = dict_parts[generegion][key]
            for ogeneregion in dict_parts:
                if generegion == ogeneregion: 
                    continue
                for okey in dict_parts[ogeneregion]:
                    otag = ogeneregion + "." + str(okey)
                    opart = dict_parts[ogeneregion][okey]
                    f_o, b_o = overlap_proportion(part, opart, ranges)
                    if max(f_o, b_o) > 0.333:
                        part_graph.add_edge(tag, otag)
    C = list(nx.find_cliques(part_graph))
    groups = [sorted(c) for c in sorted(C, key=lambda x: sorted(x)[0])]
    a = sort_parts_list(C)
    return C, a


def sort_parts_list(list_c):
    if len(list_c) == 1: 
        return list_c
    G = nx.Di_graph()
    for i, c in enumerate(list_c):
        for j, d in enumerate(list_c):
            if lex_less(c, d): 
                G.add_edge(i, j)
    return [list_c[i] for i in nx.topological_sort(G)]


def lex_less(item1, item2):
    if item1 == item2: 
        return False
    d1 = {}
    d2 = {}
    for i in item1:
        m1, o1 = parse_id(i)
        d1[m1] = o1
    for j in item2:
        m2, o2 = parse_id(j)
        d2[m2] = o2
    isec = set(d1.keys()).intersection(set(d2.keys()))
    if all([d1[i] == d2[i] for i in isec]): 
        return False
    return all([int(d1[i]) <= int(d2[i]) for i in isec])


def overlap_proportion(a1, a2, ranges=True):
    a1r = range(min(a1), max(a1) + 1) if ranges else a1
    a2r = range(min(a2), max(a2) + 1) if ranges else a2
    if len(a1r) * len(a2r) == 0: 
        return 0, 0
    isec = 1.0 * len(set(a1r).intersection(a2r))
    return isec / len(a1r), isec / len(a2r)


########################################
# Gene part operations
########################################

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


########################################
# Id parsing
########################################

def parse_id(str_id):
    generegion = re.sub(r"([^.]*)\.(.*)", r"\1", str_id)
    option = int(re.sub(r"([^.]*)\.([0-9]*)", r"\2", str_id))
    return generegion, option


def get_gr(str_in):
    return re.sub(r".*(generegion_[0-9]+).*", r"\1", str_in)


########################################
# Alignments
########################################

def align(f_in, f_out, double=False, safety=False):
    align_vanilla(f_in, f_out)


def align_vanilla(f_in, f_out):
    call_function("linsi --quiet " + f_in + " > " + f_out)


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


def cds_to_aa_locs(locs):
    return list(set([a / 3 for a in locs] + [(a + ((3 - a) % 3)) / 3 for a in locs]))


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


########################################
# Feature finding
########################################

def find_starts(proto_exon, cds, left=True, right=True, central=True, happy_with_central=True):
    """Find starts in either direction.
    """
    pos, boundary, frame = proto_exon
    c = cds[pos - 1:pos + 2]
    cds = cds.lower()
    out = []
    if central and frame == 0 and is_start_codon(c):
        out = [pos]
        if happy_with_central: 
            return out
    if right: 
        out += walk_termini(cds, pos + frame, 3, lambda c: is_stop_codon(c), lambda x: x + 1 > boundary,
                                  lambda c: is_start_codon(c), lambda x, c: c[x - 1:x + 2])
    if left: 
        out += walk_termini(cds, pos + frame, -3, lambda c: is_stop_codon(c), lambda x: x - 1 < 0,
                                 lambda c: is_start_codon(c), lambda x, c: c[x - 1:x + 2])
    return out


def find_stops(proto_exon, cds):
    """ Find the first stop codon in the right hand direction.
    """
    left, right, frame = proto_exon
    pos = right + ((3 - (1 + right - left - frame)) % 3) - 3
    cds = cds.lower()
    c = cds[pos - 3:pos]
    if is_stop_codon(c):
        return [int(pos)]
    else:
        return walk_termini(cds, pos, 3, lambda c: False, lambda x: x + 1 > len(cds), lambda c: is_stop_codon(c),
                            lambda x, c: c[x - 3:x])


def walk_termini(cds, pos, step, esc_cod, esc_pos, win_fn, cds_slice):
    while True:
        pos = pos + step
        c = cds_slice(pos, cds)
        if esc_cod(c) or esc_pos(pos): 
            break
        if win_fn(c):
            return [int(pos)]
    return []


def contains_stop(sequence, frame=0):
    return any(is_stop_codon(sequence[i:i + 3]) for i in range(frame, len(sequence), 3))


def walk_splice(pos, cds, step, esc_pos, esc_cod, nframe, items, direction, cds_slice, cod_test):
    frames_done = []
    while True:
        nframe = (nframe + direction * step) % 3
        pos = pos + step
        if esc_pos(pos): 
            break
        c = cds_slice(pos, cds)
        if cod_test(c):
            if esc_cod(pos, nframe): 
                break
            if not nframe in frames_done: 
                items[nframe].append(int(pos))
            frames_done += [nframe]
            if len(set(frames_done)) == 3: 
                break
    return items


def walk_splice_left(pos, cds, step, esc_pos, esc_cod, nframe, lefts):
    return walk_splice(pos, cds, step, esc_pos, esc_cod, nframe, lefts, -1, lambda x, c: c[x - 3:x - 1],
                       lambda x: is_acceptor(x))


def walk_splice_right(pos, cds, step, esc_pos, esc_cod, nframe, rights):
    return walk_splice(pos, cds, step, esc_pos, esc_cod, nframe, rights, 1, lambda x, c: c[x:x + 2],
                       lambda x: is_donor(x))


def find_splice_sites(proto_exon, cds, left=True, right=True, happy_with_central=True, first=True):
    # LHS splice sites need to lend the right frame.
    frame = proto_exon[2]
    cds = cds.lower()
    lefts = idict([0, 1, 2], [])
    if left:
        pos, boundary = proto_exon[0:2]
        c = cds[pos - 3:pos - 1]
        good = is_acceptor(c)
        if good: 
            lefts[frame].append(int(pos))
        if (not good) or (good and not happy_with_central):
            lefts = walk_splice_left(pos, cds, -1, lambda x: x - 3 < 0,
                                     lambda x, f: contains_stop(cds[x - ((3 - f) % 3) - 1:x + frame], 0), frame, lefts)
            lefts = walk_splice_left(pos, cds, 1, lambda x: x >= boundary, lambda x, f: False, frame, lefts)
    # The frames here are donor frames.
    rights = idict([0, 1, 2], [])
    if right:
        boundary, pos, oframe = proto_exon
        c = cds[pos:pos + 2]
        good = is_donor(c)
        frame = (pos - boundary + 1 - proto_exon[2]) % 3
        if good: 
            rights[frame].append(int(pos))
        if (not good) or (good and not happy_with_central):
            rights = walk_splice_right(pos, cds, -1, lambda x: x < boundary, lambda x, f: False, frame, rights)
            rights = walk_splice_right(pos, cds, 1, lambda x: x + 1 > len(cds) - 1,
                                       lambda x, f: contains_stop(cds[boundary - ((3 - oframe) % 3) - 1:x], 0), frame,
                                       rights)
    return lefts, rights


def split_at_stops(ko, string_context, min_exon=20):
    happy = []
    waiting = ko
    original = ko[:]
    while waiting:
        waiting_next = []
        for line in waiting:
            cds = string_context[int(line[3]) - 1:int(line[4])]
            frame = int(line[7])
            scanstart = frame
            scanend = frame + 3
            safe = True
            while scanend <= len(cds):
                codon = cds[scanstart:scanend]
                if is_stop_codon(codon):
                    newline1 = line[:]
                    newline2 = line[:]
                    newline1[4] = int(line[3]) + scanstart - 1
                    newline2[3] = int(line[3]) + scanend
                    newline2[7] = 0
                    waiting_next += [newline1, newline2]
                    safe = False
                    break
                scanstart += 3
                scanend += 3
            if safe: 
                happy.append(line)
        waiting = waiting_next
    # Throw out any bad exons. Also throw out any of the new split bits
    # if they are tiny: these bits are probably rubbish.
    happy = [i for i in happy if int(i[4]) >= int(i[3])]
    happy = [i for i in happy if [a for a in original if gtf_line_equals(a, i)] or int(i[4]) - int(i[3]) >= min_exon]
    return happy


##########################################
# File-based gtf/fasta operations
##########################################

def translate_gtf_to_file(gtf, cds, tag, path_out):
    aa_parts = translate_gtf(gtf, cds, tag)
    if aa_parts:
        write_seqs(aa_parts, path_out)
    else:
        blank_file(path_out)


def translate_gtf(gtf, cds, tag):
    for line in gtf:
        cdspart = cds[int(line[3]) + int(line[7]) - 1: int(line[4])]
        aapart = translate_part(cdspart)
        if aapart: 
            yield make_seq(tag + ".b" + str(line[3]) + "e" + str(line[4]), aapart)


def framify(path_gtf):
    """Framify assumes that the first exon is in frame 0
    """
    data = read_csv(path_gtf)
    gene_names = list(set([a[1] for a in data]))
    genes = [[a for a in data if a[1] == gene] for gene in gene_names]
    results = []
    for g in genes:
        gs = sorted([e for e in g if e[2].lower() == "cds"], key=lambda x: int(x[3]))
        lastframe = 0
        for cds in gs:
            cdsc = cds[:]
            cdsc[7] = lastframe
            results.append(cdsc)
            nextframe = (3 - (int(cds[4]) - int(cds[3]) + 1 - lastframe)) % 3
            lastframe = nextframe
    write_csv(results, path_gtf)


def split_to_abs_frame(gtf):
    gtf_framed = idict([0, 1, 2], [])
    for line in gtf:
        start = int(line[3])
        frame_rel = int(line[7])
        frame_abs = (start + frame_rel) % 3
        gtf_framed[frame_abs].append(line)
    return gtf_framed


def merge_stranded(path_in, path_out, str_infmt="gtf"):
    # GTF format will have strand in column 7, bed format is in col 6.
    col = "6" if str_infmt == "bed" else "7"
    # Split by strand
    call_function("awk '$" + col + "==\"+\"' " + path_in + " > " + path_in + ".pos")
    call_function("awk '$" + col + "==\"-\"' " + path_in + " > " + path_in + ".neg")
    # Merge separately
    call_function(
        "sort -k4,4n " + path_in + ".pos | bedtools merge -i - | cut -f1-3 | sed -r \"s/$/\\t+/g\" > " + path_out)
    call_function(
        "sort -k4,4n " + path_in + ".neg | bedtools merge -i - | cut -f1-3 | sed -r \"s/$/\\t-/g\" >> " + path_out)
    # Kill the files
    delete_if_present(path_in + ".pos")
    delete_if_present(path_in + ".neg")


def merge_friendly(path_gtf, path_out, sort=True):
    """Merge elements of a gtf file, making sure only to merge
       which are frame-compatible. Requires frame info.
    """
    gtf = read_csv(path_gtf)
    gtf_framed = split_to_abs_frame(gtf)
    out = []
    for i in gtf_framed:
        if not any(gtf_framed[i]): 
            continue
        path_gtftmp = path_gtf + ".f" + str(i) + ".tmp"
        write_csv(gtf_framed[i], path_gtftmp)
        merge_stranded(path_gtftmp, path_gtftmp + ".merged", "gtf")
        merged = read_csv(path_gtftmp + ".merged")
        model = gtf_framed[i][0]
        for line in merged:
            newline = model[:]
            newline[1] = "."
            newline[8] = "source_id \"any\""
            newline[3] = int(line[1]) + 1
            newline[4] = line[2]
            newline[7] = (i - int(newline[3])) % 3
            out.append(newline)
    if sort: 
        out = sorted(out, key=lambda x: int(x[3]))
    write_csv(out, path_out)


def clean_gtf(path_gtf_o, path_gtf):
    call_function("sort -u " + path_gtf_o + " | sort -k4,4n > " + path_gtf)


def replace_seq(seq_object, new_seq):
    s = copy.deepcopy(seq_object)
    s.seq = new_seq
    return s


def fetch_aa(path_gtf_in, path_genome, path_cds_fasta_out, path_aa_fasta_out, int_translation_table):
    """Fetch cds and translate to aa
    """
    fetch_cds(path_gtf_in, path_genome, path_cds_fasta_out, "CDS")
    prot_sequences = [replace_seq(s, s.seq.translate(table=int_translation_table)) for s in
                      read_seqs(path_cds_fasta_out, False)]
    write_seqs(prot_sequences, path_aa_fasta_out)


def fetch_cds(path_gtf_in, path_genome, path_cds_fasta_out, str_token):
    """Fetch the cds fasta sequences for the given gtfs.
    """
    call_function(
        "infile=" + path_gtf_in + "; outfile=" + path_cds_fasta_out + "; genome=" + path_genome + "; token=" + str_token + """
        tf=`mktemp -d`
        gtf_cds="$tf/gtf_cds"
        gtf_bed="$tf/gtf_bed"

        #Prepare the gtf
        grep -v_p "^$" $infile | awk -v token="$token" 'toupper($3)==toupper(token)' > $gtf_cds
        cut -f1-8 $gtf_cds > $gtf_bed.1
        sed -r "s/.*transcript_id[ =]\\"?([^\\";]*)\\"?;?.*/\\1/g" $gtf_cds > $gtf_bed.2
        paste $gtf_bed.1 $gtf_bed.2 | perl -ne 'chomp; @l=split; printf "$l[0]\\t%s\\t$l[4]\\t$l[8]\\t.\\t$l[6]\\n", $l[3]-1' | sort -u | sort -k1,1V -k2,2n > $gtf_bed

        #Negative strand
        awk '$6=="-"' $gtf_bed > $gtf_bed.neg
        bedtools getfasta -name -s -full_header -fi $genome -fo $gtf_bed.neg.tab -bed $gtf_bed.neg -tab
        tac $gtf_bed.neg.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtf_bed.neg.fa

        #Then positive strand
        awk '$6=="+"' $gtf_bed > $gtf_bed.pos
        #cat $gtf_bed.pos
        bedtools getfasta -name -s -full_header -fi $genome -fo $gtf_bed.pos.tab -bed $gtf_bed.pos -tab
        cat $gtf_bed.pos.tab | awk '{a[$1]=a[$1]""$2} END {for (i in a) {print ">"i"\\n"a[i]}}' > $gtf_bed.pos.fa

        cat $gtf_bed.pos.fa $gtf_bed.neg.fa | sed -r "s/^>(.*)$/£££>\\1###/g" | sed -r \"s/$/###/g\" | tr '\\n' ' ' | sed -r "s/£££/\\n/g" | sed -r "s/### //g" | grep -v XXX | grep -v "\*[A-Z]" | grep -v "###$" | sed -r "s/###/\\n/g" | grep -v_p "^$" > $outfile

        rm -r $tf
        """)


def get_base(path_gtf, int_slop, dict_chr_sizes, path_base, path_base_rel, path_base_tight, sequence_id):
    """Get the base for a gtf file, with margins of a specified size
    """
    gtf = read_csv(path_gtf)
    # Use the data from the first line, but using min and max location coordinates.
    gene_end = max([int(b[4]) for b in gtf])
    gene_start = min([int(b[3]) for b in gtf])
    # Grab the gene and transcript names where possible
    str_gene_id = re.sub(r'.*gene_id \"([^\"]*)\".*', r'\1', [a[8] for a in gtf if "gene_id" in a[8]][0])
    str_transcript_id = re.sub(r'.*transcript_id \"([^\"]*)\".*', r'\1',
                               [a[8] for a in gtf if "transcript_id" in a[8]][0])
    # Slop the entry, ensuring that the slopped coordinates are still viable within the genome
    entry = gtf[0][0:9]
    chr_name = entry[0]
    chr_size = dict_chr_sizes[chr_name]
    descr = "transcript_id \"" + str_transcript_id + ".t1\"; gene_id \"" + str_gene_id + "\""
    entry = make_gtf_line(entry, max(int(gene_start) - int_slop, 1), min(int(gene_end) + int_slop, chr_size), chr_name,
                          "", descr, "base")
    # Construct both slopped and unslopped versions of the base.
    entry.append(sequence_id)
    entry_tight = make_gtf_line(entry, gene_start, gene_end, "", "", "", "base_tight")
    write_csv([entry], path_base)
    write_csv([entry_tight], path_base_tight)
    # Construct the relative gtf for the slopped base
    entry_rel = make_gtf_line(entry, 1, entry[4] - entry[3] + 1, "relative", "", "", "")
    write_csv([entry_rel], path_base_rel)


def make_gtf_line(refline, pstart, pend, chrom, strand, descr, tag):
    line = refline[:]
    if chrom: 
         line[0] = chrom
    if tag: 
           line[2] = tag
    if pstart: 
        line[3] = pstart
    if pend: 
          line[4] = pend
    if strand: 
        line[6] = strand
    if descr: 
         line[8] = descr
    return line


def toggle_relative(path_gtf_in, path_gtf_out, path_gtf_base, cds_only=False, tag="", dirn="fwd"):
    # Read in the data
    data_gtf_in = read_csv(path_gtf_in)
    data_gtf_base = read_csv(path_gtf_base)[0]
    if cds_only: 
        data_gtf_in = [a for a in data_gtf_in if a[2] == "CDS"]
    strand = data_gtf_base[6]
    g_start_o = int(data_gtf_base[3])
    g_end_o = int(data_gtf_base[4])
    chrom = data_gtf_base[0]
    descr = data_gtf_base[8]
    # Go through each line and flip it.
    data_gtf_flipped = []
    for oline in data_gtf_in:
        if dirn == "fwd":
            l_end = g_end_o - int(oline[4]) + 1 if strand == "-" else int(oline[3]) - g_start_o + 1
            r_end = g_end_o - int(oline[3]) + 1 if strand == "-" else int(oline[4]) - g_start_o + 1
        else:
            l_end = g_end_o - int(oline[4]) + 1 if strand == "-" else g_start_o + int(oline[3]) - 1
            r_end = g_end_o - int(oline[3]) + 1 if strand == "-" else g_start_o + int(oline[4]) - 1
        line = make_gtf_line(oline, l_end, r_end, chrom, "+" if dirn == "fwd" else strand, descr, "")
        if tag != "": 
            line[1] = tag
        data_gtf_flipped.append(line)
    write_csv(data_gtf_flipped, path_gtf_out)


########################################
# Setup functions
########################################

def get_genome_info_indi(path_genome):
    fas = pyfaidx.Fasta(path_genome)
    return dict((record.name, len(record)) for record in fas)


def get_genome_info(dict_seq_info):
    """Get the sizes of all chromosomes in each genome
    """
    genome_files = list(set([dict_seq_info[b]["genome"] for b in dict_seq_info]))
    genome_sizes = dict((path_genome, get_genome_info_indi(path_genome)) for path_genome in genome_files)
    return genome_sizes


def read_input_locations(path_inf, path_ref):
    """Read CSV file containing the locations for the sequence files, etc.
       Each row is a sequence. The columns are [gtf, genome].
    """
    sprint("loading and checking input data locations...")
    data = read_csv(path_inf, ignore_hash=True)
    dict_seq_info = read_input_file(read_csv(path_inf, ignore_hash=True), {}, False)
    if path_ref: 
        dict_seq_info = read_input_file(read_csv(path_ref, ignore_hash=True), dict_seq_info, True)
    return dict_seq_info


def process_codon_option(co, length):
    s_co = co.split(",")
    for a in s_co:
        a = re.sub(r" ", r"", a)
        if not len(a) == length or any(not s.lower() in "acgtu" for s in a):
            sys.exit("Invalid start codon/splice site choice: " + str(a))
        yield re.sub(r"u", r"t", a.lower())


def process_codon_options(do, ac, sc):
    s_do = list(set(process_codon_option(do, 2)))
    s_ac = list(set(process_codon_option(ac, 2)))
    s_sc = list(set(process_codon_option(sc, 3)))
    return s_do, s_ac, s_sc


def read_input_file(data, dict_seq_info, keep=False):
    next_index = len(dict_seq_info)
    for i, line in enumerate(data):
        # Check files exist
        path_gtf = check_file_exists(line[0])
        path_genome = check_file_exists(line[1])
        # TODO CUPCAKES:
        # need to check chromosomes, coordinates, and that each gtf refers to a single gene
        # Also need to check that each genome id is unique.
        seq_id = "sequence_" + str(next_index + i)
        dict_seq_info[seq_id] = {"id": seq_id, "gtf": path_gtf, "genome": path_genome,
                                 "gtf_id": os.path.basename(path_gtf), "genome_id": os.path.basename(path_genome),
                                 "keep": keep}
    return dict_seq_info


def prepare_output_folder(path_out_dir):
    path_out_dir = make_if_absent(path_out_dir if path_out_dir else generate_dir())
    return make_if_absent(path_out_dir + "/results"), make_if_absent(path_out_dir + "/working")


def generate_dir():
    return tempfile.mkdtemp(prefix="OMGene_out_" + datetime.datetime.now().strftime("%y%m%d") + "_Run_id_", dir="..")


def prepare_sequences(path_w_dir, dict_seq_info, dict_genome_info, int_num_cores, int_slop_amount):
    pool = multiprocessing.Pool(int_num_cores)
    path_s_dir = make_if_absent(path_w_dir + "/sequences")
    for seq_id, ds in dict_seq_info.items():
        # Get the base gtf and its CDS. Organise them nicely into folders.
        path_s_dir_seqid = make_if_absent(path_s_dir + "/" + seq_id)
        ds["gtf_base"] = path_base = path_s_dir_seqid + "/" + seq_id + ".base"
        ds["gtf_base_tight"] = path_base_tight = path_base + ".tight"
        ds["gtf_base_rel"] = path_base_rel = path_base + ".relative"
        ds["cds_base"] = path_cds_fasta_out_base = path_base + ".cds.fasta"
        # Get gtf files
        path_gtf_o = ds["gtf"]
        path_genome = ds["genome"]
        dict_chr_sizes = dict_genome_info[path_genome]
        # Clean up the gtf file and move it to the working directory.
        ds["gtf"] = path_gtf = path_s_dir_seqid + "/" + seq_id + ".gtf"
        clean_gtf(path_gtf_o, path_gtf)
        pasync(pool, grab_bases, args=(
            path_gtf, int_slop_amount, dict_chr_sizes, path_base, path_base_rel, path_base_tight, path_genome,
            path_cds_fasta_out_base, seq_id))
        # Get the cds and aa for the genes
        ds["aa"] = path_aa_fasta_out = path_s_dir_seqid + "/" + seq_id + ".aa.fasta"
        ds["cds"] = path_cds_fasta_out = path_s_dir_seqid + "/" + seq_id + ".cds.fasta"
        pasync(pool, fetch_aa, args=(path_gtf, path_genome, path_cds_fasta_out, path_aa_fasta_out, 1))
    pool.close()
    pool.join()
    for ds in dict_seq_info.values():
        ds["sequence_aa"] = str(read_seqs(ds["aa"])[0].seq).strip("-")
        ds["sequence_cds"] = str(read_seqs(ds["cds"])[0].seq).strip("-")


def grab_bases(path_gtf, int_slop_amount, dict_chr_sizes, path_base, path_base_rel, path_base_tight, path_genome,
               path_cds_fasta_out_base, sequence_id):
    get_base(path_gtf, int_slop_amount, dict_chr_sizes, path_base, path_base_rel, path_base_tight, sequence_id)
    fetch_cds(path_base, path_genome, path_cds_fasta_out_base, "base")


def write_gene_region_file(line, generegion, path_generegion):
    nline = [line[0], "orthofxer", "generegion", line[1], line[2],
             ".", line[3], ".", "transcript_id \"" + generegion + "\"; gene_id \"" + generegion + "\"", generegion]
    write_csv([nline], path_generegion)
    return path_generegion


def sequence_matches_region(seq_gtf, region_gtf):
    return [l for l in grab_lines("bedtools intersect -a " + region_gtf + " -b " + seq_gtf) if ''.join(l).strip()]


def prepare_generegions(dict_seq_info, dict_genome_info, path_w_dir, int_num_cores, int_slop_amount):
    dict_generegions = {}
    path_m_dir = make_if_absent(path_w_dir + "/generegions")
    path_g_dir = make_if_absent(path_m_dir + "/genomes")
    pool = multiprocessing.Pool(int_num_cores)
    for path_genome in list(set([a["genome"] for a in dict_seq_info.values()])):
        # Join all tight base files together and then perform a bedtools merge
        genome_id = os.path.basename(path_genome)
        bases_tight = [a["gtf_base_tight"] for a in dict_seq_info.values() if a["genome"] == path_genome]
        path_all_base_loci = path_g_dir + "/" + genome_id + ".bases"
        path_merged_base_loci = path_g_dir + "/" + genome_id + ".bases.merged"
        concat_files(bases_tight, path_all_base_loci)
        merge_stranded(path_all_base_loci, path_merged_base_loci, "gtf")
        sequences = [copy.deepcopy(a) for a in dict_seq_info.values() if a["genome"] == path_genome]
        for line in read_csv(path_merged_base_loci):
            # Initialise the generegion in the holder dict
            generegion = "generegion_" + str(len(dict_generegions.keys()))
            dg = dict_generegions[generegion] = {}
            dg["m_dirl"] = path_m_dir_l = make_if_absent(path_m_dir + "/" + str(generegion))
            dg["genome"] = path_genome
            dg["generegion"] = generegion
            # The path of the generegion base descriptor
            path_generegion = write_gene_region_file(line, generegion, path_m_dir_l + "/" + str(generegion) + ".gtf")
            # Figure out which sequences correspond to this generegion.
            matchingseqs = []
            for s in sequences:
                if sequence_matches_region(s["gtf"], path_generegion):
                    dict_seq_info[s["id"]]["generegion"] = generegion
                    matchingseqs.append(s["id"])
            dg["sequences"] = matchingseqs
            # Prepare all the lovely file names
            path_generegion_stub = path_m_dir_l + "/" + str(generegion)
            dg["base"] = path_generegion_slopped = path_generegion_stub + ".base"
            dg["base_rel"] = path_generegion_slopped_rel = path_generegion_stub + ".base_rel"
            dg["base_tight"] = path_generegion_slopped_tight = path_generegion_stub + ".base_tight"
            dg["cds_base"] = path_generegion_cds = path_generegion_stub + ".base.cds.fasta"
            dict_chr_sizes = dict_genome_info[path_genome]
            grab_bases(path_generegion, int_slop_amount, dict_chr_sizes, \
                       path_generegion_slopped, path_generegion_slopped_rel, path_generegion_slopped_tight, \
                       path_genome, path_generegion_cds, generegion)
            # Sort out any "keep" gene models.
            keepseqs = [sequence for sequence in dg["sequences"] if dict_seq_info[sequence]["keep"]]
            changeseqs = [sequence for sequence in dg["sequences"] if not sequence in keepseqs]
            if len(keepseqs) == 0:
                dg["keep"] = False
                continue
            else:
                ccc = copy.deepcopy(dg)
                if len(changeseqs) == 0:
                    del dict_generegions[generegion]
                else:
                    dict_generegions[generegion]["sequences"] = changeseqs
                    dict_generegions[generegion]["keep"] = False
                for k in keepseqs:
                    ngeneregion = "generegion_" + str(len(dict_generegions.keys()))
                    dict_generegions[ngeneregion] = copy.deepcopy(ccc)
                    dict_generegions[ngeneregion]["sequences"] = [k]
                    dict_generegions[ngeneregion]["keep"] = True
    pool.close()
    pool.join()
    # Add the cds string as a data item.
    for dg in dict_generegions.values():
        dg["cds_base_data"] = str(read_seqs(dg["cds_base"])[0].seq)
        dg["baselen"] = len(dg["cds_base_data"])
    return dict_generegions


def write_output(dict_generegions, dict_seq_info, res, path_results_dir):
    """Write out all the various bits of output
    """
    path_gtf_dir = make_if_absent(path_results_dir + "/gtf/")
    path_aa_dir = make_if_absent(path_results_dir + "/aa/")
    path_cds_dir = make_if_absent(path_results_dir + "/cds/")
    all_aa = []
    for generegion, dg in dict_generegions.items():
        for sequence in dg["sequences"]:
            gtf_id = dict_seq_info[sequence]["gtf_id"]
            gtf_out = path_gtf_dir + gtf_id + ".fixed.gtf"
            aa_out = path_aa_dir + gtf_id + ".fixed.aa.fa"
            cds_out = path_cds_dir + gtf_id + ".fixed.cds.fa"
            # Prepare gtf
            gtf = [part["gtfline"] for part in res[generegion]] if generegion in res else []
            write_csv(gtf, gtf_out + ".tmp")
            toggle_relative(gtf_out + ".tmp", gtf_out, dict_generegions[generegion]["base"], False, "", "rev")
            os.remove(gtf_out + ".tmp")
            dict_seq_info[sequence]["res_gtf"] = read_csv(gtf_out)
            dict_seq_info[sequence]["res_gtfpath"] = gtf_out
            # Prepare cds
            cds_string = "".join([dg["cds_base_data"][int(line[3]) - 1: int(line[4])] for line in gtf])
            dict_seq_info[sequence]["res_cds"] = cds_string
            write_seqs([make_seq(gtf_id, cds_string)], cds_out)
            # Prepare aa
            aa_string = translate_part(cds_string)
            dict_seq_info[sequence]["res_aa"] = aa_string
            write_seqs([make_seq(gtf_id, aa_string)], aa_out)
            all_aa.append(make_seq(gtf_id, aa_string))
    write_seqs(all_aa, path_results_dir + "/all.fa")
    align(path_results_dir + "/all.fa", path_results_dir + "/all.aln")
    return dict_generegions


def relativise_sequences(dict_generegions, dict_seq_info):
    for generegion, dg in dict_generegions.items():
        path_gtf_base = dg["base"]
        path_m_dir_l = dg["m_dirl"]
        dg["boundaries"] = path_boundaries = path_m_dir_l + "/" + generegion + ".boundaries"
        with open(path_boundaries, "w") as o:
            for seq_id in dg["sequences"]:
                path_gtf = dict_seq_info[seq_id]["gtf"]
                dict_seq_info[seq_id][
                    "gtf_rel"] = path_gtf_rel = path_m_dir_l + "/" + seq_id + ".rel." + generegion + ".gtf"
                toggle_relative(path_gtf, path_gtf_rel, path_gtf_base, True, seq_id, "fwd")
                # Can't guarantee that the gff will have frame information in it, so provide that information.
                framify(path_gtf_rel)
                with(open(path_gtf_rel, "r")) as b: o.write(b.read())


########################################
# Exonerate functions
########################################

def go_exonerate(path_w_dir, dict_generegions, dict_seq_info, int_num_cores):
    sprint("Performing first round exonerate...")
    # Prepare output directories
    path_e_dir = make_if_absent(path_w_dir + "/exonerate")  # qe
    path_e_dir1 = make_if_absent(path_e_dir + "/first_round")  # qe
    path_e_dir2 = make_if_absent(path_e_dir + "/second_round")  # qe
    path_e_dir3 = make_if_absent(path_e_dir + "/third_round")  # qe
    # Perform first round exonerate
    exonerate_first_round(dict_generegions, dict_seq_info, int_num_cores, path_e_dir1)  # qe
    # Run exonerate again using all the new potential gene parts.
    sprint("Performing second round exonerate...")
    exonerate_second_round(dict_generegions, int_num_cores, path_e_dir2)  # qe
    combine_exonerate1and2(dict_generegions, path_e_dir)  # qe
    # Run a third time.
    sprint("Performing third round exonerate...")
    exonerate_third_round(dict_generegions, dict_seq_info, path_e_dir3)  # qe
    sprint("done getting options")


def exonerate(path_base, list_all_aa_base, path_exonerate_out):
    exonerate_lines = []
    for aa in list_all_aa_base:
        aa_name = os.path.basename(aa)
        capture = subprocess.Popen(
            "exonerate --showalignment no --showtargetgff --frameshift -10000 --model protein2genome "
            + aa + " " + path_base + "  | grep \"exonerate:protein2genome:local\" | grep -P \"(cds)|(gene)\" "
            + "| grep -v \"#\" | sed -r \"s/.*exonerate:protein2genome:local/relative\\t" + aa_name + "/g\" "
            + "| sed -r \"s/cds/CDS/g\"", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout = [x for x in capture.stdout]
        stderr = [x for x in capture.stderr]
        if not stderr: 
            exonerate_lines += [a.strip("\n").split("\t") for a in stdout]
    gene_num = 0
    result = []
    for line in exonerate_lines:
        if line[2] == "gene":
            gene_num += 1
            source = re.sub(r".*sequence ([^ ;]*) ;.*", r"\1", line[8])
        nline = line[:]
        nline[8] = source
        nline[1] = "gene_" + str(gene_num) + "." + line[1]
        if line[2].lower() == "cds": 
            result.append(nline)
    write_csv(result, path_exonerate_out)


def implement_exonerate(path_cds_base, list_all_aa, path_exonerate_out):
    exonerate(path_cds_base, list_all_aa, path_exonerate_out)
    framify(path_exonerate_out)


def exonerate_first_round(dict_generegions, dict_seq_info, int_num_cores, path_e_dir):
    sprint("Exonerating and piling sequences against sequence regions...")
    for generegion, dg in dict_generegions.items():
        dg["exonerated"] = path_exonerate_out = path_e_dir + "/" + generegion + ".exonerate1"
        # There should only be one sequence in each keep generegion.
        if dg["keep"]:
            copyfile(dict_seq_info[dg["sequences"][0]]["gtf_rel"], path_exonerate_out)
        else:
            list_all_aa = [dict_seq_info[a]["aa"] for a in dict_seq_info if not a in dg["sequences"]]
            path_cds_base = dg["cds_base"]
            path_base_rel = dg["base_rel"]
            implement_exonerate(path_cds_base, list_all_aa, path_exonerate_out)
        # Merge the exonerate output
        dg["exonerated_merged"] = path_exonerate_out_merged = path_exonerate_out + ".merged"
        merge_friendly(path_exonerate_out, path_exonerate_out_merged)
        # Deal with any in-frame stops
        gtf = read_csv(path_exonerate_out_merged)
        cds = dg["cds_base_data"]
        if not dg["keep"]: 
            gtf = split_at_stops(gtf, cds)
        write_csv(gtf, path_exonerate_out_merged)
        # Split out the individual aa parts for each gtf line.
        dg["path_gene_parts_fasta"] = path_parts = path_exonerate_out_merged + ".parts"
        translate_gtf_to_file(gtf, cds, generegion, path_parts)


def exonerate_second_round(dict_generegions, int_num_cores, path_e_dir2):
    for generegion, dg in dict_generegions.items():
        dg["exonerated2"] = path_second_round = path_e_dir2 + "/" + generegion + ".exonerated2"
        path_cds = dg["cds_base"]
        if dg["keep"]:
            copyfile(dg["exonerated"], path_second_round)
        else:
            othergene_parts = [dict_generegions[ogeneregion]["path_gene_parts_fasta"] for ogeneregion in
                               dict_generegions if ogeneregion != generegion]
            implement_exonerate(path_cds, othergene_parts, path_second_round)
            path_second_round_merged = dg["exonerated2_merged"] = path_second_round + ".merged"
            merge_friendly(path_second_round, path_second_round_merged)


def exonerate_third_round(dict_generegions, dict_seq_info, path_e_dir3):
    """This function aims to clean up any regions in the exonerated output
       that overlap but are in incompatible frames.
    """
    # Generate possible versions and choose between them.
    for generegion, dg in dict_generegions.items():
        # Get possible combinations of bits of gene from the exonerate output.
        options = get_options(read_csv(dg["exonerated_all_merged"]))
        # Pick the best of those options and write them to file.
        dg["path_options"] = make_if_absent(dg["m_dirl"] + "/options/") + "/" + generegion + ".options.fasta"
        choose_options(dict_generegions[generegion], options)
    # Run exonerate on the cleaned up genes again.
    for generegion, dg in dict_generegions.items():
        options = flatten_lol(
            [[o["fasta"] for o in dict_generegions[omet]["options"].values()] for omet in dict_generegions if
             not omet == generegion])
        path_base = dict_generegions[generegion]["cds_base"]
        path_ex3 = path_e_dir3 + "/" + generegion + ".exonerated3"
        # Implement exonerate with this new set of data
        if not dg["keep"]:
            implement_exonerate(path_base, options, path_ex3)
        else:
            copyfile(dg["exonerated_merged"], path_ex3)
        gtf = read_csv(path_ex3)
        cds = dg["cds_base_data"]
        # Be careful about stops.
        if not dg["keep"]: 
            gtf = split_at_stops(gtf, cds)
        write_csv(gtf, path_ex3)
        path_all = dg["exonerated_all_merged"]
        # Add them all up.
        allgtf = read_csv(path_all) + read_csv(path_ex3)
        write_csv(allgtf, path_all)
        merge_friendly(path_all, path_all)
    # Choose between the resulting options.
    for generegion in dict_generegions:
        # Add in any bits of the input that do not overlap the exonerate-found regions.
        gtf_o = safe_gtf(flatten_lol(
            [read_csv(dict_seq_info[oseq]["gtf_rel"]) for oseq in dict_generegions[generegion]["sequences"]]))
        gtf = read_csv(dict_generegions[generegion]["exonerated_all_merged"])
        gtf = add_non_overlappers(gtf, gtf_o)
        # Hack on lazy code
        for i in gtf: i[0] = "relative"
        options = get_options(gtf)
        choose_options(dict_generegions[generegion], options)


def add_non_overlappers(list_gtf_lines, gtf_o):
    res = copy.deepcopy(list_gtf_lines)
    framed_r = split_to_abs_frame(res)
    framed_o = split_to_abs_frame(gtf_o)
    for frame in framed_o:
        for line in framed_o[frame]:
            if not any([gtf_lines_overlap(line, i) for i in framed_r[frame]]):
                res += [line]
    return safe_gtf(res)


def combine_exonerate1and2(dict_generegions, path_e_dir):
    for generegion, dg in dict_generegions.items():
        # Set up files
        dg["exonerated_all"] = path_all = path_e_dir + "/" + generegion + ".all"
        dg["exonerated_all_merged"] = path_all_merged = path_all + ".merged"
        # Perform the merge
        if dg["keep"]:
            copyfile(dg["exonerated_merged"], path_all)
        else:
            concat_files([dg["exonerated_merged"], dg["exonerated2_merged"]], path_all)
        merge_friendly(path_all, path_all_merged)
        # Split into parts and split at stops.
        gtf = read_csv(path_all_merged)
        cds = dg["cds_base_data"]
        if not dict_generegions[generegion]["keep"]:
            gtf = split_at_stops(gtf, cds)
        write_csv(gtf, path_all_merged)
        translate_gtf_to_file(gtf, cds, generegion, path_all_merged + ".parts")


def get_options(list_gtf_lines):
    xsafe, goverlaps = get_overlaps(list_gtf_lines)
    options = []
    if not goverlaps: 
        return [xsafe]
    # The first set of options should involve simply removing one or other of the offending overlappers.
    good_from_overlapping = get_non_overlapping_subsets(goverlaps)
    for suboption in good_from_overlapping:
        options.append(clean_gtf_live(xsafe + suboption))
    # The second set of options should involve chopping lines along their join point, and throwing away parts that end up too tiny.
    xoverlaps = goverlaps[:]
    for pair in itertools.combinations(goverlaps, 2):
        l_safe, l_overlaps = get_overlaps(pair)
        if l_overlaps:
            first_first = int(l_overlaps[0][3]) <= int(l_overlaps[1][3])
            first = l_overlaps[int(not first_first)]
            second = l_overlaps[int(first_first)]
            # Only want to consider things of this form:
            # ----------------
            #	   -----------------
            if int(first[4]) >= int(second[4]): 
                continue
            # Find the centroid
            overlap_start = int(second[3])
            overlap_end = int(first[4])
            breakpoint = (overlap_start + overlap_end) / 2
            # Generate new gene parts
            newfirst = first[:]
            newfirst[4] = breakpoint - 1
            newsecond = second[:]
            newsecond[3] = breakpoint
            newsecond[7] = (3 - (breakpoint - int(second[3]) - int(second[7]))) % 3
            # Don't want to include really tiny ones!!
            if min(gtf_length(newfirst), gtf_length(newsecond)) < 20: 
                continue
            xoverlaps += [newfirst, newsecond]
    good_from_overlapping = get_non_overlapping_subsets(xoverlaps)
    for suboption in good_from_overlapping:
        options.append(clean_gtf_live(xsafe + suboption))
    # The final set of options should involve simply splitting up into chunks, and throwing away parts that end up too tiny.
    #	xoverlaps=goverlaps[:]
    #	print "==============4"
    #	for pair in itertools.combinations(goverlaps, 2):
    #		l_safe, l_overlaps = get_overlaps(pair)
    #		if l_overlaps:
    #			first_first = int(l_overlaps[0][3]) <= int(l_overlaps[1][3])
    #			first = l_overlaps[int(not first_first)]
    #			second = l_overlaps[int(first_first)]
    #			firstchunk1 = first[:]
    #			firstchunk1[4] = int(second[3]) - 1
    #			firstchunk2 = first[:]
    #			firstchunk2[3] = int(second[3])
    #			firstchunk2[7] = (3 - (int(second[3]) - int(first[3]) - int(first[7]))) % 3
    #			secondchunk1 = second[:]
    #			secondchunk1[4] = int(first[4])
    #			secondchunk2 = second[:]
    #			secondchunk2[3] = int(first[4]) + 1
    #			secondchunk2[7] = (3 - (int(first[4]) + 1  - int(second[3]) - int(second[7]))) % 3
    #			for opt in [firstchunk1, firstchunk2, secondchunk1, secondchunk2]:
    #				if gtf_length(firstchunk1) >= 20: xoverlaps += [opt]
    #	good_from_overlapping = get_non_overlapping_subsets(xoverlaps)
    #	for suboption in good_from_overlapping:
    #		options.append(clean_gtf_live(xsafe + suboption))
    return options


def choose_options(dg, final):
    dg["options"] = {}
    for i, ens in enumerate(final):
        # Construct entries for each option.
        e = sorted(ens, key=lambda x: int(x[3]))
        path_option_out = dg["path_options"] + ".option_" + str(i)
        write_csv(e, path_option_out)
        # Grab the fasta sequences.
        path_option_fasta = path_option_out + ".aa.fasta"
        dg["options"][i] = {"gtf": path_option_out, "fasta": path_option_fasta, "parts": {}}
        cds = dg["cds_base_data"]
        with open(path_option_fasta, "w") as f:
            f.write(">" + dg["generegion"] + "." + str(i) + "\n")
            for j, line in enumerate(e):
                cdspart = cds[int(line[3]) + int(line[7]) - 1: int(line[4])]
                aapart = translate_part(cdspart)
                f.write(aapart)
                dg["options"][i]["parts"][j] = {"part": aapart, "partlength": len(aapart), "gtfline": line}
            f.write("\n")


########################################
# Stats
########################################

def get_stats(dict_generegions, dict_seq_info, res, path_results_dir, path_w_dir):
    # Get a gene tree.
    path_s_dir = make_if_absent(path_results_dir + "/stats/")
    path_t_dir = make_if_absent(path_results_dir + "/trees/")
    # call_function("iqtree -redo -s "+path_results_dir+"/all.aln -pre "+path_t_dir+"/all.aln.tree -m TEST -quiet")
    # Get original fasta alignment score. If there are any empty sequences, take them out.
    oseqs = read_seqs(path_w_dir + "/or.aln")
    oscore = alignment_score(oseqs, omit_empty=True)
    # Get new alignment score. If there are any empty sequences, take them out.
    nseqs = read_seqs(path_results_dir + "/all.aln")
    nscore = alignment_score(nseqs, omit_empty=True)
    # Write global stats to file:
    path_stats_global = path_s_dir + "/global.stats"
    with open(path_stats_global, "w") as f:
        f.write("oldscore\tnewscore\n")
        f.write(str(oscore) + "\t" + str(nscore))
    stats_sh = tempfile.mktemp() + ".sh"
    stats_script(stats_sh)
    # Write detailed stats to file...
    for sequence in dict_seq_info:
        ogtf = dict_seq_info[sequence]["gtf"]
        ngtf = dict_seq_info[sequence]["res_gtfpath"]
        # sprint("comparing " + ogtf + " and " + ngtf)
        sprint("Gene coordinates written to " + ngtf)
        call_function("bash " + stats_sh + " " + ogtf + " " + ngtf + " > " + path_s_dir + "/" + dict_seq_info[sequence][
            "gtf_id"] + ".stats")


def stats_script(tmpf):
    with open(tmpf, "w") as f:
        f.write(
            '#!/bin/bash\n\n# Input: the unfixed gtf, and the fixed gtf.\n\n# Possible events:\n#\texon removal - in-frame, simple\n#\texon addition - in-frame, simple\n#\tintron removal - in-frame, simple\n#\tintron addition - in-frame, simple\n#\tmoved start codon - in-frame, simple\n#\tmoved start codon - in-frame, intron-containing\n#\t\tmoved intron start - in-frame, simple\n#\t\tmoved intron end - in-frame, simple\n#\tmoved start codon - frame change,simple\n#\tmoved start codon - frame change, intron-containing.\n#\t\tmoved intron start - frame change, simple\n#\t\tmoved intron end - frame change, simple\n#\texon remains exactly the same\n#\texon addition - frame change, simple\n#\texon removal - frame change, simple\n#\tintron addition - frame change, simple\n#\tintron removal - frame change, simple\n#\t\tcomplex events (all other events).\n#\n# Output is in the format\n#\tevent type...inframe?...exoncontaining...changedbases\n\n#echo "====================="\n#echo $1\n\ngene=$1\nfixed=$2\n\not=`mktemp -d `\ngene_cds="$ot/gene_cds"\nfixed_cds="$ot/fixed_cds"\n\nawk \'$3 == "CDS"\' $gene | awk -F "\\t" \'BEGIN{OFS="\\t"}{$9="name_here"; print}\' | sort -k4,4n > $gene_cds\nawk \'$3 == "CDS"\' $fixed | awk -F "\\t" \'BEGIN{OFS="\\t"}{$9="name_here"; print}\' | sort -k4,4n > $fixed_cds\n\n# If the fixed gene is empty, just class this as gene removal and continue\n# Get length of original gene.\nif [[ ! -s $fixed_cds ]]; then \n\tremoved=`cat $gene_cds | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print 0-a}\'`\n\techo -e "removed_gene\\tinframe\\tcomplex\\t-\\t$removed"\n\texit\nfi\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\n\n# Get the total amount added and subtracted from the original.\n\nremoved=`bedtools subtract -a $gene_cds -b $fixed_cds.adj | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print 0-a}\'`\nadded=`bedtools subtract -b $gene_cds.adj -a $fixed_cds | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\necho -e "total_added\\t-\\t-\\t$added"\necho -e "total_removed\\t-\\t-\\t$removed"\n\nstrand=`head -n1 $gene_cds | cut -f7`\n\nif [[ $strand == "-" ]]; then\n\tmaxval=`cat $gene_cds $fixed_cds | cut -f4,5 | sed -r "s/\\t/\\n/g" | sort -n | tail -n1`\n\tawk -v m="$maxval" \'BEGIN{OFS="\\t"} {b=m+1-$4; a=m+1-$5; $4=a; $5=b; print}\' $gene_cds | sort -k4,4n > $gene_cds.tmp\n\tawk -v m="$maxval" \'BEGIN{OFS="\\t"} {b=m+1-$4; a=m+1-$5; $4=a; $5=b; print}\' $fixed_cds | sort -k4,4n > $fixed_cds.tmp\n\tmv $fixed_cds.tmp $fixed_cds\n\tmv $gene_cds.tmp $gene_cds\nfi\n\n# Do the start codon stuff first\n\nstarto=`head -n1 $gene_cds | cut -f 4`\nstartx=`head -n1 $fixed_cds | cut -f 4`\n\n# Check for exons that have remained exactly the same\nbedtools intersect -b $gene_cds  -a $fixed_cds -f 1 -F 1 | sed -r "s/.*/nochange\\tinframe\\tsimple\\t0/g"\nbedtools intersect -b $gene_cds  -a $fixed_cds -f 1 -F 1  > $ot/identical\n\nif [[ $starto -ne $startx ]]; then\n\tstart_chunk="$ot/start_chunk"\n\t# For bedtools and its annoying subtraction problem\n\tif [[ $starto -ge $startx ]]; then\n\t\thead -n1 $gene_cds | awk -v a=$starto -v b=$startx \'BEGIN{OFS="\\t"} {$4 = b; $5 = a -1; print}\' > $start_chunk\n\t\tbedtools intersect -b $start_chunk -a $fixed_cds > $start_chunk.is\n\t\tawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $start_chunk > $start_chunk.adj\n\t\tbedtools subtract -a $fixed_cds -b $start_chunk.adj > $fixed_cds.headless\n\t\tbedtools subtract -a $gene_cds -b $start_chunk.adj > $gene_cds.headless\n\t\tchange=`cat $start_chunk.is | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\n\t\t# Check frame\n\t\tif [[ $diff_frame -eq 0 ]]; then\n\t\t\tmessage1="moved_start\\tinframe"\n\t\telse\n\t\t\tmessage1="moved_start\\tframeshift"\n\t\tfi\n\t\t# Check whether any introns have been introduced by the new start codon.\n\t\tif cmp -s "$start_chunk.is" "$start_chunk"; then\n\t\t\tmessage2="simple"\n\t\telse\n\t\t\tmessage2="introncontaining"\n\t\tfi\n\t\techo -e "$message1\\t$message2\\t$change"\n\telse\n\t\thead -n1 $gene_cds | awk -v a=$starto -v b=$startx \'BEGIN{OFS="\\t"}{$4 = a; $5 = b -1; print}\' > $start_chunk\n\t\tbedtools intersect -b $start_chunk -a $gene_cds > $start_chunk.is\n\t\tawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $start_chunk > $start_chunk.adj\n\t\tbedtools subtract -a $fixed_cds -b $start_chunk.adj > $fixed_cds.headless\n\t\tbedtools subtract -a $gene_cds -b $start_chunk.adj > $gene_cds.headless\n\t\tchange=`cat $start_chunk.is | awk \'BEGIN{a=0} {a=a+$5-$4+1} END{print a}\'`\n\t\tif [[ $diff_frame -eq 0 ]]; then\n\t\t\tmessage1="moved_start\\tinframe"\n\t\telse\n\t\t\tmessage1="moved_start\\tframeshift"\n\t\tfi\n\t\t# Check whether any introns have been introduced by the new start codon.\n\t\tif cmp -s "$start_chunk.is" "$start_chunk"; then\n\t\t\tmessage2="simple"\n\t\telse\n\t\t\tmessage2="introncontaining"\n\t\tfi\n\t\techo -e "$message1\\t$message2\\t-$change"\n\tfi\n\tgene_cds=$gene_cds.headless\n\tfixed_cds=$fixed_cds.headless\nfi\n\n# Check for simple exon removal events\nloneexons_o="$ot/lone_o"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\nbedtools subtract -A -a $gene_cds -b $fixed_cds.adj > $loneexons_o\nawk \'{a=$5 - $4+ 1; if(a % 3 == 0) {print "removed_exon\\tinframe\\tsimple\\t-"a} else {print "removed_exon\\tframeshift\\tsimple\\t-"a}}\' $loneexons_o\n\n# Check for simple exon addition  events\nloneexons_x="$ot/lone_x"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nbedtools subtract -A -b $gene_cds.adj -a $fixed_cds > $loneexons_x\nawk \'{a=$5 - $4+ 1; if(a % 3 == 0) {print "added_exon\\tinframe\\tsimple\\t"a} else {print "added_exon\\tframeshift\\tsimple\\t"a}}\' $loneexons_x\n\n# Adjust our copy of the x\'d file once we\'ve acknowledged the changes.\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $loneexons_x > $loneexons_x.adj\nbedtools subtract -A -a $fixed_cds -b $loneexons_x.adj > $fixed_cds.tmp\ncat $fixed_cds.tmp $loneexons_o > $fixed_cds\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\'  $ot/identical >  $ot/identical.adj\nbedtools subtract -A -a  $fixed_cds -b $ot/identical.adj > $fixed_cds.tmp\nbedtools subtract -A -a  $gene_cds -b $ot/identical.adj > $gene_cds.tmp\nmv $fixed_cds.tmp $fixed_cds\nmv $gene_cds.tmp $gene_cds\n\n#For intron checking -- invert.\n#cat $fixed_cds\n#cat $gene_cds\n\n# Check for simple intron removal events\nbase_gene="$ot/base_gene"\nbase_fixed="$ot/base_fixed"\n\nstarto=`sort -k4,4n $gene_cds | head -n1 | cut -f 4`\nendo=`sort -k4,4n $gene_cds | tail -n1 | cut -f 5`\nstartx=`sort -k4,4n $fixed_cds | head -n1 | cut -f 4`\nendx=`sort -k4,4n $fixed_cds | tail -n1 | cut -f 5`\n\nhead -n1 $gene_cds | awk -v a=$starto -v b=$endo \'BEGIN{OFS="\\t"} {$4 = a; $5 = b; print}\' > $base_gene\nhead -n1 $fixed_cds | awk -v a=$startx -v b=$endx \'BEGIN{OFS="\\t"} {$4 = a; $5 = b; print}\' > $base_fixed\n\ngene_introns="$ot/gene_introns"\nfixed_introns="$ot/fixed_introns"\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds > $gene_cds.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds > $fixed_cds.adj\n\nbedtools subtract -a $base_gene -b $gene_cds.adj > $gene_introns\nbedtools subtract -a $base_fixed -b $fixed_cds.adj > $fixed_introns\n\n# Check for simple intron addition events\nloneintrons_o="$ot/lonei_o"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_introns > $fixed_introns.adj\nbedtools subtract -A -a $gene_introns -b $fixed_introns.adj > $loneintrons_o\nawk \'{a=$5 - $4 + 1; if(a % 3 == 0) {print "removed_intron\\tinframe\\tsimple\\t"a} else {print "removed_intron\\tframeshift\\tsimple\\t"a}}\' $loneintrons_o\n\n# Check for simple intron removal events\nloneintrons_x="$ot/lonei_x"\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_introns > $gene_introns.adj\nbedtools subtract -A -b $gene_introns.adj -a $fixed_introns > $loneintrons_x\nawk \'{a=$5 - $4 + 1; if(a % 3 == 0) {print "added_intron\\tinframe\\tsimple\\t-"a} else {print "added_intron\\tframeshift\\tsimple\\t-"a}}\' $loneintrons_x\n\n# Check for complex events and adjust our copy of the x\'d file once we\'ve acknowledged the changes.\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $loneintrons_x > $loneintrons_x.adj\nbedtools subtract -A -a $fixed_introns -b $loneintrons_x.adj > $fixed_introns.tmp\ncat $fixed_introns.tmp $loneintrons_o > $fixed_introns\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_introns > $fixed_introns.adj\nbedtools subtract -a  $base_fixed -b $fixed_introns.adj > $fixed_cds\n\nbedtools intersect -a $gene_cds -b $fixed_cds -wa -wb > $ot/cds_intersection\ncut -f1-9 $ot/cds_intersection | sort | uniq -u > $gene_cds.good\ncut -f1-9 $ot/cds_intersection | sort | uniq -d > $gene_cds.junk\ncut -f10-18 $ot/cds_intersection | sort | uniq -u > $fixed_cds.good\ncut -f10-18 $ot/cds_intersection | sort | uniq -d > $fixed_cds.junk\n\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds.good > $fixed_cds.good.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds.good > $gene_cds.good.adj\nbedtools subtract -b $gene_cds.good.adj -a $fixed_cds.good | awk \'{a=$5 - $4 +1; if(a % 3 == 0) {print "exon_extension\\tinframe\\tsimple\\t"a} else {print "exon_extension\\tframeshift\\tsimple\\t"a}}\'\nbedtools subtract -a $gene_cds.good -b $fixed_cds.good.adj | awk \'{a=$5 - $4 +1; if(a % 3 == 0) {print "exon_contraction\\tinframe\\tsimple\\t-"a} else {print "exon_contraction\\tframeshift\\tsimple\\t-"a}}\'\n\n# Junk\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $gene_cds.junk  > $gene_cds.junk.adj\nawk \'BEGIN{OFS="\\t"} {$5 = $5 + 1; print}\' $fixed_cds.junk > $fixed_cds.junk.adj\njunk_removed=`bedtools subtract -a $gene_cds.junk.adj -b $fixed_cds.junk | awk \'BEGIN{a=0} {a=a+$5-$4 +1} END{print 0-a}\'`\njunk_added=`bedtools subtract -b $gene_cds.junk -a $fixed_cds.junk.adj | awk \'BEGIN{a=0} {a=a+$5-$4 +1} END{print a}\'`\necho -e "other_added\\t-\\t-\\t$junk_added"\necho -e "other_removed\\t-\\t-\\t$junk_removed"\n\n#cat $gene_cds\n#cat $fixed_cds\n\n\nrm -r $ot\n')


########################################
# Choice function
########################################

def most_coherent_set(sequences, path_w_dir):
    # We take consensus slices to have something
    # to compare things to.
    sequences = list(sequences)
    sbm = group_sequences_by_generegion(sequences)
    css, order = get_best_slices(sbm)
    write_slices(css, order, path_w_dir + "/slices.aln")
    # Remember where all the sequences are.
    generegionlookup, all_sequences = get_generegion_lookup(order, sbm)
    checkpoints = get_checkpoints(generegionlookup)
    subsets = get_subset_strings(css, all_sequences, generegionlookup)
    winners = random_best_subsets(subsets, checkpoints)
    # If there are multiple winning sets, pick the best ones.
    options = {}
    best_score = -1000
    best_id = 0
    for winner in winners:
        winningseqs = {}
        for i, w in enumerate(bin(winner).split("b")[1].zfill(len(sequences))):
            if int(w) == 1:
                k = order[generegionlookup[i]]
                winningseqs[k] = winningseqs.get(k, []) + [all_sequences[i]]
        if len(winners) == 1 and all([len(winningseqs[k]) == 1 for k in winningseqs]):
            return dict((k, winningseqs[k][0]) for k in winningseqs)
        else:
            possibilities = itertools.product(*winningseqs.values())
            for p in possibilities:
                score = alignment_score(p)
                option_id = len(options)
                options[option_id] = {"seqs": dict((order[i], p[i]) for i in range(len(order))), "score": score}
                if score > best_score:
                    best_id = option_id
                    best_score = score
    return options[best_id]["seqs"]


def write_slices(css, order, path_out):
    sequences = dict((k, "".join([j[0][i] for j in css])) for i, k in enumerate(order))
    write_seqs([make_seq(k, sequences[k]) for k in sequences], path_out)


def random_best_subsets(subsets, checkpoints):
    full = "1" * len(subsets[0])
    ss = [int(s, 2) for s in subsets if not s == full]
    if not ss:
        return [int(full, 2)]
    checks_b = [int(s, 2) for s in checkpoints]
    results = []
    for i in range(1000):
        # Pick a random starting point. This will be naturally
        # weighted by the most abundant entries.
        st = ss[:]
        while True:
            r = random.choice(st)
            st = [bin_and(r, i) for i in st if bin_compat(i, r, checks_b)]
            sq = [s for s in st if not s == r]
            if sq:
                st = sq[:]
            else:
                results += st
                break
    ress = set(results)
    maxsup = 0
    winners = []
    for i in ress:
        sup = support(i, ss)
        if sup > maxsup:
            winners = [i]
            maxsup = sup
        elif sup == maxsup:
            winners += [i]
        return winners


def get_generegion_lookup(order, sbm):
    all_sequences = []
    generegionlookup = {}
    for i, k in enumerate(order):
        for sequence in sbm[k]:
            all_sequences.append(sequence)
            generegionlookup[len(all_sequences) - 1] = i
    return generegionlookup, all_sequences


def get_subset_strings(css, all_sequences, generegionlookup):
    subsets = []
    for posi, pos in enumerate(css):
        for option in pos:
            posstring = ""
            for i, sequence in enumerate(all_sequences):
                aa = sequence[posi]
                cs_aa = option[generegionlookup[i]]
                posstring += str(int(aa == cs_aa))
            subsets += [posstring]
    return subsets


def get_checkpoints(generegionlookup):
    checkpoints = dict((k, ["0"] * len(generegionlookup)) for k in generegionlookup.values())
    for k in generegionlookup:
        checkpoints[generegionlookup[k]][k] = "1"
    checkpoints = ["".join(a) for a in checkpoints.values()]
    return checkpoints


def slice_alignment(sbm):
    l = len(sbm[sbm.keys()[0]][0])
    dicty = dict((k, {}) for k in sbm)
    for i in range(0, l):
        for d in sbm:
            dicty[d][i] = [k[i] for k in sbm[d]]
    return dicty


def group_sequences_by_generegion(sequences):
    dicty = {}
    for s in sequences:
        generegion = re.sub(r"(generegion_[0-9]*).*", r"\1", s.id)
        dicty[generegion] = dicty.get(generegion, []) + [s]
    return dicty


def get_best_slices(sbm):
    ssm = slice_alignment(sbm)
    css = []
    slen = len(sbm[sbm.keys()[0]][0])
    order = [k for k in ssm]
    for i in range(0, slen):
        sslice = [ssm[k][i] for k in order]
        c = consensus(sslice)
        if not c:
            sprint("this slice doesn't yield>..")
            sprint(sslice)
            sys.exit()
        css.append(consensus(sslice))
    return css, order


def consensus(slices):
    slc = [set(i) for i in slices]
    cp = consensus_patterns(slc)
    cpall = "".join(cp)
    slc = [[a for a in x if a in cpall] for x in slc]
    return grab_representatives(slc, cp)


def consensus_patterns(slices):
    """Choose the best subset, one from each slice.
       Input is a list of lists of AAs.
    """
    slicetypes = {}
    for s in slices:
        slicetype = sorted(set(s))
        slicestring = "".join(slicetype)
        slicetypes[slicestring] = slicetypes.get(slicestring, 0) + 1
    # Need to check whether these are equivalent for my purposes
    st_choices = dict((s, [slicetypes[s] * i for i in s]) for s in slicetypes)
    # Get all the choices.
    choices = set(["".join(sorted("".join(c))) for c in itertools.product(*st_choices.values())])
    scores = dict((c, col_score(c)) for c in choices)
    maxscore = max(scores.values())
    maxscorers = [s for s in scores if scores[s] == maxscore]
    return maxscorers


def get_or_aln(dict_generegions, dict_seq_info, path_or_fa, path_or_aln):
    write_seqs(gather_sequences(dict_generegions, dict_seq_info), path_or_fa)
    align(path_or_fa, path_or_aln)


def gather_sequences(dict_generegions, dict_seq_info):
    for generegion, dg in dict_generegions.items():
        for sequence in dg["sequences"]:
            yield make_seq("original_" + generegion + "." + sequence, dict_seq_info[sequence]["sequence_aa"])


def grab_representatives(slc, cp):
    reps = Counter("".join(["".join(s) for s in slc]))
    ret = []
    for c in cp:
        poss = Counter(c)
        essential = []
        for s in slc:
            res = s
            for x in s:
                if reps[x] == poss[x]: 
                    res = [x]
            essential.append(res)
        # I haven't been able to construct an example where the product of /essential is
        # not the set we want. I'm sure they must exist, but for now I'm just going to
        # keep them all (double checking that they match the required pattern), and leave
        # the special cases to their special ways.
        for i in itertools.product(*essential):
            # Note the consensus patterns must be sorted
            if "".join(sorted(i)) == c:
                ret.append(i)
            else:
                pass
    return ret


########################################
# Part combinations
########################################

def get_combinations_fasta(partslist, combis, cds, generegion, minintron, minexon, check_start=True):
    labels = {}
    if not partslist: 
        return labels
    for j, c in enumerate(combis):
        partstring = ""
        so_cds = ""
        clist = []
        parts = []
        for i in range(0, len(c) / 2):
            part = partslist[i]
            gtf = part["gtfline"]
            nb = c[2 * i]
            ne = c[2 * i + 1]
            partc = copy.deepcopy(part)
            partc["gtfline"][3] = nb
            partc["gtfline"][4] = ne
            # Have to extract the frame carefully here
            # partc["gtfline"][7] = [a for a in partc["options_left"] if c[2*i] in partc["options_left"][a]][0]
            partc["gtfline"][7] = (3 - len(so_cds)) % 3
            clist += [nb, ne]
            parts.append(partc)
            # Add to the cumulative gene
            so_cds += cds[int(nb) - 1: int(ne)]
            partstring += "b" + str(nb) + "e" + str(ne)
        # The beginning might not be a start codon, so the frame isn't necessarily zero.
        # Use the frame of the gtf part to figure out the frame of the variant.
        firstpos = clist[0]
        firstpart = partslist[0]
        localframe = [a for a in firstpart["options_left"] if firstpos in firstpart["options_left"][a]][0]
        aa = translate_part(so_cds[localframe:])
        # Do a bunch of checks to make extra sre the gene is well-formed
        # Check that:
        # The gtf parts don't overlap
        # There are no stop codons that aren't at the end of the gene.
        # The thing has a viable start codon (if the first exon becomes very short the start codon can get destroyed).
        start_check = (not check_start) or (is_start_codon(so_cds[0:3]) and localframe == 0)
        contains_n = "n" in so_cds[localframe:].lower()
        no_stop_check = not "*" in aa[:len(aa) - 2]
        well_ordered = all(clist[p] <= clist[p + 1] for p in range(len(clist) - 1))
        good_introns = all(clist[(2 * p) + 1] + minintron < clist[2 * (p + 1)] for p in range(len(clist) / 2 - 1))
        # If the gene has made the cut, then happy days: continue.
        if well_ordered and good_introns and no_stop_check and start_check and not contains_n:
            label = generegion + ".superopt_" + str(j) + "." + partstring
            labels[label] = {"parts": parts, "aa": aa, "seq": make_seq(label, aa), "generegion": generegion}
    return labels


def get_part_combinations(partslist, cds, check_init_term=False, reject=[]):
    boundaries = flatten_lol([[part["left_flat"], part["right_flat"]] for part in partslist])
    boundaries = [[a for a in b if not a in reject] for b in boundaries]
    combis = list(itertools.product(*boundaries))
    combiscopy = combis[:]
    for c in combis:
        for i in range(0, len(c) / 2 - 1):
            donor = partslist[i]
            acceptor = partslist[i + 1]
            framesd = [a for a in donor["options_right"] if c[2 * i + 1] in donor["options_right"][a]]
            framesr = [a for a in acceptor["options_left"] if c[2 * i + 2] in acceptor["options_left"][a]]
            # Sometimes we may try to fuse a non-donor or non-acceptor (i.e. an initial or terminal
            # exon) to something else. We can't do this.
            if not (is_donor_site(c[2 * i + 1], cds) and is_acceptor_site(c[2 * i + 2], cds)):
                if c in combiscopy: 
                    combiscopy.remove(c)
            if not any((framed + framer) % 3 == 0 for (framed, framer) in itertools.product(framesd, framesr)):
                if c in combiscopy: 
                    combiscopy.remove(c)
    return combiscopy


def fetch_by_label(lbls, wnrs):
    return dict((gr, [] if wnrs[gr].id == gr + ".blank" else lbls[gr][wnrs[gr].id]["parts"]) for gr in wnrs)


########################################
# Basic fixing
########################################

def fix_it(adj, parts, d_gr, path_w_dir, minintron, minexon, path_winners_aln):
    """This function repeatedly adds in small chunks of gene and picks the best combination at
       each stage. Each stage is wiggled recursively to aim for optimal results.
    """
    # Run through each set of adjacent gene parts, sequentially adding them in.
    res = idict(d_gr, [])
    for a, ar in enumerate(adj):
        for r in ar:
            generegion, partid = parse_id(r)
            if not generegion in res: 
                res[generegion] = []
            latestpart = parts[generegion][partid]
            already_there = [x for x in res[generegion] if int(x["id"]) == int(partid)]
            if already_there:
                already_there[0]["status"] = "fresh"
            if not already_there:
                latestpart["status"] = "fresh"
                # Make sure the previous part is waiting.
                # Can't be terminal
                if res[generegion]:
                    res[generegion][-1]["status"] = "waiting_a"
                    res[generegion][-1]["terminal"] = False
                res[generegion].append(latestpart)
            first_part_id = min([int(x["id"]) for x in res[generegion]])
            for x in res[generegion]:
                if x["id"] == first_part_id:
                    x["initial"] = True
        path_fix = make_if_absent(path_w_dir + "/fix/fix_" + str(a))
        sprint("Fixing " + str(a))
        res = incremental_fix_recursive(res, path_fix, minintron, minexon, path_winners_aln, d_gr)
    sprint("finishing genes.")
    path_final = make_if_absent(path_w_dir + "/fix/fix_final")
    res = incremental_fix_recursive(res, path_final, minintron, minexon, path_winners_aln, d_gr, t_terminal=True)
    return res


def incremental_fix(parts, p_l_dir, mi, mx, path_winners_aln, d_gr, refine=False, t_terminal=False):
    """Takes a set of parts and wiggles the ends around.
    """
    path_fasta = p_l_dir + "/options.all.fasta"
    # Wiggle each generegion separately and then consider together at the end.
    labels, options = prepare_part_sets(parts, path_winners_aln, t_terminal, mi, mx, p_l_dir, d_gr)
    # Finally, compare the options to choose the winners.
    return process_labels(labels, options, p_l_dir, path_winners_aln, False, refine)


def incremental_fix_recursive(parts, path_f_dir, mi, mx, path_winners_aln, d_gr, t_terminal=False):
    iteration = 0
    iterate = True
    inparts = copy.deepcopy(parts)
    results = idict(parts, [])
    stringsets = []
    resses = {}
    prev_alns = {}
    while iterate:
        path_i_dir = make_if_absent(path_f_dir + "/iteration_" + str(iteration))
        res, prev_aln = incremental_fix(inparts, path_i_dir, mi, mx, path_winners_aln, d_gr, t_terminal=t_terminal)
        partstrings = dict((a, get_part_string(a, res[a])) for a in res)
        for i, e in enumerate(stringsets):
            if all([partstrings[k] == e[k] for k in partstrings]):
                iterate = False
                choice_indices = range(0, iteration)
        resses[iteration] = res
        prev_alns[iteration] = prev_aln + ".re"
        align(prev_aln, prev_aln + ".re")
        stringsets += [partstrings]
        iteration += 1
        inparts = copy.deepcopy(res)
    # Pick whichever res has the best alignment score.
    # If there's only one, just return that.
    winner = choose_iteration(prev_alns, choice_indices)
    res = refine_statuses(resses[winner])
    return res


def status_check(data, list_good):
    # Can accept both strings and part lists...
    if type(data) == dict and "status" in data:
        status = data["status"]
    elif type(data) == str:
        status = data
    else:
        raise Exception
    return any(status.startswith(a) for a in list_good)


def kill_bad_parts(res, parts_list):
    parts = copy.deepcopy(res)
    for part in parts:
        if not (part["left_flat"] and part["right_flat"]):
            parts_list = [a for a in parts_list if int(a["id"]) != int(part["id"])]
            prev_ids = [int(a["id"]) for a in parts_list if int(a["id"]) < int(part["id"])]
            if prev_ids:
                for opart in parts_list:
                    if opart["id"] == max(prev_ids):
                        opart["status"] = "waiting_c"
                        opart["gtfline"][4] = opart["ogtfline"][4]
    return parts, parts_list


def scan_left(i, parts_list, gtf, part, status, cds):
    # If the status is fresh, lock the RHS.
    if status_check(status, ["fresh"]): 
        part["options_right"][get_end_frame(gtf)] = [int(gtf[4])]
    # If the part is at the beginning, look for a start codon.
    if i == 0:
        part["options_left"][0] += find_starts(get_protoexon(gtf), cds, True, True, True, False)
    # Otherwise, look for an acceptor
    else:
        part["options_left"] = \
            find_splice_sites(get_protoexon(gtf), cds, left=True, right=False, happy_with_central=False)[0]


def scan_right(i, parts_list, gtf, part, status, cds, t_terminal):
    # Try to find donors, too.
    if status_check(status, ["waiting"]) and not (t_terminal and i == len(parts_list) - 1):
        part["options_right"] = \
            find_splice_sites(get_protoexon(gtf), cds, left=False, right=True, happy_with_central=False)[1]
    # If this is the last part, also try adding in a terminal option.
    if status_check(status, ["waiting"]) and i == len(parts_list) - 1:
        stops = find_stops(get_protoexon(gtf), cds)
        if t_terminal: 
            part["options_right"] = idict([0, 1, 2], [])
        if stops: 
            part["options_right"][0] += [stops[0]]
        # Allow for microstops. We can't tell at this point but we'll double check later.
        if gtf_length(gtf) <= 3: 
            part["options_right"][0] += [int(gtf[4]) + 1]


def prepare_parts(parts_list_in, cds, t_terminal=False):
    parts_list = copy.deepcopy(parts_list_in)
    finished = False
    while not finished:
        res = []
        for i, part in enumerate(parts_list):
            part = initialise_options(part)
            gtf = part["gtfline"]
            status = part["status"]
            if not "ogtfline" in part.keys(): 
                part["ogtfline"] = gtf[:]
            if status_check(status, ["fresh", "double"]): 
                  scan_left(i, parts_list, gtf, part, status, cds)
            if status_check(status, ["waiting", "double"]): 
                scan_right(i, parts_list, gtf, part, status, cds,
                                                                       t_terminal)
            # flatten options and add to result.
            res.append(flatten_options(part))
        # Some parts will prevent the gene from being constructed properly. Check for these.
        # If we can't find boundaries for a part, we have to kill it.
        if all([part["left_flat"] and part["right_flat"] for part in res]):
            finished = True
        else:
            res, parts_list = kill_bad_parts(res, parts_list)
    return res


def get_reworked_combinations(plist, cds, generegion, minintron, minexon, t_terminal=False):
    plist = prepare_parts(plist, cds, t_terminal)
    return get_combinations_fasta(plist, get_part_combinations(plist, cds), cds, generegion, minintron,
                                  minexon) if plist else {}


def no_fix(parts):
    pc = parts[:]
    for p in pc: p["status"] = "done"
    return pc


def tinypart(plist):
    # Check if any of the most recent parts are too small to be useful.
    if len(plist) < 2: 
        return False
    lastpartgtf = plist[-2]["gtfline"]
    return (int(lastpartgtf[4]) - int(lastpartgtf[3]) < 40)


def check_integrity(plist, cds, gr, mi, mx, labels_l, tryfix, t_terminal):
    if tryfix and ((not labels_l) or tinypart(plist)):
        plist1, plist2 = rework_ends(plist)
        labels_l1 = get_reworked_combinations(plist1, cds, gr, mi, mx, t_terminal=t_terminal)
        labels_l2 = get_reworked_combinations(plist2, cds, gr, mi, mx, t_terminal=t_terminal)
        for l in labels_l1: labels_l[l] = labels_l1[l]
        for l in labels_l2: labels_l[l] = labels_l2[l]
    return labels_l


def check_stops(plist, cds, gr, mi, mx, labels_l):
    plist1 = copy.deepcopy(plist[:-1])
    plist1[-1]["terminal"] = True
    plist1[-1]["status"] = "waiting_q"
    plist1 = prepare_parts(plist1, cds, t_terminal=True)
    labels_t = get_reworked_combinations(plist1, cds, gr, mi, mx, t_terminal=True)
    for l in labels_t: labels_l[l] = labels_t[l]
    return labels_l


def prepare_part_sets(parts, path_winners_aln, t_terminal, minintron, minexon, path_f_dir, d_gr):
    labels = idict(parts, {})
    options = idict(parts, [])
    write(parts, path_f_dir + "/prevparts.tpy")
    for gr in parts:
        if parts[gr]:
            cds = d_gr[gr]["cds_base_data"]
            # Prepare the initial round of parts
            tryfix = not d_gr[gr]["keep"]
            plist = prepare_parts(parts[gr] if tryfix else no_fix(parts[gr]), cds, t_terminal=t_terminal)
            write(plist, path_f_dir + "/plist.tpy")
            if plist:
                write_seqs([make_seq(str(p["id"]), p["part"]) for p in parts[gr]],
                           path_f_dir + "/parts." + gr + ".fasta")
                # Now get some part combinations
                labels_l = get_combinations_fasta(plist, get_part_combinations(plist, cds), cds, gr, minintron, minexon)
                # Check if any of the combis are empty. If so, try and rearrange the most recent bits.
                # At best we should be able to get the previous lot back.
                labels_l = check_integrity(plist, cds, gr, minintron, minexon, labels_l, tryfix, t_terminal=t_terminal)
                # If this is the terminal step, double check that removing the last part of
                # each gene does not improve things.
                if t_terminal and len(plist) > 1: 
                    labels_l = check_stops(plist, cds, gr, minintron, minexon, labels_l)
                # Get the options
                labels[gr], options[gr] = labels_l, [l["seq"] for l in labels_l.values()]
    return labels, options


def refine_statuses(parts):
    partsnew = idict(parts, [])
    for generegion in parts:
        parts_dict = dict((x["id"], x) for x in parts[generegion])
        part_ids = [int(x["id"]) for x in parts[generegion]]
        for i, part_id in enumerate(part_ids):
            part = parts_dict[part_id]
            part["status"] = "waiting_d" if status_check(part, ["fresh", "double"]) or (
                    not part["terminal"] and int(part_id) == max(part_ids)) else "done"
            partsnew[generegion].append(part)
    return partsnew


def set_status_by_id(plist, pid, status):
    for p in plist:
        if p["id"] == pid: 
            p["status"] = status


def set_terminal_by_id(plist, pid):
    for p in plist:
        if p["id"] == pid: 
            p["terminal"] = True


def rework_ends(partslist):
    # Release two variants on the partslist
    # 1. Removing the penultimate part
    # 2. Removing the last part.
    p1 = copy.deepcopy(partslist)
    p2 = copy.deepcopy(partslist)
    fresh_ids = [a["id"] for a in partslist if status_check(a, ["fresh"])]
    # We can find ourselves without any freshids, if the end reworking is taking place
    # after a terminal part has been removed.
    last_id = int(fresh_ids[0]) if fresh_ids else max([int(a["id"]) for a in partslist])
    prev_ids = [int(a["id"]) for a in partslist if int(a["id"]) < last_id]
    # If there is no previous id, simply delete this gene part and return empty.
    if not prev_ids: 
        return [], []
    # If there is only one previous id, try both omitting the fresh id and omitting the original.
    # In this case we have to set the fresh id to initial.
    prev_id = max(prev_ids)
    if len(prev_ids) == 1:
        p1 = [v for v in p1 if not int(v["id"]) == prev_id]
        set_status_by_id(p1, prev_id, "double_z")
        p2 = [v for v in p1 if not int(v["id"]) == last_id]
    # If there are multiple previous ids, we must attempt to unravel the pen-penultimate one.
    elif len(prev_ids) > 1:
        p1 = [v for v in p1 if not int(v["id"]) == prev_id]
        prevprev_id = max([int(a["id"]) for a in partslist if int(a["id"]) < prev_id])
        set_status_by_id(p1, prevprev_id, "waiting_e")
        set_status_by_id(p1, last_id, "double_p")
        p2 = [v for v in p1 if not int(v["id"]) == last_id]
    o_last_part = [v for v in partslist if int(v["id"]) == last_id][0]
    if o_last_part["terminal"]:
        if p1: 
            set_terminal_by_id(p1, max([int(a["id"]) for a in p1]))
        if p2: 
            set_terminal_by_id(p2, max([int(a["id"]) for a in p2]))
    return p1, p2


def initialise_options(part):
    for bag in ["options_left", "options_right"]:
        part[bag] = idict([0, 1, 2], [])
    gtf = part["gtfline"]
    part["options_left"][int(gtf[7])] = [int(gtf[3])]
    part["options_right"][get_end_frame(gtf)] = [int(gtf[4])]
    return part


def flatten_options(part):
    sl = part["options_left"]
    sr = part["options_right"]
    part["left_flat"] = sl[0] + sl[1] + sl[2]
    part["right_flat"] = sr[0] + sr[1] + sr[2]
    return part


def choose_iteration(prev_alns, choice_indices):
    winning_score = -float("inf")
    for i in choice_indices:
        score = get_score(prev_alns[i])
        if score > winning_score:
            winner, winning_score = i, score
    return winner


########################################
# Filtering
########################################

def filter_changes(adj, c_parts, c_coords, path_fltrs, dict_cds_bases, path_aln, minintron, minexon, dict_generegions,
                   c_mcount):
    # Sort the genes into piles based on gene region.
    compare = {}
    for gene in c_parts:
        generegion = get_gr(gene)
        compare[generegion] = compare.get(generegion, []) + [gene]
    # For each gene region, find introns and exons that only exist in the
    # new gene. These come out as generegion coordinates.
    nov_ex, nov_in, nov_fl = get_novel_regions(compare, c_parts, c_coords, dict_generegions)
    # Grab the intervals
    alnseqs = read_seqs(path_aln, "fasta")
    in_ex, in_in, in_fl = get_inspection_intervals(nov_ex, nov_in, nov_fl, c_parts, c_coords, alnseqs, compare,
                                                   path_aln)
    # Filter the intervals
    reject = filter_intervals(in_ex, in_in, in_fl, c_parts, c_coords, alnseqs, compare)
    # Remove anything from the reject pile that was present in any of the original gene.
    reject = allow_originals(reject, c_parts, compare)
    return compare_parts(adj, c_parts, c_mcount, path_fltrs, dict_cds_bases, path_aln, minintron, minexon,
                         dict_generegions, d_reject=reject)


def get_intervals(nov, c_parts, compare, alnlen, c_coords, path_aln, grab):
    feature_intervals = {}
    for gene in nov:
        o_genes = [g for g in compare[get_gr(gene)] if not g == gene]
        feature_intervals[gene] = flatten_lol(
            grab(feature, c_parts, gene, o_genes, alnlen, c_coords, path_aln) for feature in nov[gene])
    return feature_intervals


def get_inspection_intervals(nov_ex, nov_in, nov_fl, c_parts, c_coords, alnseqs, compare, path_aln):
    call_function("echo \"\" > " + path_aln + ".fil")
    exon_intervals = {}
    intron_intervals = {}
    flank_intervals = {}
    alnlen = len(alnseqs[0])
    exon_intervals = get_intervals(nov_ex, c_parts, compare, alnlen, c_coords, path_aln, grab_exon_intervals)
    intron_intervals = get_intervals(nov_in, c_parts, compare, alnlen, c_coords, path_aln, grab_intron_intervals)
    flank_intervals = get_intervals(nov_fl, c_parts, compare, alnlen, c_coords, path_aln, grab_flank_intervals)
    return exon_intervals, intron_intervals, flank_intervals


def grab_exon_intervals(exon, c_parts, gene, orig_genes, alnlen, c_coords, path_aln):
    # This 100% should exist and should be unique.
    r_part_id = [a for a in c_parts[gene] if
                 int(c_parts[gene][a]["gtfline"][3]) == exon[0] and int(c_parts[gene][a]["gtfline"][4]) == exon[1]][0]
    rgtf = c_parts[gene][r_part_id]["gtfline"]
    interval_probes = []
    # Go through each of the original genes, inspecting only the junctions
    # that appear to have changed.
    # LHS first
    r_prev_ids = [a for a in c_parts[gene] if a < r_part_id]
    r_interval = [int(c_parts[gene][max(r_prev_ids)]["gtfline"][4]) if r_prev_ids else 0,
                  int(c_parts[gene][r_part_id]["gtfline"][3])]
    for o in orig_genes:
        o_part_ids = [k for k in c_parts[o] if gtf_lines_overlap(c_parts[o][k]["gtfline"], rgtf)]
        o_prev_ids = [k for k in c_parts[o] if all(k < l for l in o_part_ids)]
        o_interval = [int(c_parts[o][max(o_prev_ids)]["gtfline"][4]) if o_prev_ids else 0,
                      int(c_parts[o][min(o_part_ids)]["gtfline"][3])] if o_part_ids else []
        intervals = grab_aln_intervals_ex(o_interval, r_interval, gene, o, c_parts, c_coords, path_aln + ".fil", alnlen,
                                          "l")
        if intervals: 
            interval_probes.append(
            grab_interval_probe(r_prev_ids, [r_part_id], intervals, gene, o, flavour="ex", direction="right"))
    # RHS next.
    r_next_ids = [a for a in c_parts[gene] if a > r_part_id]
    r_interval = [int(c_parts[gene][r_part_id]["gtfline"][4]),
                  int(c_parts[gene][min(r_next_ids)]["gtfline"][3]) if r_next_ids else -1]
    for o in orig_genes:
        o_part_ids = [k for k in c_parts[o] if gtf_lines_overlap(c_parts[o][k]["gtfline"], rgtf)]
        o_next_ids = [k for k in c_parts[o] if all(k > l for l in o_part_ids)]
        o_interval = [int(c_parts[o][max(o_part_ids)]["gtfline"][4]),
                      int(c_parts[o][min(o_next_ids)]["gtfline"][3]) if o_next_ids else -1] if o_part_ids else []
        intervals = grab_aln_intervals_ex(o_interval, r_interval, gene, o, c_parts, c_coords, path_aln + ".fil", alnlen,
                                          "r")
        if intervals: 
            interval_probes.append(
            grab_interval_probe([r_part_id], r_next_ids, intervals, gene, o, flavour="ex", direction="left"))
    return interval_probes


def grab_flank_intervals(flank, c_parts, gene, orig_genes, alnlen, c_coords, path_aln):
    # Each flank represents a single end part.
    # This should 100% exist and should be unique
    cpg = c_parts[gene].keys()
    r_right = [cpg[i] for i in cpg if int(c_parts[gene][cpg[i]]["gtfline"][3]) == flank[1] + 1]
    r_left = [cpg[i] for i in cpg if int(c_parts[gene][cpg[i]]["gtfline"][4]) == flank[0] - 1]
    interval_probes = []
    if r_right:
        r_part_id = r_right[0]
        rgtf = c_parts[gene][r_part_id]["gtfline"][:]
        rgtf[3] = flank[0]
        rgtf[4] = c_parts[gene][r_part_id]["gtfline"][3] - 1
        r_interval = [rgtf[3], rgtf[4]]
        for o in orig_genes:
            oflank = int(c_parts[o][min(c_parts[o].keys())]["gtfline"][3]) - 1
            o_interval = [flank[0], oflank]
            intervals = grab_aln_intervals_fl(o_interval, r_interval, gene, o, c_parts, c_coords, path_aln + ".fil",
                                              alnlen, "r")
            if intervals: 
                interval_probes.append(
                grab_interval_probe([], [r_part_id], intervals, gene, o, flavour="fl", direction="right"))
    elif r_left:
        r_part_id = r_left[0]
        rgtf = c_parts[gene][r_part_id]["gtfline"][:]
        rgtf[4] = flank[1]
        rgtf[3] = int(c_parts[gene][r_part_id]["gtfline"][4]) + 1
        [rgtf[3], rgtf[4]]
        r_interval = [rgtf[3], rgtf[4]]
        for o in orig_genes:
            oflank = int(c_parts[o][max(c_parts[o].keys())]["gtfline"][4]) + 1
            o_interval = [oflank, flank[1]]
            intervals = grab_aln_intervals_fl(o_interval, r_interval, gene, o, c_parts, c_coords, path_aln + ".fil",
                                              alnlen, "l")
            if intervals: 
                interval_probes.append(
                grab_interval_probe([r_part_id], [], intervals, gene, o, flavour="fl", direction="left"))
    return interval_probes


def grab_intron_intervals(intron, c_parts, gene, orig_genes, alnlen, c_coords, path_aln):
    # Rather than being a single part here, each intron represents a pair of gene parts.
    # This should 100% exist and should be unique
    cpg = c_parts[gene].keys()
    r_part_ids = [[cpg[i], cpg[i + 1]] for i in cpg[:-1] if
                  int(c_parts[gene][cpg[i]]["gtfline"][4]) == intron[0] - 1 and int(
                      c_parts[gene][cpg[i + 1]]["gtfline"][3]) == intron[1] + 1][0]
    # Construct an intron gtfline for this intron.
    rgtf = c_parts[gene][r_part_ids[0]]["gtfline"][:]
    rgtf[3] = c_parts[gene][r_part_ids[0]]["gtfline"][4] + 1
    rgtf[4] = c_parts[gene][r_part_ids[1]]["gtfline"][3] - 1
    r_interval = [rgtf[3], rgtf[4]]
    interval_probes = []
    # Go through each of the original genes, inspecting only the introns
    # that appear to have changed.
    for o in orig_genes:
        introngtfs = get_introns(c_parts[o])
        o_intron_ids = [k for k in introngtfs if gtf_lines_overlap(introngtfs[k]["gtfline"], rgtf)]
        # If there are no intron ids here, just carry on.
        # This Ross will be sorted out anyway.
        if not o_intron_ids: 
            continue
        o_interval = [introngtfs[min(o_intron_ids)]["gtfline"][3], introngtfs[max(o_intron_ids)]["gtfline"][4]]
        intervals = grab_aln_intervals_in(o_interval, r_interval, gene, o, c_parts, c_coords, path_aln + ".fil", alnlen,
                                          force=len(o_intron_ids) > 1)
        if intervals: 
            interval_probes.append(
            grab_interval_probe([r_part_ids[0]], [r_part_ids[1]], intervals, gene, o, flavour="in", direction="both"))
    return interval_probes


def get_features(parts, ubound):
    part_ids = sorted(parts.keys(), key=lambda x: int(parts[x]["gtfline"][3]))
    exons = []
    introns = []
    flanks = []
    for i in range(len(part_ids) - 1):
        gtf1 = parts[i]["gtfline"]
    gtf2 = parts[i + 1]["gtfline"]
    if i == 0: 
        exons += [[int(gtf1[3]), int(gtf1[4])]]
    exons += [[int(gtf2[3]), int(gtf2[4])]]
    introns += [[int(gtf1[4]) + 1, int(gtf2[3]) - 1]]
    # Add the flanking regions. These will be processed similarly to introns.
    flanks += [[0, int(parts[part_ids[0]]["gtfline"][3]) - 1]]
    flanks += [[int(parts[part_ids[-1]]["gtfline"][4]) + 1, ubound - 1]]
    return exons, introns, flanks


def get_interval_juncs(interval_group, c_parts, resgene):
    res = []
    # Interval group should contain either one item or two.
    for i in interval_group:
        if "r" in i and i["direction"] in ["right", "both"]: 
            res += [c_parts[resgene][i["r"]]["gtfline"][3]]
        if "l" in i and i["direction"] in ["left", "both"]: 
            res += [c_parts[resgene][i["l"]]["gtfline"][4]]
    return list(set(res))


def get_interval_group(tag1, tag2, dirn, p_ids, intervals, a, ival_groups):
    partners = [i for i in intervals if tag1 in i and p_ids.index(i[tag1]) == p_ids.index(a[tag2]) + dirn \
                and i["ogene"] == a["ogene"] and not i == a]
    if partners:
        if not anyperm([a, partners[0]], ival_groups):
            return [[a, partners[0]]]
        elif not [a] in ival_groups:
            return [[a]]
    return []


def get_interval_groups(c_parts, resgene, intervals):
    # They either should all be exonic or all not.
    if any([a["flavour"] == "ex" for a in intervals]):
        interval_groups = []
        part_ids = c_parts[resgene].keys()
        for a in intervals:
            # grab any adjacent parts, to be considered together
            if "r" in a: 
                interval_groups += get_interval_group("l", "r", -1, part_ids, intervals, a, interval_groups)
            if "l" in a: 
                interval_groups += get_interval_group("r", "l", 1, part_ids, intervals, a, interval_groups)
        return interval_groups
    else:
        return de_dup([[a] for a in intervals])


def allow_originals(reject, c_parts, compare):
    # Removes any coordinates in the reject pile that existed in the original genes.
    # We only want to exclude new stuff that we've found
    res = {}
    for generegion in compare:
        if not generegion in reject: 
            continue
        orig_genes = [a for a in compare[generegion] if not "result" in a]
        res[generegion] = flatten_lol([a for a in reject[generegion].values()])
        for o in orig_genes:
            coords = flatten_lol([[int(i["gtfline"][3]), int(i["gtfline"][4]), ] for i in c_parts[o].values()])
            res[generegion] = [r for r in res[generegion] if not r in coords]
    return res


def aq(aln_nc, rgene, ogene, compare, chop):
    for grpair in itertools.combinations([a for a in compare if not a == get_gr(rgene)], 2):
        region = chop if not aln_nc or len(aln_nc[0]) == 0 else aln_nc
        oseqs = [a for a in region if (not "result" in a.id and get_gr(a.id) in grpair) or a.id == ogene]
        rseqs = [a for a in region if ("result" in a.id and get_gr(a.id) in grpair) or a.id == rgene]
        if not rseqs or not len(rseqs[0]): 
            return False
        nonzerocols = [i for i in range(len(rseqs[0])) if not all(k[i] == "-" for k in rseqs)]
        # Judge the changed and unchanged regions separately.
        og = [a for a in aln_nc if a.id == ogene][0]
        rg = [a for a in aln_nc if a.id == rgene][0]
        changed = [i for i in range(len(og)) if og[i] != rg[i]]
        unchanged = [i for i in range(len(og)) if og[i] == rg[i]]
        changed_nz = [i for i in changed if not all(k[i] == "-" for k in rseqs)]
        unchanged_nz = [i for i in unchanged if not all(k[i] == "-" for k in rseqs)]
        # The flanking alignment needs to be good as a start.
        if unchanged_nz:
            if not alignment_score(chop_alignment(aln_nc, unchanged_nz)) / len(unchanged_nz) > 3:
                continue
        # If the flanking alignment is okay, then we just need to check that the insides are okay also.
        # This is quite a strict criterion.
        if not changed_nz:
            return False
        if changed_nz:
            if alignment_score(chop_alignment(aln_nc, changed_nz)) / len(changed_nz) > 3:
                return False
    return True


def minimal_gains(aln_nc, rgene, ogene, compare):
    if not aln_nc: 
        return True
    for grpair in itertools.combinations([a for a in compare if not a == get_gr(rgene)], 2):
        # Approach is as follows. Take clumps of three genes. Get the nonconstant portions
        # of their alignment in the region of interest.
        ogenes = [a for a in aln_nc if (not "result" in a.id and get_gr(a.id) in grpair) or a.id == ogene]
        rgenes = [a for a in aln_nc if ("result" in a.id and get_gr(a.id) in grpair) or a.id == rgene]
        og = [a for a in aln_nc if a.id == ogene][0]
        rg = [a for a in aln_nc if a.id == rgene][0]
        # Grab the scores
        oscore = alignment_score(
            [a for a in aln_nc if (not "result" in a.id and get_gr(a.id) in grpair) or a.id == ogene])
        rscore = alignment_score([a for a in aln_nc if ("result" in a.id and get_gr(a.id) in grpair) or a.id == rgene])
        # If contiguous in both cases and less than or equal to 3 aa change, don't allow.
        # If contiguous and longer, require a total increase of 4 or more
        if not any((og[i] == "-" and rg[i] == "-") for i in range(len(rg))):
            if len(aln_nc[0]) < 4: 
                return True
            if rscore > oscore + 4: 
                return False
        # If noncontiguous, require a total increase of 2 or more
        else:
            if rscore > oscore + 2: 
                return False
    return True


def parallels(aln_nc, rgene, ogene, compare):
    if not aln_nc or len(aln_nc[0]) == 0: 
        return False
    for grpair in itertools.combinations([a for a in compare if not a == get_gr(rgene)], 2):
        oseqs = [a for a in aln_nc if (not "result" in a.id and get_gr(a.id) in grpair) or a.id == ogene]
        rseqs = [a for a in aln_nc if ("result" in a.id and get_gr(a.id) in grpair) or a.id == rgene]
        # Check to see if there are any columns that are all empty in the original and not in the new.
        if not oseqs: 
            continue
        for i in range(len(oseqs[0])):
            if all(a[i] == "-" for a in oseqs) and any(a[i] != "-" for a in rseqs): 
                return True
    return False


def get_novel_regions(compare, c_parts, c_coords, dict_generegions):
    # Extract all the exons and introns of the result gene as regions.
    # We look in the original gene set to see if there is anything that is present
    # in the original bunch that is not in the new gene.
    nov_ex = {}
    nov_in = {}
    nov_fl = {}
    for generegion in compare:
        results = [a for a in compare[generegion] if "result" in a]
        if not results: 
            continue
        # There should be at most one result.
        result = results[0]
        # Grab the introns and exons from the old and new.
        ubound = dict_generegions[generegion]["baselen"]
        r_exons, r_introns, r_flanks = get_features(c_parts[result], ubound)
        origs = [a for a in compare[generegion] if not a == result]
        o_features = [get_features(c_parts[o], ubound) for o in origs]
        o_exons = de_dup(flatten_lol(([o[0] for o in o_features])))
        o_introns = de_dup(flatten_lol(([o[1] for o in o_features])))
        o_flanks = de_dup(flatten_lol(([o[2] for o in o_features])))
        # Check which of these are novel.
        nov_ex[result] = [e for e in r_exons if not e in o_exons]
        nov_in[result] = [i for i in r_introns if not i in o_introns]
        nov_fl[result] = [i for i in r_flanks if not i in o_flanks]
    return nov_ex, nov_in, nov_fl


def filter_intervals(in_ex, in_in, in_fl, c_parts, c_coords, alnseqs, compare):
    # Filter the intervals based on various criteria.
    reject = {}
    dontreject = {}
    for in_feature in [in_ex, in_in, in_fl]:
        for resgene in in_feature:
            # If any the features are exonic, we try to group them into
            # pairs that we can compare together. Otherwise consider them individually.
            intervals = in_feature[resgene]
            if not intervals: 
                continue
            interval_groups = get_interval_groups(c_parts, resgene, intervals)
            # Now consider each of the interval groups. There should be
            # at most two items in each interval group. When there are two
            # items in an interval group, there will be two parts to the junction, so
            # consider them both.
            for i in interval_groups:
                ival = list(set(flatten_lol([interv["interval"] for interv in i])))
                rjuncs = get_interval_juncs(i, c_parts, resgene)
                chop = chop_alignment(alnseqs, ival)
                generegion = re.sub(r"(generegion_[0-9]*).*", r"\1", resgene)
                # Make containers for the results of the different tests.
                if not generegion in reject: 
                    reject[generegion] = {}
                if not generegion in dontreject: 
                    dontreject[generegion] = {}
                ogene = i[0]["ogene"]
                ###############################
                # Alignment improvement check.
                ###############################
                # Grab the alignment scores for the original alignment and the new
                seqo = [a for a in alnseqs if a.id == ogene][0]
                seqr = [a for a in alnseqs if a.id == resgene][0]
                nonconstantregion = [a for a in ival if 0 <= a < len(seqo) and seqo[a] != seqr[a]]
                nonconstantslopped = flatten_lol([range(a - 10, a + 10) for a in nonconstantregion])
                aln_nc_slop = chop_alignment(alnseqs, sorted(set(nonconstantslopped)))
                oldscore = alignment_score([a for a in aln_nc_slop if not "result" in a.id and gr_pick(a.id, ogene)])
                newscore = alignment_score([a for a in aln_nc_slop if "result" in a.id])
                # Check whether the new score is better than the old score.
                # Form specific rejections groups for specific tests.
                if not "alnimprove" in reject[generegion]: 
                    reject[generegion]["alnimprove"] = []
                if not "alnimprove" in dontreject[generegion]: 
                    dontreject[generegion]["alnimprove"] = []
                if not oldscore <= newscore:
                    reject[generegion]["alnimprove"] += rjuncs
                else:
                    dontreject[generegion]["alnimprove"] += rjuncs
                ################################
                # Alignment quality check.
                ################################
                # Check the alignment quality of the region of interest. If the
                # alignment isn't very good to start with, or if the resulting alignment
                # isn't very good, reject.
                if not "aq" in reject[generegion]: 
                    reject[generegion]["aq"] = []
                if not "aq" in dontreject[generegion]: 
                    dontreject[generegion]["aq"] = []
                if aq(aln_nc_slop, resgene, ogene, compare, chop):
                    reject[generegion]["aq"] += rjuncs
                else:
                    dontreject[generegion]["aq"] += rjuncs
                    ###############################
                    # Minimal gains check.
                    ###############################
                # Minimal gains: nonconstant portion must exhibit an improvement of >5
                # in one of the triplet scores.
                if not "minimalgains" in reject[generegion]: 
                    reject[generegion]["minimalgains"] = []
                if not "minimalgains" in dontreject[generegion]: 
                    dontreject[generegion]["minimalgains"] = []
                chopo = [a for a in chop if a.id == ogene][0]
                chopr = [a for a in chop if a.id == resgene][0]
                nonconstantregion = [a for a in range(len(chopo)) if chopo[a] != chopr[a]]
                aln_nc = chop_alignment(chop, nonconstantregion)
                if minimal_gains(aln_nc, resgene, ogene, compare):
                    reject[generegion]["minimalgains"] += rjuncs
                else:
                    dontreject[generegion]["minimalgains"] += rjuncs
                ################################
                # Parallel junctions check
                ################################
                # Check for regions where a triplet was blank before, but
                # has become at least partially not blank, including GOI.
                if not "parallels" in reject[generegion]:
                    reject[generegion]["parallels"] = []
                if not "parallels" in dontreject[generegion]:
                    dontreject[generegion]["parallels"] = []
                if parallels(aln_nc, resgene, ogene, compare):
                    reject[generegion]["parallels"] += rjuncs
                else:
                    dontreject[generegion]["parallels"] += rjuncs

    for generegion in reject:
        for test in reject[generegion]:
            reject[generegion][test] = list(
                set([a for a in reject[generegion][test] if not a in dontreject[generegion][test]]))
    return reject


def gr_pick(agene, ogene):
    return agene == ogene or not get_gr(agene) == get_gr(ogene)


def get_introns(pdict):
    res = {}
    keys = sorted(pdict.keys())
    for i in range(len(pdict) - 1):
        gtf = pdict[i]["gtfline"][:]
        gtf[3] = int(pdict[i]["gtfline"][4]) + 1
        gtf[4] = int(pdict[i + 1]["gtfline"][3]) - 1
        res[i] = {"gtfline": gtf, "ids": [keys[i], keys[i + 1]]}
    return res


def grab_interval_probe(leftright, rightleft, intervals, gene, o, flavour="", direction=""):
    res = {"interval": intervals, "gene": gene, "ogene": o}
    if flavour: 
        res["flavour"] = flavour
    if leftright: 
        res["l"] = max(leftright)
    if rightleft: 
        res["r"] = min(rightleft)
    res["direction"] = direction
    return res


def grab_aln_intervals_in(o_interval, r_interval, gene, o, c_parts, c_coords, path_fil, alnlen, force=False):
    i_locs = []
    if not o_interval == r_interval or force:
        cpg = c_parts[gene].keys()
        opg = c_parts[o].keys()
        r_part_ids = \
            [[cpg[i], cpg[i + 1]] for i in cpg[:-1] if int(c_parts[gene][cpg[i]]["gtfline"][4]) == r_interval[0] - 1 \
             and int(c_parts[gene][cpg[i + 1]]["gtfline"][3]) == r_interval[1] + 1][0]
        ri_loc = c_coords[gene][r_part_ids[0]][-15:-1] + c_coords[gene][r_part_ids[1]][0:15]
        if o_interval:
            o_part_lid = [i for i in opg[:-1] if int(c_parts[o][i]["gtfline"][4]) == o_interval[0] - 1][0]
            o_part_rid = [i for i in opg[1:] if int(c_parts[o][i]["gtfline"][3]) == o_interval[1] + 1][0]
            oi_loc = c_coords[o][o_part_lid][-15:-1] + c_coords[o][o_part_rid][0:15]
            i_locs = interpolate(ri_loc + oi_loc)
        else:
            i_locs = interpolate(ri_loc)
        i_locs = cds_to_aa_locs(i_locs)
        filstring = "".join(["Y" if i in i_locs else "-" for i in range(alnlen)])
        call_function("echo \">" + gene + "." + o + "\n" + filstring + "\">> " + path_fil)
    return i_locs


def grab_aln_intervals_fl(o_interval, r_interval, gene, o, c_parts, c_coords, path_fil, alnlen, dirn, force=False):
    i_locs = []
    ik1, ik2 = [0, 15] if dirn == "r" else [-15, -1]
    if not o_interval == r_interval or force:
        cpg = c_parts[gene].keys()
        opg = c_parts[o].keys()
        if dirn == "l":
            r_part_id = [i for i in cpg if int(c_parts[gene][i]["gtfline"][4]) == r_interval[0] - 1][0]
            o_part_id = [i for i in opg if int(c_parts[o][i]["gtfline"][4]) == o_interval[0] - 1][0]
            ri_loc = c_coords[gene][r_part_id][ik1:ik2]
            oi_loc = c_coords[o][o_part_id][ik1:ik2]
            i_locs = interpolate(ri_loc + oi_loc)
        else:
            r_part_id = [i for i in cpg if int(c_parts[gene][i]["gtfline"][3]) == r_interval[-1] + 1][0]
            o_part_id = [i for i in opg if int(c_parts[o][i]["gtfline"][3]) == o_interval[-1] + 1][0]
            ri_loc = c_coords[gene][r_part_id][ik1:ik2]
            oi_loc = c_coords[o][o_part_id][ik1:ik2]
            i_locs = interpolate(ri_loc + oi_loc)
        i_locs = cds_to_aa_locs(i_locs)
        filstring = "".join(["T" if i in i_locs else "-" for i in range(alnlen)])
        call_function("echo \">" + gene + "." + o + "\n" + filstring + "\">> " + path_fil)
    return i_locs


def grab_aln_intervals_ex(o_interval, r_interval, gene, o, c_parts, c_coords, path_fil, alnlen, direction, force=False):
    switch_key = 1 if direction == "l" else 0
    gtf_key = 3 if direction == "l" else 4
    ik1, ik2 = [0, 15] if direction == "l" else [-15, -1]
    cipher = "D" if direction == "l" else "R"
    i_locs = []
    if not o_interval == r_interval or force:
        rpart = [a for a in c_parts[gene] if int(c_parts[gene][a]["gtfline"][gtf_key]) == r_interval[switch_key]][0]
        ri_loc = c_coords[gene][rpart][ik1:ik2]
        oi_loc = []
        if o_interval:
            opart = [a for a in c_parts[o] if int(c_parts[o][a]["gtfline"][gtf_key]) == o_interval[switch_key]][0]
            oi_loc = c_coords[o][opart][ik1:ik2]
            i_locs = interpolate(oi_loc + ri_loc)
        else:
            i_locs = ri_loc
        i_locs = cds_to_aa_locs(i_locs)
        filstring = "".join([cipher if i in i_locs else "-" for i in range(alnlen)])
        call_function("echo \">" + gene + "." + o + "\n" + filstring + "\">> " + path_fil)
    return i_locs


########################################
# All the other functions
########################################

def prepare_winner_alignment(path_w_dir, winners, dict_generegions, or_fa=""):
    path_winners_out = blank_file(path_w_dir + "/all.options.winners")  # qe
    path_winners_aln = path_winners_out + ".aln"  # qe
    for i in winners:
        generegion, option = parse_id(i)
        path_fasta = dict_generegions[generegion]["options"][int(option)]["fasta"]
        call_function("cat " + path_fasta + " | sed -r \"s/\*/X/g\" >> " + path_winners_out)
    if or_fa:
        call_function("cat " + or_fa + " >> " + path_winners_out)
    align(path_winners_out, path_winners_aln)
    return path_winners_aln


def bound_align_manual(parts_info, cds, path_stub, slookup_rev):
    """Align a set of boundaries from the same region.
    """
    all_gtf = flatten_lol([[q["gtfline"] for q in parts_info[p]["parts"]] for p in parts_info])
    write_csv(all_gtf, path_stub + ".all.gtf")
    merge_friendly(path_stub + ".all.gtf", path_stub + ".regions.gtf", sort=True)
    regions = dict(enumerate(read_csv(path_stub + ".regions.gtf")))
    # Get a cds alignment, to convert later to an aa alignment.
    # Do this separately for each part in the merge, so that non-compatible parts
    # are not placed on top of each other.
    regionsblank = dict((r, "-" * (int(regions[r][4]) - int(regions[r][3]) + 1)) for r in regions)
    partsout = []
    partstrings = {}
    aaout = []
    for p in parts_info:
        # Each part will be contained in exactly one of the regions.
        regionlines = copy.deepcopy(regionsblank)
        regionsstring = ""
        for part in parts_info[p]["parts"]:
            regionposs = [r for r in regions if overlap_in_frame(regions[r], part["gtfline"])]
            if regionposs:
                region = regionposs[0]
                line = "-" * (int(part["gtfline"][3]) - int(regions[region][3])) \
                       + cds[int(part["gtfline"][3]) - 1: int(part["gtfline"][4])] \
                       + "-" * (-int(part["gtfline"][4]) + int(regions[region][4]))
                regionlines[region] = add_region_line(regionlines[region], line)
        for region in sorted(regions.keys()):
            regionsstring += regionlines[region]
        partstrings[p] = regionsstring
        partsout += [[">" + slookup_rev[p]], [regionsstring]]
    # The method assumes that each gtf is properly written (and starts in-frame), as it should be.
    # Might need to keep an eye on this though...
    aastrings = {}
    for p in partstrings:
        aalength = len(partstrings) / 3
        aa_raw = parts_info[p]["aa"]
        aastring = ""
        nucpos = 0
        for i, e in enumerate(partstrings[p]):
            if i % 3 == 0:
                if e == "-":
                    aastring += "-"
                else:
                    # Ends aren't always full.
                    if nucpos / 3 < len(aa_raw):
                        aastring += aa_raw[nucpos / 3]
            if e != "-":
             nucpos = nucpos + 1
        aastrings[p] = aastring
    aalength = max([len(aastrings[aa]) for aa in aastrings])
    # There can sometimes be slight length imbalances, due to hanging codons. Fill them in with blanks.
    for aa in aastrings:
        aastring = aastrings[aa]
        if not len(aastring) == aalength:
            aastring = aastring + "-" * (aalength - len(aastring))
        aaout += [[">" + slookup_rev[aa]], [aastring]]
    write_csv(partsout, path_stub + ".regions")
    write_csv(aaout, path_stub + ".regions.aa")
    return path_stub + ".regions.aa"


def add_region_line(existing, new):
    # Add a region over the existing one, such that adding "-" counts for nothing.
    if len(existing) != len(new): 
        raise Exception()
    return string.join([existing[i] if new[i] == "-" else new[i] for i in range(len(existing))].join(""))


def prepare_choose_winners(dict_options, path_w_dir):
    path_all_fasta = path_w_dir + "/all.options.fasta"
    path_all_aln = path_all_fasta + ".aln"
    # Sort out the sequences for processing
    seqlookup = {}
    seqlookup_rev = {}
    towrite = {}
    opts_copy = copy.deepcopy(dict_options)
    for generegion in opts_copy:
        seqlookup[generegion] = {}
        towrite[generegion] = []
        seqlookup_rev[generegion] = {}
        for option in opts_copy[generegion]:
            seqid = option.id
            shortid = generegion + ".option_" + str(len(seqlookup[generegion]))
            seqlookup[generegion][shortid] = seqid
            seqlookup_rev[generegion][seqid] = shortid
            option.id = shortid
            towrite[generegion].append(option)
    return path_all_fasta, path_all_aln, seqlookup, seqlookup_rev, towrite


def choose_winners_vanilla(dict_options, path_w_dir, try_blank={}):
    path_all_fasta, path_all_aln, seqlookup, seqlookup_rev, towrite = prepare_choose_winners(dict_options, path_w_dir)
    if not try_blank and all(len(a) <= 1 for a in dict_options.values()):
        return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
    with open(path_all_fasta, "w") as f:
        for generegion in towrite: write_seqs(towrite[generegion], f)
    align(path_all_fasta, path_all_aln, safety=False)
    return choose_them(path_all_aln, try_blank, seqlookup, seqlookup_rev, path_w_dir)


def choose_winners_ref(dict_options, path_w_dir, try_blank={}, path_ref_aln="", double_check=False, spanning_band=False,
                       or_aln=False, parallel_check=False):
    path_all_fasta, path_all_aln, seqlookup, seqlookup_rev, towrite = prepare_choose_winners(dict_options, path_w_dir)
    if not try_blank and all(len(a) <= 1 for a in dict_options.values()):
        return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
    # initial microexons can cause zero-length aa strings to be formed, which can't be aligned.
    nonempty = [a for a in flatten_lol(towrite.values()) if len(a) != 0]
    empty = [a for a in flatten_lol(towrite.values()) if len(a) == 0]
    write_seqs(nonempty, path_all_fasta)
    # Store the realigned reference so that we can compare the original.
    path_ref_out = path_all_aln + ".ref"
    align_ref(path_all_fasta, path_all_aln, path_ref_aln, path_ref_out)
    # Add blank sequences to the alignment file.
    a = read_seqs(path_all_aln)
    write_seqs(a + [make_seq(t.id, "-" * len(a[0])) for t in empty], path_all_aln)
    # If double_check is on, we extract the nonconstant portion from the alignment and
    # align that using l-ins-i (it'll be a relatively narrow column so won't eat up
    # too much compute power). We then compare the alignment score for the two and
    # go with which ever one is better.
    if double_check: 
        path_all_aln, path_ref_out = realign_band(path_all_aln, path_ref_out, spanning_band=spanning_band)
    # If this option is turned on, check for "parallel junction expansion": that is,
    # events that force gaps to be created in the original alignment. Simply reject any such events.
    #	if or_aln and parallel_check: path_all_aln = reject_parallels(path_all_aln, path_ref_out)
    return choose_them(path_all_aln, try_blank, seqlookup, seqlookup_rev, path_w_dir)


def reject_parallels(path_all_aln, path_ref_aln):
    # Reject anything which has something when the original alignment has all gaps.
    # The ref and all alignments should be the same length.
    res_seqs = read_seqs(path_all_aln)
    ref_seqs = [a for a in read_seqs(path_ref_aln) if "orig" in a.id]
    sbm = group_sequences_by_generegion(res_seqs)
    path_noparallels = path_all_aln + ".nopara"
    to_keep = get_nc(sbm)
    if not res_seqs or not ref_seqs:
        write_seqs([], path_noparallels)
        return path_noparallels
    if not len(res_seqs[0]) == len(ref_seqs[0]): 
        raise
    reject = []
    res_seqs = chop_alignment(res_seqs, to_keep)
    ref_seqs = chop_alignment(ref_seqs, to_keep)
    for i in range(len(res_seqs[0])):
        if all(a[i] == "-" for a in ref_seqs):
            for a in res_seqs:
                if not a[i] == "-": 
                    reject.append(a.id)
    write_seqs([a for a in res_seqs if not a.id in reject], path_noparallels)
    write_seqs(res_seqs + ref_seqs, path_noparallels + ".check")
    return path_noparallels


def realign_band(path_in_aln, path_ref_aln, thebuffer=20, spanning_band=False):
    # Find the window in the alignment where the interesting stuff is happening.
    # Then realign this with l-ins-i, and double check the results. Since the
    # aligner may have made a mistake (particularly if the alignment is messy),
    # there might be more than one such region.
    path_column = path_in_aln + ".column.fasta"
    path_columnref = path_in_aln + ".column.fasta.ref"
    # Don't go any further if there are no input sequences
    inseqs = read_seqs(path_in_aln)
    if not inseqs: 
        return path_in_aln, path_ref_aln
    # Remove any all-blanks...
    seqlen = len(inseqs[0])
    removecols = [i for i in range(seqlen) if all(s[i] == "-" for s in inseqs)]
    inseqs_nonblank = chop_alignment(inseqs, removecols, True)
    refseqs_nonblank = chop_alignment(read_seqs(path_ref_aln), removecols, True)
    # Get to work. Want to chop out the windows where everything is fine: these
    # are not interesting to us.
    sbm = group_sequences_by_generegion(inseqs_nonblank)
    if len(sbm.keys()) < 2: 
        return path_in_aln, path_ref_aln
    # The squeezemap is a binary representation of whether a column should be considered.
    # Grab the coordinates of the nonuniform bits, with a buffer added.
    to_keep = get_nc(sbm, thebuffer)
    # There's a choice between taking the span of all of the keep regions and
    # Considering them separately. The spanning region is more accurate but
    # Could be costly if the regions are distant
    if not to_keep: 
        return path_in_aln, path_ref_aln
    to_keep = range(min(to_keep), max(to_keep) + 1) if spanning_band else sorted(set(to_keep))
    sbmalt = dict((k, chop_alignment(sbm[k], to_keep)) for k in sbm)
    refseqs_keep = chop_alignment(refseqs_nonblank, to_keep)
    # Write out the newly found columns.
    for s in refseqs_keep: s.id = "reference." + s.id
    write_seqs(flatten_lol(sbmalt.values()), path_column)
    write_seqs(refseqs_keep, path_columnref)
    # Retain the lateral context of the band that we've taken out.
    # This involves just adding in the bits we ignored.
    squeeze_fa = grab_squeeze_map(inseqs, removecols, to_keep, sbm)
    # Align and compare with the original alignment.
    path_column_all = path_column + ".all"
    path_column_aln = path_column + ".aln"
    write_seqs(flatten_lol(sbmalt.values()) + refseqs_keep, path_column_all)
    write_seqs(squeeze_fa, path_column + ".squeeze.fa")
    # Align the bastard.
    align(path_column_all, path_column_aln)
    oldscore = alignment_score([a for a in read_seqs(path_column) if not "ref" in a.id])
    newscore = alignment_score([a for a in read_seqs(path_column_aln) if not "ref" in a.id])
    # If the realigned sequences are better, use those. If not use the original ref alignment.
    path_sequences = path_column if oldscore > newscore else path_column_aln
    write_seqs([a for a in read_seqs(path_sequences) if not "ref" in a.id], path_column + ".bestaln")
    write_seqs([a for a in read_seqs(path_sequences) if "ref" in a.id], path_column + ".bestaln.ref")
    return path_column + ".bestaln", path_column + ".bestaln.ref"


def get_nc(sbm, thebuffer=0):
    seqlen = len(sbm[sbm.keys()[0]][0])
    squeeze_map = [0 if all(len(set([a[i] for a in sbm[k]])) == 1 for k in sbm) else 1 for i in range(seqlen)]
    to_keep = flatten_lol(range(i - thebuffer, i + thebuffer + 1) for (i, e) in enumerate(squeeze_map) if e == 1)
    return [i for i in to_keep if i >= 0 and i < seqlen]


def grab_squeeze_map(oseqs, removecols, to_keep, sbm):
    squeeze_fa = []
    sbmo = group_sequences_by_generegion(oseqs)
    for k in sbm:
        if not sbmo[k]: 
            continue
        # The sequences by definition are identical when chopped, so just take the first.
        s = sbmo[k][0]
        s_adj = [s[i] for i in range(len(s)) if not i in removecols]
        s_adjt = [s[i] for i in range(len(s)) if i in removecols]
        ts = "".join([("@" if i in to_keep else s_adj[i]) for i in range(len(s_adj))] + s_adjt)
        squeeze_fa.append(replace_seq(s, Seq(ts)))
    return squeeze_fa


def choose_winners_seeded(dict_options, path_w_dir, seed="", seed_info={}, cds_bases={}, try_blank={}):
    path_all_fasta, path_all_aln, seqlookup, seqlookup_rev, towrite = prepare_choose_winners(dict_options, path_w_dir)
    if not try_blank and all(len(a) <= 1 for a in dict_options.values()):
        return dict((k, dict_options[k][0]) for k in dict_options if dict_options[k])
    if seed == "auto":
        seeds = []
        for generegion in towrite:
            path_d_fa = path_all_fasta + "." + generegion + ".fa"
            write_seqs(towrite[generegion], path_d_fa)
            seeds.append(path_d_fa)
        align_seeded(seeds, path_all_aln, prealigned=False, safety=False)
    elif seed == "manual":
        # Manually align the seed sequences. This will only work for sequences from
        # the same gene region.
        seeds = []
        for generegion in seed_info:
            if seed_info[generegion]:
                seeds.append(
                    bound_align_manual(seed_info[generegion], cds_bases[generegion], path_all_fasta + "." + generegion,
                                       seqlookup_rev[generegion]))
        align_seeded(seeds, path_all_aln, prealigned=True, safety=True)
    return choose_them(path_all_aln, try_blank, seqlookup, seqlookup_rev, path_w_dir)


def choose_them(path_all_aln, try_blank, seqlookup, seqlookup_rev, path_w_dir):
    # Gather the sequences and pick the best set.
    sequences = [a for a in read_seqs(path_all_aln) if not (a.id).startswith("prev.")]
    if not sequences:
        return {}

    # The propercounts option ensures we fill in absent options with blanks.
    seqlen = len(sequences[0].seq)
    for k in [a for a in try_blank if try_blank[a]]:
        sequences.append(make_seq(k + ".blank", "-" * seqlen))
        if not k in seqlookup:
            seqlookup[k] = {}
            seqlookup_rev[k] = {}
        seqlookup[k][k + ".blank"] = k + ".blank"
        seqlookup_rev[k][k + ".blank"] = k + ".blank"
    write_seqs(sequences, path_w_dir + "/options.aln")
    sbm = group_sequences_by_generegion(sequences)
    if len(sbm.keys()) == 1:
        return {sbm.keys()[0]: make_seq(sbm.keys()[0] + ".blank", "")}
    winner = most_coherent_set(sequences, path_w_dir)
    res = {}
    for k, w in winner.items():
        w.id = seqlookup[k][winner[k].id]
        res[k] = w
    return res


def get_winner_aln_coords(winners, path_winners_aln, dict_generegions, path_out):
    aligned = [a for a in read_seqs(path_winners_aln) if not "orig" in a.id]
    orders = {}
    gtfs = {}
    for w in winners:
        generegion, option = parse_id(w)
        winnerpick = dict_generegions[generegion]["options"][int(option)]
        dict_generegions[generegion]["winner"] = winnerpick
        wparts = winnerpick["parts"]
        part_ids = sorted(wparts.keys())
        dict_generegions[generegion]["winner"] = wparts
        gtfd = dict((a, wparts[a]["gtfline"]) for a in part_ids)
        # By design this gtf should be clean and safe.
        order = gtfd.keys()
        orders[w] = order
        gtfs[w] = [gtfd[o] for o in order]
        for key in part_ids:
            parth = wparts[key]
            parth["initial"] = (int(key) == 0)
            parth["terminal"] = (int(key) == len(part_ids) - 1)
            parth["id"] = int(key)
    locs_aa, locs, good_gtfs = flatten_alignment(aligned, gtfs, pre_sorted=True, path_out=path_out)
    dict_part_coords = {}
    for w in winners:
        generegion, option = parse_id(w)
        dict_part_coords[generegion] = {}
        for i, o in enumerate(orders[w]):
            l = locs[w][i]
            dict_part_coords[generegion][o] = [min(l), max(l)]
    return dict_part_coords


def prepare_options(dict_generegions):
    options = idict(dict_generegions.keys(), [])
    for generegion in dict_generegions:
        dict_options = dict_generegions[generegion]["options"]
        options[generegion] = flatten_lol([read_seqs(option["fasta"]) for option in dict_options.values()])
    winners = choose_winners_vanilla(options, path_w_dir).values()
    return [a.id for a in winners]


def prepare_compare(res, dict_generegions, dict_seq_info, path_dbl):
    # Extract all the information from a res necessary to compare
    # Its parts with the original. We're going to do the comparison
    # A couple of times so it's good to have this separate.
    comparator_gtf = {}
    comparator_aa = {}
    comparator_parts = {}
    comparator_generegions = {}
    comparator_sequences = {}
    comparator_mcount = {}
    # Work through all the gene regions and add in the genes.
    for generegion in dict_generegions:
        cds = dict_generegions[generegion]["cds_base_data"]
        comparator_mcount[generegion] = 1
        # Do all the original sequences first.
        for sequence in dict_generegions[generegion]["sequences"]:
            seqid = generegion + "-original-" + sequence
            # Our result shouldn't contain non-viable genes. HOWEVER, the original input
            # could well do. We must protect ourselves against this.
            # In the case that there is an invalid gene input to begin with and no other result
            # is acquired, output for said gene will simply be blank.
            local_cds = dict_seq_info[sequence]["sequence_cds"]
            aa = dict_seq_info[sequence]["sequence_aa"]
            if contains_stop(local_cds[:-3]): continue
            comparator_gtf[seqid] = read_csv(dict_seq_info[sequence]["gtf_rel"])
            # Make sure it ends properly!
            comparator_aa[seqid] = re.sub(r"([A-Z])$", r"\1*", aa)
            comparator_generegions[seqid] = generegion
            comparator_sequences[seqid] = sequence
            comparator_mcount[generegion] += 1
        # Finally add in the result. Careful though, it might not exist!
        if generegion in res:
            newid = generegion + "-result"
            comparator_gtf[newid] = [p["gtfline"] for p in res[generegion]]
            comparator_aa[newid] = translate_gtf_live(comparator_gtf[newid], cds)
            comparator_generegions[newid] = generegion
            comparator_sequences[newid] = "result"
    # Grab the parts
    comparator_parts = grab_comparator_parts(comparator_gtf, comparator_generegions, comparator_aa, path_dbl)
    sequences = [make_seq(c, comparator_aa[c]) for c in comparator_aa]
    # Make an alignment
    path_sequences = path_dbl + "/comparator.fa"
    path_aln = path_dbl + "/comparator.aln"
    write_seqs(sequences, path_sequences)
    align(path_sequences, path_aln)
    # Grab the local coords for the alignment.
    alignedseqs = dict((s.id, s) for s in read_seqs(path_aln))
    locs_aa, locs, good_gtfs = flatten_alignment(alignedseqs.values(), comparator_gtf, path_aln + ".flat")
    return locs, comparator_parts, comparator_mcount, path_aln


def revive_casualties(res, d_gr, d_si, path_w_dir):
    # Designed to revive any genes that have been killed.
    casualties = [k for k in d_gr if not k in res]
    safe = [k for k in d_gr if k in res]
    if casualties:
        options = dict((k, [res[k]]) for k in safe)
        for k in casualties: options[k] = [make_mock_part_sets(d_si[s]) for s in d_gr[k]["sequences"]]
        return scrape_bucket(options, path_w_dir, d_gr)
    else:
        return res


def scrape_bucket(options, path_w_dir, d_gr):
    # If there is only one sequences per generegion, return this list. Otherwise choose between.
    if all(len(s) == 1 for s in options.values()):
        return dict((k, options[k][0]) for k in options)
    else:
        cas_dir = make_if_absent(path_w_dir + "/casualties")
        labels = idict(options, {})
        aas = idict(options, [])
        for k in options:
            for l, s in enumerate(options[k]):
                labels[k][k + "." + str(l)] = s
                aas[k].append(make_seq(k + "." + str(l), translate_part(
                    get_cds_live([a["gtfline"] for a in s], d_gr[k]["cds_base_data"]))))
        winners = choose_winners_vanilla(aas, cas_dir)
        resn = {}
        for k in winners:
            resn[k] = labels[k][winners[k].id]
        return resn


def check_aq_improvement(res_o, d_gr, d_si, path_w_dir, or_score):
    seqs = {}
    aas = [make_seq(k + "." + str(random.random()),
                    translate_part(get_cds_live([a["gtfline"] for a in res_o[k]], d_gr[k]["cds_base_data"]))) for k in
           res_o]
    write_seqs(aas, path_w_dir + "/res.check.fa")
    align(path_w_dir + "/res.check.fa", path_w_dir + "/res.check.aln")
    alnseqs = read_seqs(path_w_dir + "/res.check.aln")
    if alignment_score(alnseqs) < or_score:
        for k in d_gr:
            seqs[k] = [make_mock_part_sets(d_si[s]) for s in d_gr[k]["sequences"]]
        return scrape_bucket(seqs, path_w_dir, d_gr)
    else:
        return res_o


def make_mock_part_sets(s):
    # All we really need is the gtflines.
    gtf = read_csv(s["gtf_rel"])
    gtf = safe_gtf(gtf)
    partslist = [{"gtfline": line} for line in gtf]
    return partslist


def compare_with_original(res, dict_generegions, dict_seq_info, path_w_dir):
    # Align the original and the new sequences.
    # Group them into overlap regions.
    # Choose winners for each as usual.
    sprint("Double-checking against original sequences...")
    path_dbl = make_if_absent(path_w_dir + "/doublecheck")
    dict_cds_bases = dict((k, dict_generegions[k]["cds_base_data"]) for k in dict_generegions)
    # Prepare containers for all the information we're about to gather.
    # Grab the dna sequences for each generegion.
    c_coords, c_parts, c_mcount, path_aln = prepare_compare(res, dict_generegions, dict_seq_info, path_dbl)
    groups, adj = group_parts(c_coords, ranges=True)
    # Do a first round of refinements by making a nice patchwork of old and new bits
    # based on what are the best choices.
    res_refined = compare_parts(adj, c_parts, c_mcount, path_dbl, dict_cds_bases, path_aln, minintron, minexon,
                                dict_generegions)
    # The second round of refinements applies various filters to make sure nothing
    # too exciting is going on. (Aq, Mg, Pl filters.)
    sprint("Applying filters...")
    path_fltrs = make_if_absent(path_w_dir + "/filters")
    c_coords, c_parts, c_mcount, path_aln = prepare_compare(res_refined, dict_generegions, dict_seq_info, path_fltrs)
    groups, adj = group_parts(c_coords, ranges=True)
    res_filtered = filter_changes(adj, c_parts, c_coords, path_fltrs, dict_cds_bases, path_aln, minintron, minexon,
                                  dict_generegions, c_mcount)
    # We then do one more run of compare parts to double check that the filtered parts
    # work together to form viable genes. Shouldn't be too much to check here.
    return res_filtered


def grab_comparator_parts(comparator_gtf, comparator_generegions, comparator_aa, path_dbl):
    comparator_parts = {}
    for newid in comparator_gtf:
        comparator_parts[newid] = {}
        full_aa = comparator_aa[newid]
        progress = 0
        gtflength = len(comparator_gtf[newid])
        for i, line in enumerate(sorted(comparator_gtf[newid], key=lambda x: int(x[3]))):
            part = {}
            cdslength = int(line[4]) - int(line[3]) - int(line[7]) + 1
            partlength = int(math.ceil(cdslength / 3.0))
            part["gtfline"] = line
            part["part"] = full_aa[progress: progress + partlength]
            part["part_id"] = newid + "." + str(i)
            part["initial"] = (i == 0)
            part["terminal"] = (i == gtflength)
            part["partlength"] = partlength
            part["generegion"] = comparator_generegions[newid]
            progress += partlength
            comparator_parts[newid][i] = part
    sequences = [make_seq(c, comparator_aa[c]) for c in comparator_aa]
    path_sequences = path_dbl + "/comparator.fa"
    return comparator_parts


def initialise_end_options(part):
    part["options_left"] = idict([0, 1, 2], [])
    part["options_right"] = idict([0, 1, 2], [])
    part["left_flat"] = []
    part["right_flat"] = []
    part["initial_sites"] = []
    part["terminal_sites"] = []
    part["donor_sites"] = []
    part["acceptor_sites"] = []
    return part


def initialise_probe(parts, probe_id):
    probe = copy.deepcopy(parts[0])
    probe["gtfline"][3] = min([p["gtfline"][3] for p in parts])
    probe["gtfline"][4] = max([p["gtfline"][4] for p in parts])
    probe["status"] = "fresh"
    probe["evidence"] = 0
    probe["id"] = probe_id
    probe = initialise_end_options(probe)
    return probe


def add_part_to_probe(probe, part, cds):
    gtf = part["gtfline"]
    frame = int(gtf[7])
    endframe = get_end_frame(gtf)
    probe["ogtfline"] = part["gtfline"][:]
    probe["options_left"][frame] = l = list(set(probe["options_left"][frame] + [int(gtf[3])]))
    probe["options_right"][endframe] = r = list(set(probe["options_right"][endframe] + [int(gtf[4])]))
    probe["evidence"] += 1
    probe["constituent_ids"] = probe.get("constituent_ids", []) + [part["part_id"]]
    # Find out which bits are terminal, initial, donor, acceptor.
    if is_initial_part(gtf, cds): 
         probe["initial_sites"] += [int(gtf[3])]
    if is_terminal_part(gtf, cds): 
        probe["terminal_sites"] += [int(gtf[4])]
    if is_acceptor_part(gtf, cds): 
        probe["acceptor_sites"] += [int(gtf[3])]
    if is_donor_part(gtf, cds): 
           probe["donor_sites"] += [int(gtf[4])]
    # Add the flat options
    probe["left_flat"] = list(set(l + probe["left_flat"]))
    probe["right_flat"] = list(set(r + probe["right_flat"]))
    # Allow the options to be restored later if necessary
    probe["lbackupflat"] = copy.deepcopy(probe["left_flat"])
    probe["rbackupflat"] = copy.deepcopy(probe["right_flat"])
    probe["lbackup"] = copy.deepcopy(probe["options_left"])
    probe["rbackup"] = copy.deepcopy(probe["options_right"])
    return probe


def create_part(q, status, startframe=-1, pid=""):
    p = copy.deepcopy(q)
    gtf = p["gtfline"]
    sf = startframe if startframe != -1 else int(gtf[7])
    p["status"] = status
    # Back up the options from before
    p["rbackupflat"] = p["right_flat"]
    p["lbackupflat"] = p["left_flat"]
    p["rbackup"] = p["options_right"]
    p["lbackup"] = p["options_left"]
    # Replace the actual options with the new data.
    p["options_left"] = idict([0, 1, 2], [])
    p["options_right"] = idict([0, 1, 2], [])
    p["left_flat"] = [int(gtf[3])]
    p["right_flat"] = [int(gtf[4])]
    p["options_left"][sf] = [int(gtf[3])]
    p["options_right"][get_end_frame(gtf)] = [int(gtf[4])]
    if pid: 
        p["id"] = pid
    return p


def prepare_probe(parts, i, cds):
    # Now we add in the options for each part based on what we have here.
    probe = initialise_probe(parts, i)
    for p in parts: probe = add_part_to_probe(probe, p, cds)
    # If there isn't sufficient representation at this point, we would like to
    # probe not having the part there at all.
    return probe


def compare_parts_heart(i, a, dict_parts, dict_mcount, path_ref_aln, d_reject, d_gr, res, tryempty, prev_probes, mi, mx,
                        path_chk_dir):
    # Create a set of peobes. Make sure that the only parts added in are parts that haven't already been added in.
    # If there are no parts available for a gene region, assume that neither new nor old has anything there.
    probes = {}
    for k in d_gr:
        parts = flatten_lol(
            [[copy.deepcopy(p) for p in q.values() if p["part_id"] in a and p["generegion"] == k] for q in
             dict_parts.values()])
        # Make sure we haven't already added in these parts!
        if not k in res: 
            res[k] = []
        if res[k]: 
            parts = [copy.deepcopy(part) for part in parts if not any(
            ["constituent_ids" in q and part["part_id"] in q["constituent_ids"] for q in res[k]])]
        # If there are no parts, assume both the result and the original had nothing to offer.
        if not parts: 
            continue
        probes[k] = prepare_probe(parts, i, d_gr[k]["cds_base_data"])
        tryempty[k] = probes[k]["evidence"] < dict_mcount[k]
    # Get the part combinations for option probing.
    # Need to be a little careful.
    dict_opts = {}
    labels = {}
    write(d_reject, path_chk_dir + "/reject.tpy")
    write(probes, path_chk_dir + "/probes.tpy")
    write(prev_probes, path_chk_dir + "/prev_probes.tpy")
    write(res, path_chk_dir + "/plist.tpy")
    for k in res:
        plist = copy.deepcopy(res[k]) if res[k] else []
        latest = [probes[k]] if k in probes else []
        if not plist and not latest: 
            continue
        labels[k] = careful_part_combinations(plist, latest, d_gr[k]["cds_base_data"], k, mi, mx,
                                              tryempty[k] or not latest, reject=d_reject.get(k, []))
    # Also try adding back in any probes from last time that were omitted.
    # Aims to mitigate mistakes.
    for k in res:
        all_ids = [a["part_id"] for a in res[k]]
        if k in prev_probes and k in res and not any([a in all_ids for a in prev_probes[k]["constituent_ids"]]):
            plist = copy.deepcopy(res[k]) if res[k] else []
            latest = [prev_probes[k]] + ([] if not k in probes else [probes[k]])
            labels_l = careful_part_combinations(plist, latest, d_gr[k]["cds_base_data"], k, mi, mx, tryempty[k],
                                                 reject=d_reject.get(k, []), dr=path_chk_dir)
            for l in labels_l: labels[k][l] = labels_l[l]
    dict_opts = dict((k, [l["seq"] for l in labels[k].values()]) for k in labels)
    return labels, dict_opts, probes, tryempty


def compare_parts(adj, dict_parts, dict_mcount, path_wdir, dict_cds_bases, path_ref_aln, mi, mx, d_gr, d_reject={},
                  refine=True):
    # Add in the parts as usual. At each stage pick which part works best in the context.
    # At each stage, if there are not k+1 options, where k is the total number of original
    # sequences for the generegion,  we assume removal of the part to be one of the options.
    # of the part. Note two identical options is still two options. Note also that we need to
    # pay attention to whether or not a part is initial, terminal, whatever.
    res = idict(d_gr, {})
    tryempty = idict(d_gr, False)
    probes = []
    # Now run through all the different part options in their groups, slot them in and
    # sequentially pick the best of them.
    for i, a in enumerate(adj):
        path_chk_dir = make_if_absent(path_wdir + "/check_" + str(i))
        labels, options, probes, tryempty = compare_parts_heart(i, a, dict_parts, dict_mcount, path_ref_aln,
                                                                d_reject, d_gr, res, tryempty, probes, mi, mx,
                                                                path_chk_dir)
        write(labels, path_chk_dir + "/labels.tpy")
        # Now align these lovely options and pick the best one.
        res, path_aln = process_labels(labels, options, path_chk_dir, path_ref_aln, True, refine)
        seal_done(res)
    res = check_for_starts_and_stops(res, dict_cds_bases, dict_parts, dict_mcount.keys(), path_wdir, mi, mx, d_reject)
    return res


def seal_done(res):
    for plist in res.values():
        for part in plist:
            if status_check(part["status"], ["done", "waiting"]):
                part["left_flat"] = [int(part["gtfline"][3])]
            if status_check(part["status"], ["done"]):
                part["right_flat"] = [int(part["gtfline"][4])]


def get_score(prev_aln):
    # The lengths of the bands here might not be the same size. We correct for this
    # By adding in the extra bits of sequenc from the original alignment.
    # The order in which things are added back in doesn't matter, since the alignment
    # score is columnwise.
    squeeze_map = re.sub(r"bestaln", r"squeeze.fa", prev_aln)
    if os.path.exists(squeeze_map):
        # There should be one for each gene region.
        filled = dict((get_gr(s.id), s) for s in read_seqs(prev_aln))
        for s in read_seqs(squeeze_map):
            generegion = get_gr(s.id)
            if generegion in filled:
                filled[generegion].seq += Seq(str(s.seq).replace("@", ""))
        write_seqs(filled.values(), os.path.dirname(prev_aln) + "/all.filled.fa")
        return alignment_score(filled.values(), scaled=False)
    else:
        return alignment_score(read_seqs(prev_aln), scaled=False)


def process_labels(labels, options, path_f_dir, path_ref_aln, spanning_band, refine):
    winners = choose_winners_ref(options, path_f_dir, path_ref_aln=path_ref_aln, double_check=True, or_aln=True,
                                 parallel_check=True)
    res = fetch_by_label(labels, winners)
    path_aln = path_f_dir + "/result.aln"
    write_seqs(winners.values(), path_aln)
    write(res, path_f_dir + "/results.unrefined.typ")
    if refine: 
        res = refine_statuses(res)
    write(res, path_f_dir + "/results.tpy")
    return res, path_aln


def assign_site_types(part, gtf, cds):
    if is_acceptor_part(gtf, cds): 
        part["acceptor_sites"] += [int(gtf[3])]
    if is_terminal_part(gtf, cds): 
        part["terminal_sites"] += [int(gtf[4])]
    if is_donor_part(gtf, cds): 
           part["donor_sites"] += [int(gtf[4])]
    return part


def check_for_starts(res, dict_parts, cds, gtf, generegion, minintron, minexon, reject=[]):
    # Get all the available initial parts
    parts = flatten_lol(
        [[copy.deepcopy(p) for p in q.values() if p["initial"] and p["generegion"] == generegion] for q in
         dict_parts.values()])
    labels_k = {}
    for opart in parts:
        part = copy.deepcopy(opart)
        trialer = copy.deepcopy(res[generegion])
        trialparts = [p for p in trialer if int(p["gtfline"][3]) > int(part["gtfline"][4])]
        if trialparts:
            trialparts[0]["left_flat"] = trialparts[0]["lbackupflat"]
            trialparts[0]["options_left"] = trialparts[0]["lbackup"]
            part = initialise_end_options(part)
            part = create_part(part, "fresh", startframe=0, pid="-1")
            gtf = part["gtfline"]
            # Add site type information
            part["initial_sites"] = [int(gtf[3])]
            part = assign_site_types(part, gtf, cds)
            trialparts = [part] + trialparts
            labels_l = careful_part_combinations(trialparts, [], cds, generegion, minintron, minexon, reject=reject)
            for l in labels_l: labels_k[l] = labels_l[l]
        # Also try merging the new part with any old parts that overlap it.
        part = copy.deepcopy(opart)
        trialer = copy.deepcopy(res[generegion])
        # We want to merge the part with the last part that occupies the same frame
        trialparts = [p for p in trialer if int(p["gtfline"][3]) > int(part["gtfline"][4])]
        remaining = [p for p in trialer if not (int(p["gtfline"][3]) > int(part["gtfline"][4]))]
        found = False
        for p in reversed(remaining):
            if gtf_lines_overlap(part["gtfline"], p["gtfline"]) and gtf_compatible_frame(part["gtfline"], p["gtfline"]):
                pp = copy.deepcopy(p)
                pp["status"] = "fresh"
                gtf = part["gtfline"]
                pp["gtfline"][3] = gtf[3]
                pp["options_left"][0] = [int(gtf[3])]
                pp["left_flat"] = [int(gtf[3])]
                pp["initial_sites"] = [int(gtf[3])]
                trialparts = [pp] + trialparts
                found = True
                break
            else:
                trialparts.append(p)
        if found:
            labels_l = careful_part_combinations(trialparts, [], cds, generegion, minintron, minexon, reject=reject)
            for l in labels_l: labels_k[l] = labels_l[l]
    return labels_k


def check_for_stops(dict_cds_bases, dict_parts, labels_k, cds, generegion, minintron, minexon, reject=[]):
    labels_r = {}
    parts = flatten_lol([[copy.deepcopy(p) for p in q.values() if
                          is_terminal_part(p["gtfline"], dict_cds_bases[generegion]) and p["generegion"] == generegion]
                         for q in dict_parts.values()])
    for opart in parts:
        for label in labels_k:
            trialer = copy.deepcopy(labels_k[label])
            # Try removing any end parts and shoehorning in the probe part
            part = copy.deepcopy(opart)
            trialparts = [p for p in trialer["parts"] if int(p["gtfline"][4]) < int(part["gtfline"][3])]
            if trialparts:
                trialparts[-1]["right_flat"] = trialparts[-1]["rbackupflat"]
                trialparts[-1]["options_right"] = trialparts[-1]["rbackup"]
                trialparts[-1]["status"] = "waiting_h"
                part = initialise_end_options(part)
                gtf = part["gtfline"]
                part["options_left"][0] = [int(gtf[3])]
                part["options_right"][get_end_frame(gtf)] = [int(gtf[4])]
                part["left_flat"] = [int(gtf[3])]
                part["right_flat"] = [int(gtf[4])]
                part["status"] = "fresh"
                part["lbackupflat"] = part["left_flat"]
                part["rbackupflat"] = part["right_flat"]
                part["rbackup"] = copy.deepcopy(part["options_right"])
                part["lbackup"] = copy.deepcopy(part["options_left"])
                # Add site type information
                if is_initial_part(gtf, cds): 
                    part["initial_sites"] += [int(gtf[3])]
                if is_acceptor_part(gtf, cds): 
                    part["acceptor_sites"] += [int(gtf[3])]
                part["terminal_sites"] = [int(gtf[4])]
                part["donor_sites"] = []
                part["id"] = max([t["id"] for t in trialparts]) + 1
                trialparts = trialparts + [part]
                labels_l = careful_part_combinations(trialparts, [], cds, generegion, minintron, minexon, reject=reject)
                for l in labels_l: labels_r[l] = labels_l[l]
            # Also try merging the new part with any old parts that overlap it.
            part = copy.deepcopy(opart)
            trialer = copy.deepcopy(labels_k[label])
            # We want to merge the part with the last part that occupies the same frame
            trialparts = [p for p in trialer["parts"] if int(p["gtfline"][4]) < int(part["gtfline"][3])]
            remaining = [p for p in trialer["parts"] if not (int(p["gtfline"][4]) < int(part["gtfline"][3]))]
            found = False
            for p in remaining:
                if gtf_lines_overlap(part["gtfline"], p["gtfline"]) and gtf_compatible_frame(part["gtfline"],
                                                                                             p["gtfline"]):
                    pp = copy.deepcopy(p)
                    pp["status"] = "fresh"
                    gtf = part["gtfline"]
                    pp["gtfline"][4] = gtf[4]
                    pp["options_right"][get_end_frame(gtf)] = [int(gtf[4])]
                    pp["right_flat"] = [int(gtf[4])]
                    pp["terminal_sites"] = [int(part["gtfline"][4])]
                    trialparts.append(pp)
                    found = True
                    break
                else:
                    trialparts.append(p)
            if found:
                labels_l = careful_part_combinations(trialparts, [], cds, generegion, minintron, minexon, reject=reject)
                for l in labels_l: labels_r[l] = labels_l[l]
    return labels_r


def check_for_starts_and_stops(res, dict_cds_bases, dict_parts, generegions, path_dbl_dir, minintron, minexon,
                               d_reject={}):
    # It is possible that the result of this checking process might lack a true initial or terminal exon.
    # We check at this stage whether this has taken place. We then shoehorn in all the options for initial
    # and terminal exons and pick the best ones.
    # In addition to the shoehorning, we accept extensions and contractions of existing parts.
    replace = False
    labels = idict(generegions, {})
    for k in res:
        # Decide whether the gene model needs to be replaced.
        # If the result is blank (i.e. recommend remove gene), then skip this generegion.
        if not res[k]: 
            continue
        # Check if any of the parts are initial
        partssorted = sorted(res[k], key=lambda x: int(x["gtfline"][3]))
        gtf = [p["gtfline"] for p in partssorted]
        cds = dict_cds_bases[k]
        gtfcds = get_cds_live(gtf, cds)
        firstcodon = gtfcds[:3]
        gtfcdsframe = int(gtf[0][7])
        endframe = get_end_frame(gtf[-1])
        # If the first codon is not a start codon, try adding in some starts.
        labels_k = {}
        if not (is_start_codon(firstcodon) and gtfcdsframe == 0):
            replace = True
            labels_k = check_for_starts(res, dict_parts, cds, gtf, k, minintron, minexon, d_reject.get(k, []))
        else:
            aaseq = seq_from_parts_list(res[k], cds)
            labels_k = {
                "stillgood": {"parts": res[k], "aa": aaseq, "seq": make_seq("stillgood", aaseq), "generegion": k}}
        # Check if the last part is terminal
        lastcodon = gtfcds[-3:]
        if not (is_stop_codon(lastcodon) and endframe == 0):
            replace = True
            # Get all the available terminal parts
            labels_r = check_for_stops(dict_cds_bases, dict_parts, labels_k, cds, k, minintron, minexon,
                                       d_reject.get(k, []))
            for l in labels_r: labels_k[l] = labels_r[l]
            if "stillgood" in labels_k: 
                del labels_k["stillgood"]
        labels[k] = labels_k
    # If the gene model does not need to be replaced, simply return it.
    # If it does, then choose between the options for replacement.
    # It is possible that no viable gene model will be found.
    if replace:
        dict_opts = dict((k, [l["seq"] for l in labels[k].values()]) for k in generegions)
        path_chk_dir = make_if_absent(path_dbl_dir + "/check_final")
        winners = choose_winners_seeded(dict_opts, path_chk_dir, seed="auto")
        for generegion in labels:
            if not generegion in winners: 
                continue
            w = labels[generegion][winners[generegion].id]
            for k in w: res[generegion] = w["parts"]
            write_csv([p["gtfline"] for p in res[generegion]], path_chk_dir + "/" + generegion + ".gtf")
        write_seqs(winners.values(), path_chk_dir + "/result.aln")
    return res


def careful_part_combinations(plist, latestparts, cds, generegion, mi, mx, tryempty=False, reject=[], dr=""):
    # Try not adding in the part, if relevant.
    labels_e = get_labels(plist, cds, reject, generegion, mi, mx, True) if tryempty or not latestparts else {}
    # Then try adding in the part.
    labels_l = {}
    if latestparts:
        plist += latestparts
        labels_l = get_labels(plist, cds, reject, generegion, mi, mx, True)
        # Try to foresee problems with non-viable part sets
        if len(plist) > 1:
            last = plist[-1]["gtfline"]
            penlt = plist[-2]["gtfline"]
            # If the last two lines overlap or the last is before the first,
            # try removing each of them. Note we don't need to remove the last
            # one if tryempty is on, that will already have been covered
            # Update: only try removing last part if the part in question is terminal.
            if gtf_lines_overlap(last, penlt) or int(last[4]) < int(penlt[3]):
                if (not tryempty) and plist[-1]["terminal_sites"]:
                    labels_c = get_labels(del_indices(plist, [-1]), cds, reject, generegion, mi, mx, True)
                    for l in labels_c: labels_l[l] = labels_c[l]
                labels_d = get_labels(del_indices(plist, [-2]), cds, reject, generegion, mi, mx, True)
                for l in labels_d: labels_l[l] = labels_d[l]
    # Any time we have a start codon, be sure to squash all the previous stuff
    # and try out the initial part on its own. (The part on its own should indeed
    # suffice - we shouldn't need to retain anything after it.)
    for lp in [a for a in latestparts if a["initial_sites"]]:
        labels_s = get_labels([lp], cds, reject, generegion, mi, mx, True)
        for l in labels_s: labels_l[l] = labels_s[l]
    return merge_dicts([labels_l, labels_e])


def get_labels(plistc, cds, reject, generegion, mi, mx, chst):
    opts, plist2 = get_part_combinations_safe(plistc, cds, reject)
    return get_combinations_fasta(plist2, opts, cds, generegion, mi, mx, check_start=chst)


def get_part_combinations_safe(partslist, cds, reject=[]):
    # The latest "waiting" part and any subsequent ones should be fully laid out.
    plist = copy.deepcopy(partslist)
    if [p["id"] for p in plist if p["status"].startswith("waiting")]:
        max_waiting_id = max([p["id"] for p in plist if p["status"].startswith("waiting")])
        for p in plist:
            if max_waiting_id == p["id"]:
                p["right_flat"] = copy.deepcopy(p["rbackupflat"])
                p["options_right"] = copy.deepcopy(p["rbackup"])
            if max_waiting_id < p["id"]:
                p["status"] = "double_q"
                p["right_flat"] = copy.deepcopy(p["rbackupflat"])
                p["options_right"] = copy.deepcopy(p["rbackup"])
                p["options_left"] = copy.deepcopy(p["lbackup"])
                p["left_flat"] = copy.deepcopy(p["lbackupflat"])
    return get_part_combinations(plist, cds, check_init_term=True, reject=reject), plist


#####################################################
# Execution
#####################################################

def go(path_inf, path_ref, path_results_dir, path_w_dir, minintron, minexon, int_numcores, int_slop_amount):
    """The main body of the program.
    """
    #########################################################
    # Set up basic variables and data containers. Extract
    # information from the genome files.
    #########################################################
    dict_seq_info = read_input_locations(path_inf, path_ref)  # qe
    dict_genome_info = get_genome_info(dict_seq_info)  # qe
    #########################################################
    # Get the base region for each gtf, slopping a certain
    # amount each side. Also get the fasta sequences.
    #########################################################
    sprint("Grabbing cds and aa sequences for inputted transcripts...")
    prepare_sequences(path_w_dir, dict_seq_info, dict_genome_info, int_num_cores, int_slop_amount)  # qe
    #########################################################
    # For each genome, merge together overlapping tight base
    # regions. This will allow us to, bluntly, auto-detect
    # alternative transcripts in a region
    #########################################################
    dict_generegions = prepare_generegions(dict_seq_info, dict_genome_info, path_w_dir, int_num_cores,
                                           int_slop_amount)  # qe
    ################################################################
    # Grab the relative boundaries for the original sequences
    ################################################################
    sprint("Relativising sequences...")
    relativise_sequences(dict_generegions, dict_seq_info)  # qe
    #############################################################
    # For each base, exonerate each protein sequence against it.
    # Then grab the gene parts.
    #############################################################
    go_exonerate(path_w_dir, dict_generegions, dict_seq_info, int_num_cores)  # qe
    #############################################################
    # Align all options and pick the best set
    #############################################################
    winners = prepare_options(dict_generegions)  # qe
    ######################################################
    # Align the original sequences to use as a reference
    ######################################################
    path_or_fa = path_w_dir + "/or.fasta"  # qe
    path_or_aln = path_w_dir + "/or.aln"  # qe
    get_or_aln(dict_generegions, dict_seq_info, path_or_fa, path_or_aln)  # qe
    ####################################################
    # Align the winners and get the coordinates.
    ####################################################
    path_winners_aln = prepare_winner_alignment(path_w_dir, winners, dict_generegions, path_or_fa)  # qe
    dict_part_coords = get_winner_aln_coords(winners, path_winners_aln, dict_generegions,
                                             path_winners_aln + ".flat")  # qe
    ###################################################
    # Get a connectivty graph for these bits.
    ###################################################
    groups, adj = group_parts(dict_part_coords)  # qe
    ###############################################################
    # Fix each pair in turn, as directed by the connectivity graph
    ###############################################################
    parts = dict((k, dict_generegions[k]["winner"]) for k in dict_generegions)  # qe
    res = fix_it(adj, parts, dict_generegions, path_w_dir, minintron, minexon, path_winners_aln)  # qe
    ###############################################################
    # Production code to quickly jump to the doublecheck section.
    ###############################################################
    path_p_dir = make_if_absent(path_w_dir + "/pickle")
    # Pickle.
    pickle.dump(res, open(path_p_dir + "/res.p", "w"))  # qe
    pickle.dump(dict_generegions, open(path_p_dir + "/met.p", "w"))  # qe
    pickle.dump(dict_seq_info, open(path_p_dir + "/seq.p", "w"))  # qe
    # Unpickle
    res = pickle.load(open(path_p_dir + "/res.p", "r"))  # qg
    dict_generegions = pickle.load(open(path_p_dir + "/met.p", "r"))  # qg
    dict_seq_info = pickle.load(open(path_p_dir + "/seq.p", "r"))  # qg
    ##################################################
    # Compare with the original sequences
    ##################################################
    res_o = compare_with_original(res, dict_generegions, dict_seq_info, path_w_dir)  # qg
    ##################################################
    # If the process yields nothing, retain the
    # original gene model.
    ##################################################
    res_o = revive_casualties(res_o, dict_generegions, dict_seq_info, path_w_dir)
    res_o = check_aq_improvement(res_o, dict_generegions, dict_seq_info, path_w_dir,
                                 alignment_score(read_seqs(path_or_aln)))
    ##################################################
    # Output the gtfs, fastas, etc. etc. etc.
    ##################################################
    dict_generegions = write_output(dict_generegions, dict_seq_info, res_o, path_results_dir)  # qg
    get_stats(dict_generegions, dict_seq_info, res, path_results_dir, path_w_dir)


if __name__ == '__main__':
    """OMGene - mutual improvement of gene models through gene orthology
    """
    # Check the shell...
    check_shell()

    # Read in the command line arguments
    parser = argparse.Argument_parser(description="Run OMGene")
    parser.add_argument("-i", "--info_files", metavar="infiles", dest="IN")
    parser.add_argument("-r", "--reference_files", metavar="reffiles", dest="RF", default="")
    parser.add_argument("-o", "--out_dir", metavar="outdir", dest="OD")
    parser.add_argument("-g", "--gap_penalty", metavar="outdir", dest="GP")
    parser.add_argument("--minintron", metavar="minintron", dest="MNI", default=4)
    parser.add_argument("--minexon", metavar="minexon", dest="MNE", default=4)
    parser.add_argument("-c", "--cores", metavar="numcores", dest="CO", default=1)
    parser.add_argument("-b", "--buffer", metavar="bufferamount", dest="SL", default=600)
    parser.add_argument("--safealign", action="store_true", dest="SA", default=False)
    parser.add_argument("--donors", "-do", metavar="donors", dest="DO", default="gt,gc")
    parser.add_argument("--acceptors", "-ac", metavar="acceptors", dest="AC", default="ag")
    parser.add_argument("--startcodon", "-sc", metavar="startcodons", dest="SC", default="atg")
    args = parser.parse_args()
    # CUPCAKES: dependencies: pyfaidx, exonerate, blastp
    # Check existence and non-confliction of arguments
    path_inf = args.IN
    path_ref = args.RF
    path_out_dir = args.OD
    if not path_inf: 
        sys.exit("Please specify an input file using the -i option")
    minexon = int(args.MNE)
    minintron = int(args.MNI)
    int_num_cores = int(args.CO)
    int_slop_amount = int(args.SL)
    static_safe_align = args.SA

    static_do, static_ac, static_sc = process_codon_options(args.DO, args.AC, args.SC)

    if args.GP:
        static_blosdict = blosum_dict(-1 * int(args.GP), 0)
        static_blosmat = blosum_matrix(static_blosdict)
        static_blospos = blosum_positions()

    path_results_dir, path_w_dir = prepare_output_folder(path_out_dir)
    go(path_inf, path_ref, path_results_dir, path_w_dir, minintron, minexon, int_num_cores, int_slop_amount)
