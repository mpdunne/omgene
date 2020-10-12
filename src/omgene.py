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

import multiprocessing
import datetime
import pickle
import string
import pyfaidx
import argparse
import random

from utils_todo.check_shell import *
from utils.misc import *
from utils.sequence_io import *

from utils_todo.part_groups import *
from utils.parts import *
from utils_todo.alignment import *
from utils_todo.exonerate import *
from utils_todo.choice import *
from utils_todo.fixing import *
from utils_todo.filtering import *

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
    return {**labels_l, **labels_e}


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


def process_codon_option(codon_option, expected_length):
    """
    Process a proposed codon option or set of options (separated by commas).
    This function checks that the proposed codon is valid and checks that it is the
    correc length for what is expected. Returns none if co is none.

    :param codon_option: (str) Codon option, in string format, or comma-separated list of options
    :param expected_length): (int) The expected length of the option, e.g. 3 for a start codon.
    """
    if codon_option is None:
        return None

    codon_options = codon_option.split(",")
    options = []

    for a in codon_options:
        a = re.sub(r" ", r"", a)

        if not len(a) == expected_length or any(not s.lower() in "acgtu" for s in a):
            sys.exit("Invalid start codon/splice site choice: " + str(a))

        options.append(re.sub(r"u", r"t", a.lower()))

    return options


def process_codon_options(do, ac, sc):
    s_do = process_codon_option(do, 2)
    s_ac = process_codon_option(ac, 2)
    s_sc = process_codon_option(sc, 3)
    return s_do, s_ac, s_sc


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
    parser.add_argument("-g", "--gap_penalty", metavar="outdir", dest="GP", default=1)
    parser.add_argument("--minintron", metavar="minintron", dest="MNI", default=4)
    parser.add_argument("--minexon", metavar="minexon", dest="MNE", default=4)
    parser.add_argument("-c", "--cores", metavar="numcores", dest="CO", default=1)
    parser.add_argument("-b", "--buffer", metavar="bufferamount", dest="SL", default=600)
    parser.add_argument("--safealign", action="store_true", dest="SA", default=False)
    parser.add_argument("--donors", "-do", metavar="donors", dest="DO")
    parser.add_argument("--acceptors", "-ac", metavar="acceptors", dest="AC")
    parser.add_argument("--startcodon", "-sc", metavar="startcodons", dest="SC")
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

    codon_options = process_codon_options(args.DO, args.AC, args.SC)
    initialise_features(*codon_options)
    initialise_blosum_matrix(-1 * int(args.GP), 0)

    path_results_dir, path_w_dir = prepare_output_folder(path_out_dir)
    go(path_inf, path_ref, path_results_dir, path_w_dir, minintron, minexon, int_num_cores, int_slop_amount)
