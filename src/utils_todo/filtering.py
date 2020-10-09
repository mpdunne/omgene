
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
