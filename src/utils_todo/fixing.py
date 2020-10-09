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
