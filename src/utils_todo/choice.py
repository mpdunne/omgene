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

