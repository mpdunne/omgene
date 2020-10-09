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
