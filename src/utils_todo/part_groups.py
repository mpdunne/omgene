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

