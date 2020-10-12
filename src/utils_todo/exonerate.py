import os
import subprocess
from utils_todo.gtf_tools import *


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

