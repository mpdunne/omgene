import itertools
import networkx as nx
import math
from Bio.Seq import Seq


class GtfLine:

    def __init__(self, gtf_line, updates={}):
        if type(gtf_line) is str:
            gtf_line = gtf_line.split("\t")

        if len(gtf_line) != 9:
            raise ValueError("GTF line length must be 9.")

        self.chr = updates.get('chr', gtf_line[0])
        self.src = updates.get('src', gtf_line[1])
        self.feature = updates.get('feature', gtf_line[2])
        self.start = int(updates.get('start', gtf_line[3]))
        self.end = int(updates.get('end', gtf_line[4]))
        self.score = updates.get('score', gtf_line[5])
        self.strand = updates.get('strand', gtf_line[6])
        self.frame = int(updates.get('frame', gtf_line[7]))
        self.attribute = updates.get('attribute', gtf_line[8])


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
