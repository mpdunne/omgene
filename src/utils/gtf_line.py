class GtfLine:
    """
    A class to contain a single line of a GTF.
    """
    def __init__(self, gtf_line, updates=None):
        if updates is None:
            updates = {}

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

    def list(self):
        return [
            self.chr,
            self.src,
            self.feature,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attribute
        ]

    def length(self):
        return self.end - self.start + 1

    def abs_frame(self):
        return (self.start + self.frame) % 3

    def proto_exon(self):
        return [self.start, self.end, self.frame]


def gtf_compatible_frame(gtfline1, gtfline2):
    """
    Returns whether two gtf lines are in the same absolute frame.
    """
    return gtfline1.abs_frame == gtfline2.abs_frame


def overlap_in_frame(first, second):
    """
    Returns whether the two gtf lines overlap and are in the same frame.
    """
    return gtf_lines_overlap(first, second) and first.abs_frame == second.abs_frame


def gtf_lines_overlap(first, second):
    """
    Returns whether the two gtf lines overlap.
    """
    if first.start <= second.start <= first.end:
        return True

    if second.start <= first.start <= second.end:
        return True

    return False


def gtf_line_equals(gtf1, gtf2):
    """
    Returns whether the two gtf lines are equal in terms of chromosome and coordinates.
    """
    return (gtf1.chr, gtf1.start, gtf1.end) == (gtf2.chr, gtf2.start, gtf2.end)
