from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_seqs(path_fasta_in, tolist=True):
    """
    Read in sequences from a FASTA file.

    :param path_fasta_in: The path to the FASTA file we want to read.
    :param tolist: Whether or not to convert the results into list format.
    :return: A collection of SeqRecords, in list format if requested.
    """
    seqs = SeqIO.parse(path_fasta_in, "fasta")
    return list(seqs) if tolist else seqs


def write_seqs(seqs, path_out):
    """
    Write the provided list of sequences to a FASTA file.

    :param seqs: A collection of SeqRecord objects
    :param path_out: The path of the FASTA file to write out.
    """
    SeqIO.write(seqs, path_out, "fasta")


def make_seq(label, sequence):
    """
    Make a SeqRecord object from a label and a sequence.

    :param label: (str) The label to give the sequence.
    :param sequence: (str) The sequence in string format.
    """
    return SeqRecord(id=label, name="", description="", seq=Seq(sequence))
