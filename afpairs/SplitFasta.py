import argparse
import os
import sys
from typing import TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from afpairs import FastaId


def dir_path(string):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Create FASTA files containing one sequence per file.")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file containing baits")
    parser.add_argument('-o', '--output', type=dir_path, default="",
                        help="Directory where to write FASTA files.  (default: current directory)")

    args = parser.parse_args(argv)
    split_fasta(fasta=args.input, output_dir=args.output)


def split_fasta(fasta: TextIO, output_dir: str):
    """
    Splits FASTA file into one file per sequence.

    The name of the output file will be the "{id}.fasta where" "{id}" is the id returned by
    :func:`~afpairs.FastaId.fasta_id`

    :param fasta: FASTA file containing multiple sequences
    :param output_dir: directory where to write files
    """
    sequences = parse_fasta(fasta)

    for sequence_id in sequences:
        SeqIO.write(sequences[sequence_id], os.path.join(output_dir, f"{sequence_id}.fasta"), "fasta")


def parse_fasta(fasta: TextIO) -> list[SeqRecord]:
    """
    Parses FASTA file.

    :param fasta: FASTA file
    :return: list of sequences from FASTA file
    """
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq_id = FastaId.fasta_id(record.description)
        sequences[seq_id] = record
    return sequences


if __name__ == '__main__':
    main()
