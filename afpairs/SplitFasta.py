import argparse
import os
import sys

from Bio import SeqIO
from afpairs import FastaId


def dir_path(string):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():
    parser = argparse.ArgumentParser(description="Create FASTA files containing one sequence per file.")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file containing baits")
    parser.add_argument('-o', '--output', type=dir_path, default="",
                        help="Directory where to write FASTA files.  (default: current directory)")

    args = parser.parse_args()

    sequences = parse_fasta(args.input)

    for sequence_id in sequences:
        SeqIO.write(sequences[sequence_id], os.path.join(args.output, f"{sequence_id}.fasta"), "fasta")


def parse_fasta(fasta):
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq_id = FastaId.fasta_id(record.description)
        sequences[seq_id] = record
    return sequences


if __name__ == '__main__':
    main()
