import argparse
import os
import shutil

from Bio import SeqIO
from afpairs import FastaId


def dir_path(string):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():
    parser = argparse.ArgumentParser(description="Deletes FASTA files that fails any specified limits.")
    parser.add_argument('inputs', nargs='*',
                        help="FASTA files")
    parser.add_argument('-l', '--length', type=int, default=None,
                        help="Maximum length")
    parser.add_argument('-b', '--backup', type=dir_path, default=None,
                        help="Copy FASTA files that are to be deleted to this folder for backup")

    args = parser.parse_args()

    for fasta in args.inputs:
        delete = False
        sequences = parse_fasta(fasta)
        if args.length:
            length = sum([len(sequences[seq_id]) for seq_id in sequences])
            if length > args.length:
                delete = True
        if delete:
            if args.backup:
                shutil.copyfile(fasta, os.path.join(args.backup, os.path.basename(fasta)))
            os.remove(fasta)


def parse_fasta(fasta):
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq_id = FastaId.fasta_id(record.description)
        sequences[seq_id] = record
    return sequences


if __name__ == '__main__':
    main()
