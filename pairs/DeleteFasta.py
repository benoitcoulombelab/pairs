import argparse
import os
import re
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

VALID_AMINO_ACID_REGEX = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def dir_path(string: str):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Deletes FASTA files that fails any specified limits.")
    parser.add_argument('inputs', nargs='*',
                        help="FASTA files")
    parser.add_argument('-s', '--sequence', action='store_true',
                        help="Sequence must be valid")
    parser.add_argument('-l', '--length', type=int, default=None,
                        help="Maximum length")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Verbose - show name of delete FASTA files")
    parser.add_argument('-b', '--backup', type=dir_path, default=None,
                        help="Copy FASTA files that are to be deleted to this folder for backup")

    args = parser.parse_args(argv)

    delete_fasta(inputs=args.inputs, invalid_sequence=args.sequence, length=args.length, backup=args.backup,
                 verbose=args.verbose)


def delete_fasta(inputs: list[str], invalid_sequence: bool = False, length: int = None, backup: str = None,
                 verbose: bool = False):
    """
    Deletes FASTA files if they fail any condition.

    :param inputs: FASTA files
    :param invalid_sequence: delete if sequence is invalid
    :param length: maximum length of the sum of all sequence lengths
    :param backup: FASTA files that are to be deleted will be moved to this folder instead
    :param verbose: print why files are deleted
    """
    for fasta in inputs:
        delete = False
        sequences = parse_fasta(fasta)
        if invalid_sequence:
            any_invalid = [sequence for sequence in sequences if
                           re.match(VALID_AMINO_ACID_REGEX, str(sequence.seq)) is None]
            if any_invalid:
                if verbose:
                    print(f"Invalid sequence {any_invalid[0].name} in FASTA file {fasta}")
                delete = True
        if length is not None:
            fasta_length = sum([len(sequence) for sequence in sequences])
            if fasta_length > length:
                if verbose:
                    print(f"Sequence length {fasta_length} above maximum {length} in FASTA file {fasta}")
                delete = True
        if delete:
            if backup:
                shutil.copyfile(fasta, os.path.join(backup, os.path.basename(fasta)))
            os.remove(fasta)


def parse_fasta(fasta: str) -> list[SeqRecord]:
    """
    Returns sequences read from FASTA file.

    :param fasta: FASTA file
    :return: sequences read from FASTA file
    """
    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        sequences.append(record)
    return sequences


if __name__ == '__main__':
    main()
