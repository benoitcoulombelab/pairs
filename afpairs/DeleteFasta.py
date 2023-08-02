import argparse
import os
import re
import shutil

from Bio import SeqIO

VALID_AMINO_ACID_REGEX = r"^[ARNDCQEGHILKMFPSTWYV]+$"


def dir_path(string):
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
    parser.add_argument('-b', '--backup', type=dir_path, default=None,
                        help="Copy FASTA files that are to be deleted to this folder for backup")

    args = parser.parse_args(argv)

    delete_fasta(inputs=args.inputs, invalid_sequence=args.sequence, length=args.length, backup=args.backup)


def delete_fasta(inputs: list[str] = [], invalid_sequence: bool = False, length: int = None, backup: str = None):
    """
    Deletes FASTA files if they fail any condition.

    :param inputs: FASTA files
    :param invalid_sequence: delete if sequence is invalid
    :param length: maximum length of the sum of all sequence lengths
    :param backup: FASTA files that are to be deleted will be moved to this folder instead
    """
    for fasta in inputs:
        delete = False
        sequences = parse_fasta(fasta)
        if invalid_sequence:
            delete = delete or [sequence for sequence in sequences if
                                re.match(VALID_AMINO_ACID_REGEX, str(sequence.seq)) is None]
        if length:
            fasta_length = sum([len(sequence) for sequence in sequences])
            delete = delete or fasta_length > length
        if delete:
            if backup:
                shutil.copyfile(fasta, os.path.join(backup, os.path.basename(fasta)))
            os.remove(fasta)


def parse_fasta(fasta: str) -> list:
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