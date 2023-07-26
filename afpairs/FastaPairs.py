import argparse
import os

from Bio import SeqIO, SeqRecord

from afpairs import FastaId


def dir_path(string):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(
        description="Create FASTA files, each one containing a protein pair, "
                    "one protein from baits and one protein from targets file.")
    parser.add_argument('-b', '--baits', type=argparse.FileType('r'), required=True,
                        help="FASTA file containing baits")
    parser.add_argument('-t', '--targets', type=argparse.FileType('r'), required=True,
                        help="FASTA file containing targets")
    parser.add_argument('-u', '--unique', action='store_true',
                        help="Save only one FASTA file per pair "
                             "- do not save POLR2B-POLR2A pair if POLR2A-POLR2B is also present.")
    parser.add_argument('-i', '--identity', action='store_true',
                        help="Don't save FASTA file of a protein with itself "
                             "- same protein is present in both baits and targets.")
    parser.add_argument('-o', '--output', type=dir_path, default="",
                        help="Directory where to write FASTA files.  (default: current directory)")

    args = parser.parse_args(argv)

    fasta_pairs(baits=args.baits, targets=args.targets, unique=args.unique, skip_identity=args.identity,
                output=args.output)


def fasta_pairs(baits: str, targets: str, unique: bool = False, skip_identity: bool = False, output: str = ""):
    """
    Create FASTA files, each one containing a protein pair, one protein from baits and one protein from targets file.

    :param baits: baits
    :param targets: targets
    :param unique: save only one FASTA file per unique pair -
                   do not save POLR2B-POLR2A pair if POLR2A-POLR2B is also present
    :param skip_identity: don't save FASTA file of a protein with itself -
                          if the same protein is present in both baits and targets
    :param output: where to write FASTA files
    """
    baits = parse_fasta(baits)
    targets = parse_fasta(targets)

    processed_ids = set()
    for bait in baits:
        for target in targets:
            if skip_identity and bait == target:
                continue
            reserve_id = f"{target}__{bait}"
            if unique and reserve_id in processed_ids:
                continue
            merge_id = f"{bait}__{target}"
            SeqIO.write([baits[bait], targets[target]], os.path.join(output, f"{merge_id}.fasta"), "fasta")
            processed_ids.add(merge_id)


def parse_fasta(fasta: str) -> dict[str, SeqRecord]:
    """
    Parses FASTA and returns all sequences found in file mapped by ID.

    The ID of each sequence is found using :func:`FastaId.fasta_id`

    :param fasta: FASTA file
    :return: all sequences found in file mapped by ID
    """
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq_id = FastaId.fasta_id(record.description)
        sequences[seq_id] = record
    return sequences


if __name__ == '__main__':
    main()
