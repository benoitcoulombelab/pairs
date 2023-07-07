import argparse
import os

from Bio import SeqIO
from afpairs import FastaId


def dir_path(string):
    if not string or os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main():
    parser = argparse.ArgumentParser(description="Create FASTA files containing protein pairs from a "
                                                 "baits and a targets file.")
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

    args = parser.parse_args()

    baits = parse_fasta(args.baits)
    targets = parse_fasta(args.targets)

    processed_ids = set()
    for bait in baits:
        for target in targets:
            if args.identity and bait == target:
                continue
            reserve_id = f"{target}__{bait}"
            if args.unique and reserve_id in processed_ids:
                continue
            merge_id = f"{bait}__{target}"
            SeqIO.write([baits[bait], targets[target]], os.path.join(args.output, f"{merge_id}.fasta"), "fasta")
            processed_ids.add(merge_id)


def parse_fasta(fasta):
    sequences = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seq_id = FastaId.fasta_id(record.description)
        sequences[seq_id] = record
    return sequences


if __name__ == '__main__':
    main()
