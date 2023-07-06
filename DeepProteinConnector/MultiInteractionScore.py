import argparse
import os
import re
import sys

from DeepProteinConnector import InteractionScore


def file_path(string):
    if not string or os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


def main():
    parser = argparse.ArgumentParser(description="Count interactions for multiple PDB files.")
    parser.add_argument('infile', nargs='+', type=file_path,
                        help="PDB files for which to count interactions")
    parser.add_argument('-n', '--name', default=r"(\w+)__(\w+)",
                        help="Regular expression to obtain protein/gene names based on PDB filename "
                             " (default: %(default)s)")
    parser.add_argument('-r', '--radius', type=float, default=6.0,
                        help="Distance between atoms from different residues to assume interaction "
                             "(default: %(default)s)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing counts")

    args = parser.parse_args()

    args.output.write("Bait\tTarget\tCount\n")
    for infile in args.infile:
        re_match = re.search(args.name, infile)
        if not re_match:
            raise AssertionError(f"Expression {args.name} cannot be found in filename {infile}")
        bait = re_match.group(1)
        target = re_match.group(2)
        count = interaction_score(infile, args.radius)
        args.output.write(f"{bait}\t{target}\t{count}\n")


def interaction_score(pdb: "PDB file", radius: float = 6) -> int:
    """Count number of interactions in PDB file"""
    return InteractionScore.interaction_score(pdb, radius=radius)


if __name__ == '__main__':
    main()
