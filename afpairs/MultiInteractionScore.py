import argparse
import os
import re
import sys

from afpairs import InteractionScore


def file_path(string):
    if not string or os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


def main():
    parser = argparse.ArgumentParser(description="Compute protein-protein interaction score from multiple PDB files.")
    parser.add_argument('infile', nargs='+', type=file_path,
                        help="PDB files for which to compute PPI score")
    parser.add_argument('-n', '--name', default=r"(\w+)__(\w+)",
                        help="Regular expression to obtain protein/gene names based on PDB filename "
                             " (default: %(default)s)")
    parser.add_argument('-a', '--first', default="A",
                        help="Chains of first protein separated by ','  (default: %(default)s)")
    parser.add_argument('-b', '--second', default="B",
                        help="Chains of second protein separated by ','  (default: %(default)s)")
    parser.add_argument('-r', '--radius', type=float, default=6.0,
                        help="Distance between atoms from different residues to assume interaction "
                             "(default: %(default)s)")
    parser.add_argument('-w', '--weight', action="store_true", default=False,
                        help="Normalize count by protein pair weight - "
                             "'count / log2(sum weight of both proteins)'")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing counts")

    args = parser.parse_args()
    args.first = args.first.split(',')
    args.second = args.second.split(',')

    args.output.write("Bait\tTarget\tScore\n")
    for infile in args.infile:
        re_match = re.search(args.name, infile)
        if not re_match:
            raise AssertionError(f"Expression {args.name} cannot be found in filename {infile}")
        bait = re_match.group(1)
        target = re_match.group(2)
        count = InteractionScore.interaction_score(pdb=infile, first_chains=args.first, second_chains=args.second,
                                                   radius=args.radius, weight=args.weight)
        args.output.write(f"{bait}\t{target}\t{count}\n")


if __name__ == '__main__':
    main()
