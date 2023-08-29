import argparse
import os
import re
import sys
from typing import TextIO

from afpairs import InteractionScore


def file_path(string):
    if not string or os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Compute protein-protein interaction score from multiple PDB files.")
    parser.add_argument('inputs', nargs='+', type=file_path,
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

    args = parser.parse_args(argv)
    args.first = args.first.split(',')
    args.second = args.second.split(',')
    multi_interaction_score(input_files=args.inputs, name=args.name, radius=args.radius, weight=args.weight,
                            first_chains=args.first, second_chains=args.second, output_file=args.output)


def multi_interaction_score(input_files: list[str], name: str = r"(\w+)__(\w+)",
                            radius: float = 6, weight: bool = False,
                            first_chains: list[str] = ["A"],
                            second_chains: list[str] = ["B"],
                            output_file: TextIO = sys.stdout):
    """
    Compute protein-protein interaction score from multiple PDB files.

    :param input_files: PDB files
    :param name: regular expression to obtain protein/gene names based on PDB filename
    :param radius: maximal distance between two residues' atoms to consider that the two residues interact
    :param weight: if True, normalize score by proteins' weight
    :param first_chains: chains of the first protein
    :param second_chains: chains of the second protein
    :param output_file: output file
    """
    output_file.write("Bait\tTarget\tScore\n")
    for input_file in input_files:
        re_match = re.search(name, input_file)
        if not re_match:
            raise AssertionError(f"Expression {name} cannot be found in filename {input_file}")
        bait, target = re_match.group(1, 2)
        with open(input_file, 'r') as input_in:
            score = InteractionScore.interaction_score(
                pdb=input_in, radius=radius, weight=weight,
                first_chains=first_chains, second_chains=second_chains)
        output_file.write(f"{bait}\t{target}\t{score}\n")


if __name__ == '__main__':
    main()
