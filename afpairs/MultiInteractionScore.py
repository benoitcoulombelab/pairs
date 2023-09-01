import argparse
import os
import re
import sys
from typing import TextIO

import smokesignal

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
    parser.add_argument('-P', '--partial', action="store_true", default=False,
                        help="Do not warn if all chains in PDB are not used for computing score")
    parser.add_argument('-M', '--mapping', type=argparse.FileType('r'),
                        help="Tab delimited text file used to convert names  (default: %(default)s)")
    parser.add_argument('-S', '--source_column', type=int, default='1',
                        help="Column index of source names in mapping file - 1 means first column of file" +
                             "   (default: %(default)s)")
    parser.add_argument('-C', '--converted_column', type=int, default='2',
                        help="Column index of converted names in mapping file - 1 means first column of file" +
                             "   (default: %(default)s)")

    args = parser.parse_args(argv)
    args.first = args.first.split(',')
    args.second = args.second.split(',')
    smokesignal.on(InteractionScore.MISSING_CHAIN_EVENT, InteractionScore.warn_missing_chain, max_calls=1)

    multi_interaction_score(input_files=args.inputs, name=args.name, radius=args.radius, weight=args.weight,
                            first_chains=args.first, second_chains=args.second, output_file=args.output,
                            partial=args.partial,
                            mapping_file=args.mapping, source_column=args.source_column - 1,
                            converted_column=args.converted_column - 1)


def multi_interaction_score(input_files: list[str], name: str = r"(\w+)__(\w+)",
                            radius: float = 6, weight: bool = False,
                            first_chains: list[str] = ["A"], second_chains: list[str] = ["B"],
                            output_file: TextIO = sys.stdout, partial: bool = False,
                            mapping_file: TextIO = None, source_column: int = 0, converted_column: int = 1):
    """
    Compute protein-protein interaction score from multiple PDB files.

    :param input_files: PDB files
    :param name: regular expression to obtain protein/gene names based on PDB filename
    :param radius: maximal distance between two residues' atoms to consider that the two residues interact
    :param weight: if True, normalize score by proteins' weight
    :param first_chains: chains of the first protein
    :param second_chains: chains of the second protein
    :param partial: do not warn if all chains in PDB are not used for computing score
    :param mapping_file: tab delimited text file used to convert names
    :param source_column: column index of source names in mapping file
    :param converted_column: column index of converted names in mapping file
    :param output_file: output file
    """
    mappings = {}
    if mapping_file:
        mappings = parse_mapping(mapping_file, source_column, converted_column)
    output_file.write("Bait\tTarget\tScore\n")
    for input_file in input_files:
        re_match = re.search(name, input_file)
        if not re_match:
            raise AssertionError(f"Expression {name} cannot be found in filename {input_file}")
        bait, target = re_match.group(1, 2)
        bait = mappings[bait] if bait in mappings else bait
        target = mappings[target] if target in mappings else target
        with open(input_file, 'r') as input_in:
            score = InteractionScore.interaction_score(
                pdb=input_in, radius=radius, weight=weight,
                first_chains=first_chains, second_chains=second_chains, partial=partial)
        output_file.write(f"{bait}\t{target}\t{score}\n")


def parse_mapping(mapping_file: TextIO, source_column: int = 0, converted_column: int = 1) \
        -> dict[str, str]:
    """
    Parse mapping file.

    :param mapping_file: text delimited filename
    :param source_column: index of source id columns
    :param converted_column: index of converted id columns
    :return: dictionary of source id to converted id
    """
    mappings = {}
    for line in mapping_file:
        if line.startswith('#'):
            continue
        columns = line.rstrip('\r\n').split('\t')
        source = columns[source_column]
        converted = columns[converted_column]
        if converted:
            mappings[source] = converted
    return mappings


if __name__ == '__main__':
    main()
