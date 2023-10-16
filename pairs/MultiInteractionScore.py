import argparse
import os
import re
import sys
import glob
import statistics
import json
from typing import TextIO

import smokesignal
import tqdm

from pairs import InteractionScore


class AlphafoldStatistics:
    def __init__(self, directory: str, ranked0_score: float, ranked0_confidence: float, average_confidence: float,
                 unrelaxed_average_score: float, unrelaxed_score_standard_deviation: float):
        self.directory = directory
        self.ranked0_score = ranked0_score
        self.ranked0_confidence = ranked0_confidence
        self.average_confidence = average_confidence
        self.unrelaxed_average_score = unrelaxed_average_score
        self.unrelaxed_score_standard_deviation = unrelaxed_score_standard_deviation


def file_path(string: str):
    if not string or os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Compute protein-protein interaction score from multiple PDB files.")
    parser.add_argument('inputs', nargs='+', type=file_path,
                        help="PDB files for which to compute PPI score")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing scores")
    parser.add_argument('-n', '--name', default=r"([\w-]+)__([\w-]+)",
                        help="Regular expression to obtain protein/gene names based on PDB filename "
                             " (default: %(default)s)")
    parser.add_argument('-s', '--stats', action="store_true", default=False,
                        help="Compute confidence and unrelaxed score statistics")
    parser.add_argument('-a', '--first', default="A",
                        help="Chains of first protein separated by ','  (default: %(default)s)")
    parser.add_argument('-b', '--second', default="B",
                        help="Chains of second protein separated by ','  (default: %(default)s)")
    parser.add_argument('-r', '--radius', type=float, default=6.0,
                        help="Distance between atoms from different residues to assume interaction "
                             "(default: %(default)s)")
    parser.add_argument('-w', '--weight', action="store_true", default=False,
                        help="Normalize count by protein pair weight - "
                             "'score / log2(sum weight of both proteins)'")
    parser.add_argument('-c', '--count', action="store_true", default=False,
                        help="Score is the number of residues pairs at a distance less than radius parameter")
    parser.add_argument('-p', '--progress', action="store_true", default=False,
                        help="Show progress bar")
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

    multi_interaction_score(input_files=args.inputs, output_file=args.output, name=args.name,
                            stats=args.stats,
                            radius=args.radius, weight=args.weight, count=args.count, progress=args.progress,
                            first_chains=args.first, second_chains=args.second,
                            partial=args.partial,
                            mapping_file=args.mapping, source_column=args.source_column - 1,
                            converted_column=args.converted_column - 1)


def multi_interaction_score(input_files: list[str], output_file: TextIO = sys.stdout, name: str = r"([\w-]+)__([\w-]+)",
                            stats: bool = False,
                            radius: float = 6, weight: bool = False, count: bool = False, progress: bool = False,
                            first_chains: list[str] = ["A"], second_chains: list[str] = ["B"],
                            partial: bool = False,
                            mapping_file: TextIO = None, source_column: int = 0, converted_column: int = 1):
    """
    Compute protein-protein interaction score from multiple PDB files.

    :param input_files: PDB files
    :param output_file: output file
    :param stats: compute confidence and unrelaxed score statistics
    :param name: regular expression to obtain protein/gene names based on PDB filename
    :param radius: maximal distance between two residues' atoms to consider that the two residues interact
    :param weight: if True, normalize score by proteins' weight
    :param count: if True, score is the number of residue pairs below radius
    :param progress: if True, show progress bar
    :param first_chains: chains of the first protein
    :param second_chains: chains of the second protein
    :param partial: do not warn if all chains in PDB are not used for computing score
    :param mapping_file: tab delimited text file used to convert names
    :param source_column: column index of source names in mapping file
    :param converted_column: column index of converted names in mapping file
    """
    mappings = {}
    if mapping_file:
        mappings = parse_mapping(mapping_file, source_column, converted_column)
    output_file.write("Bait\tTarget\t")
    if stats:
        output_file.write("Ranked_0 score\tRanked_0 confidence\tAverage confidence\t"
                          "Unrelaxed average score\t Unrelaxed score standard deviation\n")
    else:
        output_file.write("Score\n")
    for input_file in (tqdm.tqdm(input_files) if progress else input_files):
        re_match = re.search(name, input_file)
        if not re_match:
            raise AssertionError(f"Expression {name} cannot be found in filename {input_file}")
        bait, target = re_match.group(1, 2)
        bait = mappings[bait] if bait in mappings else bait
        target = mappings[target] if target in mappings else target
        output_file.write(f"{bait}\t{target}\t")
        if stats:
            a_s = alphafold_statistics(directory=os.path.dirname(input_file), radius=radius, weight=weight, count=count,
                                       first_chains=first_chains, second_chains=second_chains, partial=partial)
            output_file.write(f"{a_s.ranked0_score}\t{a_s.ranked0_confidence}\t{a_s.average_confidence}\t"
                              f"{a_s.unrelaxed_average_score}\t{a_s.unrelaxed_score_standard_deviation}\n")
        else:
            with open(input_file, 'r') as input_in:
                score = InteractionScore.interaction_score(
                    pdb=input_in, radius=radius, weight=weight, count=count,
                    first_chains=first_chains, second_chains=second_chains, partial=partial)
            output_file.write(f"{score}\n")


def alphafold_statistics(directory: str, radius: float = 6, weight: bool = False, count: bool = False,
                         first_chains: list[str] = ["A"], second_chains: list[str] = ["B"],
                         partial: bool = False) -> AlphafoldStatistics:
    """
    Computes statistics on AlphaFold's output.

    :param directory: AlphaFold's output
    :param radius: maximal distance between two residues' atoms to consider that the two residues interact
    :param weight: if True, normalize score by proteins' weight
    :param count: if True, score is the number of residue pairs below radius
    :param first_chains: chains of the first protein
    :param second_chains: chains of the second protein
    :param partial: do not warn if all chains in PDB are not used for computing score
    :return: statistics on AlphaFold's output
    """
    ranked0_file = os.path.join(directory, "ranked_0.pdb")
    if not os.path.isfile(ranked0_file):
        raise AssertionError(f"ranked_0.pdb does not exists in directory {directory}")
    with open(ranked0_file, 'r') as input_in:
        ranked0_score = InteractionScore.interaction_score(
            pdb=input_in, radius=radius, weight=weight, count=count,
            first_chains=first_chains, second_chains=second_chains, partial=partial)
    ranking_file = os.path.join(directory, "ranking_debug.json")
    ranked0_confidence = None
    average_confidence = None
    if os.path.isfile(ranking_file):
        with open(ranking_file, 'r') as input_in:
            rankings = json.load(input_in)
        ranked0_model = rankings["order"][0]
        ranked0_confidence = rankings["iptm+ptm"][ranked0_model]
        average_confidence = statistics.mean(rankings["iptm+ptm"].values())
    unrelaxed_files = glob.glob(os.path.join(directory, "unrelaxed_*.pdb"))
    unrelaxed_average_score = None
    unrelaxed_score_standard_deviation = None
    if unrelaxed_files:
        unrelaxed_scores = []
        for unrelaxed_file in unrelaxed_files:
            with open(unrelaxed_file, 'r') as input_in:
                unrelaxed_scores.append(InteractionScore.interaction_score(
                    pdb=input_in, radius=radius, weight=weight, count=count,
                    first_chains=first_chains, second_chains=second_chains, partial=partial))
        unrelaxed_average_score = statistics.mean(unrelaxed_scores)
        unrelaxed_score_standard_deviation = statistics.pstdev(unrelaxed_scores)
    return AlphafoldStatistics(
        directory,
        ranked0_score=ranked0_score, ranked0_confidence=ranked0_confidence,
        average_confidence=average_confidence,
        unrelaxed_average_score=unrelaxed_average_score,
        unrelaxed_score_standard_deviation=unrelaxed_score_standard_deviation)


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
