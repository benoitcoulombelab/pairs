import argparse
import sys
from typing import TextIO

import pandas


class Interaction:
    def __init__(self, bait: str, target: str, score: float):
        self.bait = bait
        self.target = target
        self.score = score


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Generate score matrix from MultiInteractionScore's output.")
    parser.add_argument('scores', nargs='*', type=argparse.FileType('r'), default=[sys.stdin],
                        help="Tab delimited files containing protein id/accession and scores")
    parser.add_argument('-u', '--unique', action='store_true',
                        help="Protein pairs are unique, so mirror the scores in matrix")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing score matrix")

    args = parser.parse_args(argv)
    score_matrix(score_files=args.scores, output_file=args.output, unique=args.unique)


def score_matrix(score_files: list[TextIO], output_file: TextIO, unique: bool = False):
    """
    Generates score matrix.
    Matrix rows contain one bait and columns contain one target.

    :param score_files: tab delimited files containing protein id/accession and scores
    :param output_file: tab delimited output file containing score matrix
    :param unique: protein pairs are unique, so mirror the scores in matrix
    """
    interactions = []
    for score_file in score_files:
        interactions.extend(parse_scores(score_file))
    matrix = interaction_matrix(interactions, unique=unique)
    matrix.to_csv(output_file, sep="\t", index_label="Target")


def interaction_matrix(interactions: list[Interaction], unique: bool = False) -> pandas.DataFrame:
    """
    Converts list of interactions into a DataFrame.

    :param unique: keep only the maximum value of bait-target and target-bait for both
    :param interactions: interactions
    :return: DataFrame containing one line per target and one column per bait
    """
    if unique:
        interactions_inverse = [Interaction(i.target, i.bait, i.score) for i in interactions]
        interactions = interactions + interactions_inverse
    interactions_df = pandas.DataFrame({"Bait": [interaction.bait for interaction in interactions],
                                        "Target": [interaction.target for interaction in interactions],
                                        "Score": [interaction.score for interaction in interactions]})
    matrix = interactions_df.pivot_table(index="Target", columns="Bait", values="Score", aggfunc="max")
    return matrix


def parse_scores(score_file: TextIO) -> list[Interaction]:
    """
    Parses score file.

    :param score_file: tab delimited files containing bait, target and score
    :return: list of interactions
    """
    interactions = []
    score_file.readline()  # Skip header
    for line in score_file:
        if line.startswith('#'):
            continue
        columns = line.rstrip('\r\n').split('\t')
        bait = columns[0]
        target = columns[1]
        score = float(columns[2])
        interactions.append(Interaction(bait, target, score))
    return interactions


if __name__ == '__main__':
    main()
