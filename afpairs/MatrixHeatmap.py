import argparse
import math
import statistics
import sys


class Bait:
    def __init__(self, name: str, weight: float):
        self.name = name
        self.weight = weight


def main():
    parser = argparse.ArgumentParser(description="Fix count matrix values for heatmap.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Matrix counts file")
    parser.add_argument('-z', '--zscore', action="store_true", default=False,
                        help="Apply Z-score to counts "
                             " (default: %(default)s)")
    parser.add_argument('-u', '--unique', action="store_true", default=False,
                        help="Protein pairs are unique, so mirror the scores in matrix "
                             "(default: %(default)s)")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing scores")

    args = parser.parse_args()

    headers = args.infile.readline()
    headers = headers.rstrip('\r\n').split('\t')
    weight = headers[1] == "Weight"
    target_start = 2 if weight else 1
    targets = headers[target_start:]

    baits = {}
    scores = {}
    for line in args.infile:
        columns = line.rstrip('\r\n').split('\t')
        bait = Bait(columns[0], float(columns[1]) if weight else 0)
        baits[bait.name] = bait
        scores[bait.name] = [float(value if value else 'nan') for value in columns[target_start:]]

    if args.zscore:
        for bait in scores:
            bait_scores = [score for score in scores[bait] if not math.isnan(score)]
            bait_mean = statistics.fmean(bait_scores)
            bait_std = statistics.stdev(bait_scores)
            scores[bait] = [(value - bait_mean) / bait_std for value in scores[bait]]

    if args.unique:
        for bait in scores:
            for target_index in range(0, len(targets)):
                target = targets[target_index]
                if bait in targets and target in scores:
                    bait_index = targets.index(bait)
                    interaction_max = max(scores[bait][target_index], scores[target][bait_index])
                    scores[bait][target_index] = interaction_max
                    scores[target][bait_index] = interaction_max

    args.outfile.write("Bait")
    if weight:
        args.outfile.write("\tWeight")
    args.outfile.write('\t')
    args.outfile.write('\t'.join(headers[target_start:]))
    args.outfile.write('\n')
    for bait in scores:
        args.outfile.write(f"{bait}\t")
        if weight:
            args.outfile.write(f"{baits[bait].weight}\t")
        args.outfile.write('\t'.join([str(score if not math.isnan(score) else '') for score in scores[bait]]))
        args.outfile.write('\n')


if __name__ == '__main__':
    main()
