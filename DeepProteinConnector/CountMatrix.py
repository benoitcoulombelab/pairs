import argparse
import math
import sys


class Mapping:
    def __init__(self, source: str, converted: str, weight: float):
        self.source = source
        self.converted = converted
        self.weight = weight


def main():
    parser = argparse.ArgumentParser(description="Generate count matrix from DeepProteinConnector counts.")
    parser.add_argument('infile', nargs='*', type=argparse.FileType('r'), default=sys.stdin,
                        help="Tab delimited text file containing protein id/accession and counts")
    parser.add_argument('-m', '--mapping', type=argparse.FileType('r'), default="mapping.txt",
                        help="Tab delimited text file containing source ids and converted ids  (default: %(default)s)")
    parser.add_argument('-s', '--source', type=int, default='1',
                        help="Column index of source ids in conversion file - 1 means first column of file" +
                             "   (default: %(default)s)")
    parser.add_argument('-c', '--converted', type=int, default='2',
                        help="Column index of converted ids in conversion file - 1 means first column of file" +
                             "   (default: %(default)s)")
    parser.add_argument('-w', '--weight', type=int,
                        help="Column index of weights used for normalization - 1 means first column of file" +
                             "   (default: %(default)s)")
    parser.add_argument('-b', '--base', type=float, default=1.0,
                        help="Base value for normalization - values are multiplied by "
                             "'base_weight / log2(sum weight of both proteins)' " +
                             "  (default: %(default)s)")
    parser.add_argument('-u', '--unique', action='store_true',
                        help="Protein pairs are unique, so mirror the counts in matrix")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Tab delimited output file containing count matrix")

    args = parser.parse_args()

    mappings = {}
    for line in args.mapping:
        if line.startswith('#'):
            continue
        columns = line.rstrip('\r\n').split('\t')
        source = columns[args.source - 1]
        converted = columns[args.converted - 1]
        weight = float(columns[args.weight - 1]) if args.weight else None
        mappings[source] = Mapping(source, converted, weight)

    matrix = {}
    target_set = set()
    for infile in args.infile:
        infile.readline()  # Skip header
        for line in infile:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\r\n').split('\t')
            bait = columns[0]
            target = columns[1]
            target_set.add(target)
            if args.unique:
                target_set.add(bait)
            count = int(columns[2])
            if bait not in matrix:
                matrix[bait] = {}
            if args.unique and target not in matrix:
                matrix[target] = {}
            if target not in matrix[bait]:
                matrix[bait][target] = count
            if args.unique and bait not in matrix[target]:
                matrix[target][bait] = count
            matrix[bait][target] = max(matrix[bait][target], count)
            if args.unique:
                matrix[target][bait] = max(matrix[target][bait], count)

    bait_list = list(matrix)
    bait_list.sort(key=lambda m: mappings[m].converted)
    target_list = list(target_set)
    target_list.sort(key=lambda m: mappings[m].converted)

    args.output.write("Bait")
    for target in target_list:
        target_mapping = mappings[target]
        args.output.write(f"\t{target_mapping.converted}")
    args.output.write('\n')
    for bait in bait_list:
        bait_mapping = mappings[bait]
        args.output.write(f"{bait_mapping.converted}")
        for target in target_list:
            target_mapping = mappings[target]
            count = matrix[bait][target] if target in matrix[bait] else None
            if count and args.weight:
                count = count * args.base / math.log2(bait_mapping.weight + target_mapping.weight)
            args.output.write(f"\t{count if count is not None else ''}")
        args.output.write('\n')


if __name__ == '__main__':
    main()
