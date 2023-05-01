import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Generate count matrix from DeepProteinConnector counts.")
    parser.add_argument('infile', nargs='*', type=argparse.FileType('r'), default=sys.stdin,
                        help="Tab delimited text file containing protein id/accession and counts")
    parser.add_argument('-m', '--mapping', type=argparse.FileType('r'), default="mapping.txt",
                        help="Tab delimited text file containing source ids and converted ids  (default: %(default)s)")
    parser.add_argument('-s', '--source_column', type=int, default='1',
                        help="Column index of source ids in conversion file - 1 means first column of file" +
                             "   (default: %(default)s)")
    parser.add_argument('-c', '--converted_column', type=int, default='2',
                        help="Column index of converted ids in conversion file - 1 means first column of file" +
                             "   (default: %(default)s)")
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
        source = columns[args.source_column - 1]
        converted = columns[args.converted_column - 1]
        mappings[source] = converted

    matrix = {}
    target_set = set()
    for infile in args.infile:
        infile.readline()  # Skip header
        for line in infile:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\r\n').split('\t')
            bait = mappings[columns[0]]
            target = mappings[columns[1]]
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
    bait_list.sort()
    target_list = list(target_set)
    target_list.sort()

    args.output.write("Bait")
    for gene in target_list:
        args.output.write(f"\t{gene}")
    args.output.write('\n')
    for bait in bait_list:
        args.output.write(f"{bait}")
        for target in target_list:
            count = matrix[bait][target] if target in matrix[bait] else None
            args.output.write(f"\t{count if count is not None else ''}")
        args.output.write('\n')


if __name__ == '__main__':
    main()
