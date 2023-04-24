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
    gene_2_set = set()
    for infile in args.infile:
        for line in infile:
            if line.startswith('#') or line.startswith('Gene1'):
                continue
            columns = line.rstrip('\r\n').split('\t')
            gene_1 = mappings[columns[0]]
            gene_2 = mappings[columns[1]]
            gene_2_set.add(gene_2)
            if args.unique:
                gene_2_set.add(gene_1)
            count = int(columns[2])
            if gene_1 not in matrix:
                matrix[gene_1] = {}
            if args.unique and gene_2 not in matrix:
                matrix[gene_2] = {}
            if gene_2 not in matrix[gene_1]:
                matrix[gene_1][gene_2] = count
            if args.unique and gene_1 not in matrix[gene_2]:
                matrix[gene_2][gene_1] = count
            matrix[gene_1][gene_2] = max(matrix[gene_1][gene_2], count)
            if args.unique:
                matrix[gene_2][gene_1] = max(matrix[gene_2][gene_1], count)

    gene_1_list = list(matrix)
    gene_1_list.sort()
    gene_2_list = list(gene_2_set)
    gene_2_list.sort()

    args.output.write("Bait")
    for gene in gene_2_list:
        args.output.write(f"\t{gene}")
    args.output.write('\n')
    for gene_1 in gene_1_list:
        args.output.write(f"{gene_1}")
        for gene_2 in gene_2_list:
            count = matrix[gene_1][gene_2] if gene_2 in matrix[gene_1] else None
            args.output.write(f"\t{count if count else ''}")
        args.output.write('\n')


if __name__ == '__main__':
    main()
