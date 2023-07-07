import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Keeps only proteins with gene regex.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file")
    parser.add_argument('-g', '--gene', default=".*",
                        help="Gene regex (default: %(default)s)")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file")

    args = parser.parse_args()

    lines = args.infile.readlines()
    for line_index in range(0, len(lines)):
        line = lines[line_index]
        line_index = line_index + 1
        if line.startswith('>'):
            search_gene = re.search("GN=" + args.gene, line)
            if search_gene:
                args.outfile.write(line)
                while line_index < len(lines) and not lines[line_index].startswith('>'):
                    line = lines[line_index]
                    line_index = line_index + 1
                    args.outfile.write(line)


if __name__ == '__main__':
    main()
