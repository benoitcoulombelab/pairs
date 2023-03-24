import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Merge fasta files.")
    parser.add_argument('infiles', nargs='*', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA files to merge")
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file")

    args = parser.parse_args()

    for infile in args.infiles:
        for line in infile:
            args.outfile.write(line)


if __name__ == '__main__':
    main()
