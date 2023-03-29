import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Prints id of FASTA file. \n\n"
                                     "If no option is used or if it fails to use the option, an accession number "
                                     "from the FASTA file will be printed.\n"
                                     "Accession is assumed to be the third element, when possible "
                                     "like P24928 in sp|P24928|RPB1_HUMAN.\n"
                                     "If there is no third element, the second is used, otherwise the first is used")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file")
    parser.add_argument('-g', '--gene', action='store_true',
                        help="Use gene name as FASTA id.")

    args = parser.parse_args()

    patterns = [r"[\w]+\|[\w]+\|([\w]+)", r"[\w]+\|([\w]+)", r"([\w]+)"]
    if args.gene:
        patterns.insert(0, r"GN=([\w]+)")
    for line in args.infile:
        if line.startswith('>'):
            for pattern in patterns:
                re_search = re.search(pattern, line)
                if re_search:
                    print(f"{re_search[1]}")
                    sys.exit(0)
    sys.exit(1, "Could not find any id")


if __name__ == '__main__':
    main()
