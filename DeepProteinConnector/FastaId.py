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

    for line in args.infile:
        if line.startswith('>'):
            sequence_name = line

    if not sequence_name:
        sys.exit(1, "Could not find sequence name starting with '>'")

    sequence_id = fasta_id(sequence_name, args.gene)
    if id:
        print(sequence_id)
        sys.exit(0)
    else:
        sys.exit(1, "Could not find any id in sequence name")


def fasta_id(sequence_name, gene=False):
    """
    Returns sequence accession/id from the sequence name.

    When possible, an accession number will be returned.
    Accession number is assumed to be the third element like P24928 in sp|P24928|RPB1_HUMAN.

    If there is no third element, the second is used, otherwise the first is used.

    :param sequence_name: Sequence name - the complete line of FASTA starting with '>'
    :param gene: If True, will return the gene name if found
    :return: Sequence id from the sequence name
    """
    patterns = [r"[\w]+\|[\w]+\|([\w]+)", r"[\w]+\|([\w]+)", r"([\w]+)"]
    if gene:
        patterns.insert(0, r"GN=([\w]+)")
    for pattern in patterns:
        re_search = re.search(pattern, sequence_name)
        if re_search:
            return f"{re_search[1]}"
    return None


if __name__ == '__main__':
    main()
