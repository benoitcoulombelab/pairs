import argparse
import re
import sys

from Bio import SeqIO


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Prints id of FASTA file. \n\n"
                                                 "If no option is used or if it fails to use the option, an accession "
                                                 "number from the FASTA file will be printed.\n"
                                                 "Accession is assumed to be the last element, when possible "
                                                 "like RPB1_HUMAN in sp|P24928|RPB1_HUMAN.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file")
    parser.add_argument('-g', '--gene', action='store_true',
                        help="Use gene name as FASTA id.")

    args = parser.parse_args(argv)

    try:
        sequence_name = next(SeqIO.parse(args.infile, "fasta")).description
    except StopIteration:
        raise AssertionError(f"Could not find any sequence in file {args.infile}")

    sequence_id = fasta_id(sequence_name, gene=args.gene)
    if id:
        print(sequence_id)
    else:
        raise AssertionError("Could not find any id in sequence name")


def fasta_id(sequence_name: str, gene: bool = False):
    """
    Returns sequence accession/id from the sequence name.

    When possible, an accession number will be returned.
    Accession number is assumed to be the last element like RPB1_HUMAN in sp|P24928|RPB1_HUMAN.

    :param sequence_name: Sequence name - the complete line of FASTA starting with '>'
    :param gene: If True, will return the gene name if found
    :return: Sequence id from the sequence name
    """
    if gene:
        pattern = r"GN=([\w]+)"
        re_search = re.search(pattern, sequence_name)
        if re_search:
            return f"{re_search[1]}"
    sequence_name = sequence_name.split(" ")[0]
    pattern = r">?\|?([^\|]+)\|?$"
    re_search = re.search(pattern, sequence_name)
    if re_search:
        sequence_name = re_search[1]
        sequence_name = re.sub(r"\W.*", "", sequence_name)
        return f"{sequence_name}"
    return None


if __name__ == '__main__':
    main()
