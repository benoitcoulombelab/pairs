import argparse
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Renames sequence of FASTA to keep gene name and protein name.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA file")
    parser.add_argument('-g', '--gene', default=".*",
                        help="Gene regex (default: %(default)s)")
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file")

    args = parser.parse_args()

    for line in args.infile:
        if not line.startswith('>'):
            args.outfile.write(line)
        else:
            search_gene = re.search(r"GN=(\w+)", line)
            if search_gene:
                gene = search_gene[1]
            search_protein = re.search(r"^>\w+\|\w+\|(\w+)", line)
            if search_protein:
                protein = search_protein[1]
            if protein:
                args.outfile.write(">")
                if gene:
                    args.outfile.write(f"{gene}_")
                args.outfile.write(protein)
                args.outfile.write("\n")
            else:
                args.outfile.write(line)


if __name__ == '__main__':
    main()
