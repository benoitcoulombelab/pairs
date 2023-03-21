import argparse
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser(description="Splits FASTA in one file per sequence. "
                                                 "The output filenames correspond to the accession number "
                                                 "in the sequence name.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="FASTA containing multiple sequences")
    parser.add_argument('-o', '--outdir', default="",
                        help="Output directory")

    args = parser.parse_args()

    if args.outdir and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    lines = args.infile.readlines()
    for line_index in range(0, len(lines)):
        line = lines[line_index]
        line_index = line_index + 1
        if line.startswith('>'):
            name = re.search(r"^>\w+\|(\w+)", line)[1]
            if not name:
                print(f"{name} cannot be empty", file=sys.stderr)
                continue
            with open(os.path.join(args.outdir, name+".fasta"), 'w') as outfile:
                outfile.write(line)
                while line_index < len(lines) and not lines[line_index].startswith('>'):
                    line = lines[line_index]
                    line_index = line_index + 1
                    outfile.write(line)


if __name__ == '__main__':
    main()
