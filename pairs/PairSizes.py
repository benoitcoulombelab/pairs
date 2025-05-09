import argparse
import os
from typing import TextIO

from Bio import SeqIO, SeqRecord

from pairs import FastaId


def dir_path(string: str):
  if not string or os.path.isdir(string):
    return string
  else:
    raise NotADirectoryError(string)


def main(argv: list[str] = None):
  parser = argparse.ArgumentParser(
      description="Creates a text file containing pair sizes.")
  parser.add_argument('-b', '--baits', type=argparse.FileType('r'),
                      required=True,
                      help="FASTA file containing baits")
  parser.add_argument('-t', '--targets', type=argparse.FileType('r'),
                      required=True,
                      help="FASTA file containing targets")
  parser.add_argument('-u', '--unique', action='store_true',
                      help="Save only one JSON file per pair "
                           "- do not save POLR2B-POLR2A pair if POLR2A-POLR2B is also present.")
  parser.add_argument('-i', '--identity', action='store_true',
                      help="Don't save JSON file of a protein with itself "
                           "- same protein is present in both baits and targets.")
  parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default="pair_sizes.txt",
                      help="Output file.  (default: %(default)s)")

  args = parser.parse_args(argv)

  pair_sizes(baits=args.baits, targets=args.targets, output=args.output,
             unique=args.unique, skip_identity=args.identity)


def pair_sizes(baits: TextIO, targets: TextIO, output: TextIO,
    unique: bool = False,
    skip_identity: bool = False):
  """
  Creates a text file containing pair sizes.

  :param baits: baits
  :param targets: targets
  :param output: output file
  :param unique: save only one JSON file per unique pair -
                 do not save POLR2B-POLR2A pair if POLR2A-POLR2B is also present
  :param skip_identity: don't save JSON file of a protein with itself -
                        if the same protein is present in both baits and targets
  """
  baits = parse_fasta(baits)
  targets = parse_fasta(targets)

  processed_ids = set()
  for bait in baits:
    for target in targets:
      if skip_identity and bait == target:
        continue
      reverse_id = f"{target}__{bait}"
      if unique and reverse_id in processed_ids:
        continue
      merge_id = f"{bait}__{target}"
      size = len(str(baits[bait].seq)) + len(str(targets[target].seq))
      output.write(f"{merge_id}\t{str(size)}\n")
      processed_ids.add(merge_id)


def parse_fasta(fasta: TextIO) -> dict[str, SeqRecord]:
  """
  Parses FASTA and returns all sequences found in file mapped by ID.

  The ID of each sequence is found using :func:`FastaId.fasta_id`

  :param fasta: FASTA file
  :return: all sequences found in file mapped by ID
  """
  sequences = {}
  for record in SeqIO.parse(fasta, "fasta"):
    seq_id = FastaId.fasta_id(record.description)
    sequences[seq_id] = record
  return sequences


if __name__ == '__main__':
  main()
