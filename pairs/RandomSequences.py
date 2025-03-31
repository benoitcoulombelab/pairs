import argparse
import random
import re
import sys
from typing import TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight

VALID_AMINO_ACID_REGEX = r"^[ACDEFGHIKLMNPQRSTVWY]+$"


def main(argv: list[str] = None):
  parser = argparse.ArgumentParser(
      description="Generates random sequences with same lengths and similar distribution of AA as FASTA input.")
  parser.add_argument('fasta', type=argparse.FileType('r'), help="FASTA file")
  parser.add_argument('-s', '--sequence', action='store_true',
                      help="Ignore sequences with an invalid amino acid")
  parser.add_argument('output', nargs="?", type=argparse.FileType('w'),
                      default=sys.stdout,
                      help="FASTA output containing random sequences")

  args = parser.parse_args(argv)

  random_sequences(fasta=args.fasta, output=args.output,
                   invalid_sequence=args.sequence)


def random_sequences(fasta: TextIO, output: TextIO,
    invalid_sequence: bool = False):
  """
  Generates random sequences with similar distribution of AA as FASTA input.

  The name of the sequences will be "DECOY_{incremented}" where "{incremented}" an incremented number.

  :param fasta: FASTA file
  :param output: FASTA output containing random sequences
  :param invalid_sequence: Ignore sequences with an invalid amino acid
  """
  input_sequences = parse_fasta(fasta)
  if invalid_sequence:
    input_sequences = [seq for seq in input_sequences if
                       re.match(VALID_AMINO_ACID_REGEX, str(seq.seq))]
  counts = aa_count(input_sequences)
  aa_sum = sum(counts.values())
  aa_probabilities = {aa: counts[aa] / aa_sum for aa in counts}
  output_sequences = []
  incremented = 1
  for seq in input_sequences:
    out_seq = generate_sequence(len(seq.seq), aa_probabilities)
    mass = int(molecular_weight(out_seq, seq_type="protein"))
    name = f"DECOY_{incremented}"
    incremented = incremented + 1
    output_sequences.append(
        SeqRecord(out_seq, id=name, description=f"MASS={mass}"))
  SeqIO.write(output_sequences, output, "fasta")


def aa_count(fasta: list[SeqIO.SeqRecord]) -> dict[str, float]:
  """
  Returns number of copies of each amino acids in FASTA file.

  :param fasta: list of sequences from a FASTA file
  :return: number of copies of each amino acids in FASTA file
  """
  counts = {}
  for seq in fasta:
    for aa in seq.seq:
      counts[aa] = counts[aa] + 1 if aa in counts else 1
  return counts


def parse_fasta(fasta: str) -> list[SeqIO.SeqRecord]:
  """
  Returns sequences read from FASTA file.

  :param fasta: FASTA file
  :return: sequences read from FASTA file
  """
  sequences = []
  for record in SeqIO.parse(fasta, "fasta"):
    sequences.append(record)
  return sequences


def generate_sequence(length: int, aa_probabilities: dict[str, float]) -> Seq:
  """
  Returns a random protein sequence of length equals to length.

  :param length: random value between 0.0 to 1.0
  :param aa_probabilities: probability of each amino acid, total of values must equal 1.0
  :return: random protein sequence of specified length
  """
  return Seq("".join(random.choices(list(aa_probabilities.keys()),
                                    list(aa_probabilities.values()), k=length)))


if __name__ == '__main__':
  main()
