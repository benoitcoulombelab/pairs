import argparse
import sys
from typing import TextIO

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Extracts protein sequences from PDB file and output them as FASTA.")
    parser.add_argument('pdb', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="PDB file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Sequence of proteins found in PDB file in FASTA format")
    parser.add_argument('-p', '--prefix', default="",
                        help="Prefix to add to chain id as sequence name in FASTA")

    args = parser.parse_args(argv)

    pdb_fasta(pdb=args.pdb, output=args.output, prefix=args.prefix)


def pdb_fasta(pdb: TextIO, output: TextIO, prefix: str = ""):
    """
    Extracts protein sequences from PDB file and output them as FASTA.

    :param pdb: PDB file
    :param output: protein sequences in FASTA
    :param prefix: prefix to add to chain id as sequence name in FASTA
    """
    parser = PDBParser()
    structure = parser.get_structure("unknown", pdb)
    sequences = []
    for chain in structure[0]:
        sequence = seq1("".join(residue.get_resname() for residue in chain.get_residues()))
        sequences.append(SeqRecord(Seq(sequence), id=f"{prefix}{chain.get_id()}", description=""))
    SeqIO.write(sequences, output, "fasta")


if __name__ == '__main__':
    main()
