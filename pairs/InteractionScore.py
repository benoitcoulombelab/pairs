import argparse
import math
import sys
from typing import TextIO

import smokesignal
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import molecular_weight, seq1

MISSING_CHAIN_EVENT = "missing_chain"
HYDROPHOBIC_AMINO_ACIDS = {"ALA", "CYS", "GLY", "ILE", "LEU", "MET", "PHE",
                           "PRO", "TRP", "TYR", "VAL"}
POLAR_AMINO_ACIDS = {"ARG": "D", "ASN": "AD", "ASP": "A", "GLN": "AD",
                     "GLU": "A", "HIS": "AD", "LYS": "D", "SER": "AD",
                     "THR": "AD", "TRP": "D", "TYR": "AD"}
CHARGED_AMINO_ACIDS = {"ARG": "+", "ASP": "-", "GLU": "-", "HIS": "+",
                       "LYS": "+"}


def main(argv: list[str] = None):
  parser = argparse.ArgumentParser(
      description="Compute protein-protein interaction score from PDB file")
  parser.add_argument('input', nargs='?', type=argparse.FileType('r'),
                      default=sys.stdin,
                      help="PDB or CIF file created by AlphaFold containing two (or more) proteins of interest")
  parser.add_argument('-a', '--first', default="A",
                      help="Chains of first protein separated by ','  (default: %(default)s)")
  parser.add_argument('-b', '--second', default="B",
                      help="Chains of second protein separated by ','  (default: %(default)s)")
  parser.add_argument('-c', '--count', action="store_true", default=False,
                      help="Score is the number of residues pairs at a distance less than radius parameter")
  parser.add_argument('-r', '--radius', type=float, default=6.0,
                      help="Distance between atoms from different residues to assume interaction "
                           "(default: %(default)s)")
  parser.add_argument('-w', '--weight', action="store_true", default=False,
                      help="Normalize count by protein pair weight - "
                           "'score / log2(sum weight of both proteins)'")
  parser.add_argument('-R', '--residues', type=argparse.FileType('w'),
                      metavar="RESIDUES",
                      help="Save residues pairs in tab separated file %(metavar)s")
  parser.add_argument('-A', '--atoms', type=argparse.FileType('w'),
                      metavar="ATOMS",
                      help="Save atoms pairs in tab separated file %(metavar)s")
  parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout,
                      help="Output file where to write PPI score")
  parser.add_argument('-P', '--partial', action="store_true", default=False,
                      help="Do not warn if all chains in PDB are not used for computing score")

  args = parser.parse_args(argv)
  args.first = args.first.split(',')
  args.second = args.second.split(',')
  smokesignal.on(MISSING_CHAIN_EVENT, warn_missing_chain, max_calls=1)

  score = interaction_score(structure_file=args.input, radius=args.radius,
                            weight=args.weight, count=args.count,
                            first_chains=args.first, second_chains=args.second,
                            residues=args.residues,
                            atoms=args.atoms, partial=args.partial)
  args.output.write(f"{score}\n")


def interaction_score(structure_file: TextIO, radius: float = 6,
    weight: bool = False, count: bool = False,
    first_chains: list[str] = ["A"],
    second_chains: list[str] = ["B"],
    residues: TextIO = None,
    atoms: TextIO = None,
    partial: bool = False) -> float:
  """
  Compute score of protein-protein interaction from PDB file.

  By default, the score is the number of residue pairs that have atoms at a distance below radius.

  :param structure_file: PDB or CIF file
  :param radius: maximal distance between two residues' atoms to consider that the two residues interact
  :param weight: if True, normalize score by proteins' weight
  :param count: if True, score is the number of residue pairs below radius
  :param first_chains: chains of the first protein
  :param second_chains: chains of the second protein
  :param residues: write residues pairs to this file, if specified
  :param atoms: write atoms pairs to this file, if specified
  :param partial: do not warn if all chains in PDB are not used for computing score
  :return: score of proteins interaction
  """
  parser = PDBParser()
  if structure_file.name.endswith(".cif"):
    parser = MMCIFParser()
  structure = parser.get_structure("unknown", structure_file)
  if not partial:
    for chain in structure[0]:
      if chain.get_id() not in first_chains + second_chains:
        smokesignal.emit(MISSING_CHAIN_EVENT, chain)
        break

  proteins_atoms = []
  for chain in first_chains + second_chains:
    proteins_atoms.extend(potential_interactor_atoms(structure[0][chain]))

  neighbor_search = NeighborSearch(proteins_atoms)
  interactions = search_interactions(neighbor_search=neighbor_search,
                                     radius=radius, level='R',
                                     first_chains=first_chains,
                                     second_chains=second_chains)

  if residues:
    write_residues(interactions, residues)

  if atoms:
    atoms_interactions = search_interactions(neighbor_search=neighbor_search,
                                             radius=radius, level='A',
                                             first_chains=first_chains,
                                             second_chains=second_chains)
    write_atoms(atoms_interactions, atoms)

  if count:
    score = len(interactions)
  else:
    score = 0
    for interaction in interactions:
      distance = minimal_distance(interaction[0], interaction[1])
      score = score + 10 / (2 ** distance)
  if weight:
    protein_pair_weight = 0
    for chain_name in first_chains + second_chains:
      protein_sequence = seq1(
          "".join(
              residue.get_resname() for residue in structure[0][chain_name]))
      protein_pair_weight = protein_pair_weight + molecular_weight(
          protein_sequence, seq_type="protein")
    return score / math.log2(protein_pair_weight)
  else:
    return score


def minimal_distance(residue_a: Residue, residue_b: Residue) -> float:
  """
  Returns distance between the closest potential interactor atoms of both residues.

  :param residue_a: first residue
  :param residue_b: second residue
  :return: distance between the closest potential interactor atoms of both residues
  """
  atoms_a = potential_interactor_atoms(residue_a)
  atoms_b = potential_interactor_atoms(residue_b)
  min_distance = atoms_a[0] - atoms_b[0]
  for atom_a in atoms_a:
    for atom_b in atoms_b:
      min_distance = min(atom_a - atom_b, min_distance)
  return min_distance


def potential_interactor_atoms(entity: Chain | Residue) -> list[Atom]:
  """
  Returns list of atoms that can interact with other residues' atoms.

  The atoms that can interact are considered as being all atoms except for:
    * Atoms in the peptide bound (N, CA, C, O, OXT)
    * Hydrogen

  :param entity: chain or residue from PDB file
  :return: list of atoms that can interact with other residues' atoms
  """
  atoms = []
  residues = [entity] if isinstance(entity, Residue) else entity.get_residues()
  for residue in residues:
    for atom in residue:
      if atom.get_name() not in ["N", "CA", "C", "O",
                                 "OXT"] and not atom.get_name().startswith("H"):
        atoms.append(atom)
  return atoms


def search_interactions(neighbor_search: NeighborSearch, radius: float,
    level: str = 'R',
    first_chains: list[str] = ["A"],
    second_chains: list[str] = ["B"],
) -> list[(Residue, Residue)] | list[(Atom, Atom)]:
  """
  Search for residue/atom pairs that are interacting.

  :param neighbor_search: instance of NeighborSearch containing atoms
  :param radius: maximum distance between atoms
  :param level: 'R' to search residues, 'A' to search atoms
  :param first_chains: chains of the first protein
  :param second_chains: chains of the second protein
  :return: list of residue/atom pairs that are interacting
  """
  interactions = neighbor_search.search_all(radius, level)
  interactions = [residue_pair for residue_pair in interactions
                  if (get_chain(residue_pair[0]).get_id() in first_chains
                      and get_chain(residue_pair[1]).get_id() in second_chains)
                  or (get_chain(residue_pair[0]).get_id() in second_chains
                      and get_chain(residue_pair[1]).get_id() in first_chains)]
  return interactions


def get_chain(entity: Residue | Atom) -> Chain:
  """
  Returns residue's/atom's chain.

  :param entity: residue or chain
  :param level: 'R' if entity is a residue, 'A' if entity is an atom
  :return: residue's/atom's chain
  """
  if isinstance(entity, Residue):
    return entity.get_parent()
  elif isinstance(entity, Atom):
    return entity.get_parent().get_parent()
  else:
    raise AssertionError(f"PDB entity '{entity}' is not one of Residue or Atom")


def write_residues(residue_pairs: list[(Residue, Residue)],
    output_file: TextIO):
  """
  Write residue pairs to output file

  :param residue_pairs: residues pairs
  :param output_file: output file
  """
  output_file.write(
      "Chain A\tResidue number A\tResidue name A\t"
      "Chain B\tResidue number B\tResidue name B\t"
      "Distance\tBond type (guess)\n")
  for residue_pair in residue_pairs:
    for i in range(0, 2):
      residue = residue_pair[i]
      chain = residue.get_parent()
      output_file.write(f"{chain.get_id()}")
      output_file.write("\t")
      output_file.write(f"{residue.get_id()[1]}")
      output_file.write("\t")
      output_file.write(f"{residue.get_resname()}")
      output_file.write("\t")
    distance = minimal_distance(residue_pair[0], residue_pair[1])
    output_file.write(f"{distance}")
    output_file.write("\t")
    if is_charged_bond(residue_pair[0], residue_pair[1]):
      output_file.write("Charged")
    elif is_hydrophobic_bond(residue_pair[0], residue_pair[1]):
      output_file.write("Hydrophobic")
    elif is_hydrogen_bond(residue_pair[0], residue_pair[1]):
      output_file.write("Hydrogen")
    output_file.write("\n")


def write_atoms(atom_pairs: list[(Atom, Atom)], output_file: TextIO):
  """
  Write atom pairs to output file

  :param atom_pairs: atom pairs
  :param output_file: output file
  """
  output_file.write("Chain A\tResidue number A\tResidue name A\tAtom A\t"
                    "Chain B\tResidue number B\tResidue name B\tAtom B\t"
                    "Distance\tBond type (guess)\n")
  for atom_pair in atom_pairs:
    for i in range(0, 2):
      atom = atom_pair[i]
      residue = atom.get_parent()
      chain = residue.get_parent()
      output_file.write(f"{chain.get_id()}")
      output_file.write("\t")
      output_file.write(f"{residue.get_id()[1]}")
      output_file.write("\t")
      output_file.write(f"{residue.get_resname()}")
      output_file.write("\t")
      output_file.write(f"{atom.get_name()}")
      output_file.write("\t")
    distance = minimal_distance(atom_pair[0], atom_pair[1])
    output_file.write(f"{distance}")
    output_file.write("\t")
    if is_charged_bond(atom_pair[0].get_parent(), atom_pair[1].get_parent()):
      output_file.write("Charged")
    elif is_hydrophobic_bond(atom_pair[0].get_parent(),
                             atom_pair[1].get_parent()):
      output_file.write("Hydrophobic")
    elif is_hydrogen_bond(atom_pair[0].get_parent(), atom_pair[1].get_parent()):
      output_file.write("Hydrogen")
    output_file.write("\n")


def is_charged_bond(first_residue: Residue, second_residue: Residue) -> bool:
  """
  Returns True if first_residue and second_residue are likely in a charged bond, False otherwise.

  :param first_residue: residue from first protein
  :param second_residue: residue from second protein
  :return: true if first_residue and second_residue are likely in a charged bond, False otherwise
  """
  if first_residue.get_resname() in CHARGED_AMINO_ACIDS and second_residue.get_resname() in CHARGED_AMINO_ACIDS:
    return CHARGED_AMINO_ACIDS[first_residue.get_resname()] != \
      CHARGED_AMINO_ACIDS[second_residue.get_resname()]
  return False


def is_hydrophobic_bond(first_residue: Residue,
    second_residue: Residue) -> bool:
  """
  Returns True if first_residue and second_residue are likely in a hydrophobic bond, False otherwise.

  :param first_residue: residue from first protein
  :param second_residue: residue from second protein
  :return: true if first_residue and second_residue are likely in a hydrophobic bond, False otherwise
  """
  return is_hydrophobic_residue(first_residue) and is_hydrophobic_residue(
      second_residue)


def is_hydrophobic_residue(residue: Residue) -> bool:
  return residue.get_resname() in HYDROPHOBIC_AMINO_ACIDS


def is_hydrogen_bond(first_residue: Residue, second_residue: Residue) -> bool:
  """
  Returns True if first_residue and second_residue are likely in a hydrogen bond, False otherwise.

  :param first_residue: residue from first protein
  :param second_residue: residue from second protein
  :return: true if first_residue and second_residue are likely in a hydrogen bond, False otherwise
  """
  if first_residue.get_resname() in POLAR_AMINO_ACIDS and second_residue.get_resname() in POLAR_AMINO_ACIDS:
    if (("A" in POLAR_AMINO_ACIDS[first_residue.get_resname()]
         and "D" in POLAR_AMINO_ACIDS[second_residue.get_resname()])
        or ("D" in POLAR_AMINO_ACIDS[first_residue.get_resname()]
            and "A" in POLAR_AMINO_ACIDS[second_residue.get_resname()])):
      return True
  return False


def warn_missing_chain(chain: Chain):
  print(f"Chain {chain.get_id()} present in PDB but not used for scoring",
        file=sys.stderr)


if __name__ == '__main__':
  main()
