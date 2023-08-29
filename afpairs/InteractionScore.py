import argparse
import math
import sys
from typing import TextIO

from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import molecular_weight, seq1


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(description="Counts the number of potentials interactions between "
                                                 "residues of two proteins.")
    parser.add_argument('input', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="PDB file created by AlphaFold containing two proteins of interest")
    parser.add_argument('-a', '--first', default="A",
                        help="Chains of first protein separated by ','  (default: %(default)s)")
    parser.add_argument('-b', '--second', default="B",
                        help="Chains of second protein separated by ','  (default: %(default)s)")
    parser.add_argument('-r', '--radius', type=float, default=6.0,
                        help="Distance between atoms from different residues to assume interaction "
                             "(default: %(default)s)")
    parser.add_argument('-w', '--weight', action="store_true", default=False,
                        help="Normalize count by protein pair weight - "
                             "'count / log2(sum weight of both proteins)'")
    parser.add_argument('-R', '--residues', type=argparse.FileType('w'), metavar="RESIDUES",
                        help="Save pairs of residues in tab separated file %(metavar)s")
    parser.add_argument('-A', '--atoms', type=argparse.FileType('w'), metavar="ATOMS",
                        help="Save pairs of atoms in tab separated file %(metavar)s")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file where to write number of potentials interactions")

    args = parser.parse_args(argv)
    args.first = args.first.split(',')
    args.second = args.second.split(',')

    score = interaction_score(pdb=args.input, radius=args.radius, weight=args.weight,
                              first_chains=args.first, second_chains=args.second, residues=args.residues,
                              atoms=args.atoms)
    args.output.write(f"{score}\n")


def interaction_score(pdb: TextIO, radius: float = 6,
                      weight: bool = False,
                      first_chains: list[str] = ["A"],
                      second_chains: list[str] = ["B"],
                      residues: TextIO = None,
                      atoms: TextIO = None) -> float:
    """
    Compute score of proteins interaction from PDB file.

    By default, the score is the number of residue pairs that have atoms at a distance below radius.

    :param pdb: PDB file
    :param radius: maximal distance between two residues' atoms to consider that the two residues interact
    :param weight: if True, normalize score by proteins' weight
    :param first_chains: chains of the first protein
    :param second_chains: chains of the second protein
    :param residues: write residues pairs to this file, if specified
    :param atoms: write atoms pairs to this file, if specified
    :return: score of proteins interaction
    """
    parser = PDBParser()
    structure = parser.get_structure("unknown", pdb)

    proteins_atoms = []
    for chain in first_chains + second_chains:
        proteins_atoms.extend(potential_interactor_atoms(structure[0][chain]))

    neighbor_search = NeighborSearch(proteins_atoms)
    interactions = search_interactions(neighbor_search=neighbor_search, radius=radius, level='R',
                                       first_chains=first_chains, second_chains=second_chains)

    if residues:
        write_residues(interactions, residues)

    if atoms:
        interactions = search_interactions(neighbor_search=neighbor_search, radius=radius, level='A',
                                           first_chains=first_chains, second_chains=second_chains)
        write_atoms(interactions, atoms)

    if weight:
        protein_pair_weight = 0
        for chain_name in first_chains + second_chains:
            protein_sequence = seq1("".join(residue.get_resname() for residue in structure[0][chain_name]))
            protein_pair_weight = protein_pair_weight + molecular_weight(protein_sequence, seq_type="protein")
        return len(interactions) / math.log2(protein_pair_weight)
    else:
        return len(interactions)


def potential_interactor_atoms(chain: Chain) -> list[Atom]:
    """
    Returns list of atoms that can interact with other residues' atoms.

    The atoms that can interact are considered as being all atoms except for:
      * Atoms in the peptide bound (N, CA, C, O, OXT)
      * Hydrogen

    :param chain: chain from PDB file
    :return: list of atoms that can interact with other residues' atoms
    """
    atoms = []
    for residue in chain:
        for atom in residue:
            if atom.get_name() not in ["N", "CA", "C", "O", "OXT"] and not atom.get_name().startswith("H"):
                atoms.append(atom)
    return atoms


def search_interactions(neighbor_search: NeighborSearch, radius: float, level: str = 'R',
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
                    if (get_chain(residue_pair[0], level).get_id() in first_chains
                        and get_chain(residue_pair[1], level).get_id() in second_chains)
                    or (get_chain(residue_pair[0], level).get_id() in second_chains
                        and get_chain(residue_pair[1], level).get_id() in first_chains)]
    return interactions


def get_chain(entity: Residue | Atom, level: str = 'R') -> Chain:
    """
    Returns residue's/atom's chain.

    :param entity: residue or chain
    :param level: 'R' if entity is a residue, 'A' if entity is an atom
    :return: residue's/atom's chain
    """
    if level == 'R':
        return entity.get_parent()
    elif level == 'A':
        return entity.get_parent().get_parent()
    else:
        raise AssertionError(f"Level '{level}' is not one of 'A' or 'R'")


def write_residues(residue_pairs: list[(Residue, Residue)], output_file: TextIO):
    """
    Write residue pairs to output file

    :param residue_pairs: residues pairs
    :param output_file: output file
    """
    output_file.write("Chain A\tResidue number A\tResidue name A\tChain B\tResidue number B\tResidue name B\n")
    for residue_pair in residue_pairs:
        for i in range(0, 2):
            residue = residue_pair[i]
            chain = residue.get_parent()
            output_file.write(f"{chain.get_id()}")
            output_file.write("\t")
            output_file.write(f"{residue.get_id()[1]}")
            output_file.write("\t")
            output_file.write(f"{residue.get_resname()}")
            output_file.write("\t" if i < 1 else "\n")


def write_atoms(atom_pairs: list[(Atom, Atom)], output_file: TextIO):
    """
    Write atom pairs to output file

    :param atom_pairs: atom pairs
    :param output_file: output file
    """
    output_file.write("Chain A\tResidue number A\tResidue name A\tAtom A\t"
                      "Chain B\tResidue number B\tResidue name B\tAtom B\n")
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
            output_file.write("\t" if i < 1 else "\n")


if __name__ == '__main__':
    main()
