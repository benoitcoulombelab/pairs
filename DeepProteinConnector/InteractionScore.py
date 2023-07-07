import argparse
import math
import os
import sys

from Bio.PDB import PDBParser, NeighborSearch
from Bio.SeqUtils import molecular_weight, seq1


def main():
    parser = argparse.ArgumentParser(description="Counts the number of potentials interactions between "
                                                 "residues of two proteins.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="PDB file created by AlphaFold containing two proteins of interest")
    parser.add_argument('-n', '--name',
                        help="Name of structure in PDB file (default: PDB filename or 'Unknown' if standard input)")
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

    args = parser.parse_args()
    if not args.name:
        if args.infile == sys.stdin:
            args.name = "Unknown"
        else:
            args.name = os.path.basename(os.path.splitext(args.infile.name)[0])

    score = interaction_score(args.infile, args.name, args.radius, args.weight, args.residues, args.atoms)
    args.output.write(f"{score}\n")


def interaction_score(pdb: "PDB file", name: "structure name" = "Unknown", radius: float = 6,
                      weight: "Normalize count by protein weight" = False,
                      residues: "File where to save pairs of residues" = None,
                      atoms: "File where to save pairs of atoms" = None) -> float:
    """Count number of interactions of PDB file"""
    parser = PDBParser()
    structure = parser.get_structure(name, pdb)

    atoms_a = potential_interactor_atoms(structure[0]["A"])
    atoms_b = potential_interactor_atoms(structure[0]["B"])

    neighbor_search = NeighborSearch(atoms_a + atoms_b)
    interactions = search_interactions(neighbor_search, radius, level='R')

    if residues:
        write_residues(interactions, residues)

    if atoms:
        interactions = search_interactions(neighbor_search, radius, level='A')
        write_atoms(interactions, atoms)

    if weight:
        protein_pair_weight = 0
        for chain_name in ["A", "B"]:
            protein_sequence = seq1("".join(residue.get_resname() for residue in structure[0][chain_name]))
            protein_pair_weight = protein_pair_weight + molecular_weight(protein_sequence, seq_type="protein")
        return len(interactions) / math.log2(protein_pair_weight)
    else:
        return len(interactions)


def potential_interactor_atoms(chain: "chain from PDB file") -> "list of atoms that can interact with other atoms":
    """Returns list of atoms that can interact with other atoms"""
    atoms = []
    for residue in chain:
        for atom in residue:
            if atom.get_name() not in ["N", "CA", "C", "O", "OXT"] and not atom.get_name().startswith("H"):
                atoms.append(atom)
    return atoms


def search_interactions(neighbor_search: NeighborSearch, radius: float, level: str) -> "list of interactions":
    """Search for interactions"""
    interactions = neighbor_search.search_all(radius, level)
    interactions = [residue_pair for residue_pair in interactions
                    if get_chain(residue_pair[0], level) != get_chain(residue_pair[1], level)]
    return interactions


def get_chain(entity: "PDB entity", level: "Level of entity - 'A' for atom or 'R' for residue" = 'R') -> \
        "entity's chain":
    if level == 'R':
        return entity.get_parent()
    elif level == 'A':
        return entity.get_parent().get_parent()
    else:
        raise AssertionError(f"Level '{level}' is not one of 'A' or 'R'")


def write_residues(interactions: "list of residue interactions", outfile: "output file"):
    """Write residues that have potential interactions to output file"""
    outfile.write("Chain A\tResidue number A\tResidue name A\tChain B\tResidue number B\tResidue name B\n")
    for residue_pair in interactions:
        for i in range(0, 2):
            residue = residue_pair[i]
            chain = residue.get_parent()
            outfile.write(f"{chain.get_id()}")
            outfile.write("\t")
            outfile.write(f"{residue.get_id()[1]}")
            outfile.write("\t")
            outfile.write(f"{residue.get_resname()}")
            outfile.write("\t" if i < 1 else "\n")


def write_atoms(interactions: "list of residue interactions", outfile: "output file"):
    """Write atoms that have potential interactions to output file"""
    outfile.write("Chain A\tResidue number A\tResidue name A\tAtom A\t"
                  "Chain B\tResidue number B\tResidue name B\tAtom B\n")
    for atom_pair in interactions:
        for i in range(0, 2):
            atom = atom_pair[i]
            residue = atom.get_parent()
            chain = residue.get_parent()
            outfile.write(f"{chain.get_id()}")
            outfile.write("\t")
            outfile.write(f"{residue.get_id()[1]}")
            outfile.write("\t")
            outfile.write(f"{residue.get_resname()}")
            outfile.write("\t")
            outfile.write(f"{atom.get_name()}")
            outfile.write("\t")
            outfile.write("\t" if i < 1 else "\n")


if __name__ == '__main__':
    main()
