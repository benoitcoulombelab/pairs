import argparse
import sys
import os

from Bio.PDB import PDBParser, NeighborSearch


def potential_interactor_atoms(chain):
    atoms = []
    for residue in chain:
        for atom in residue:
            if atom.get_name() not in ["N", "CA", "C", "O", "OXT"] and not atom.get_name().startswith("H"):
                atoms.append(atom)
    return atoms


def main():
    parser = argparse.ArgumentParser(description="Counts the number of potentials interactions between " +
                                                 "residues of two proteins.")
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="PDB file created by AlphaFold containing two proteins of interest")
    parser.add_argument('-n', '--name',
                        help="Name of structure in PDB file (default: PDB filename or 'Unknown' if standard input)")
    parser.add_argument('-r', '--radius', type=float, default=6.0,
                        help="Distance between atoms from different residues to assume interaction (default: %(default)s)")
    parser.add_argument('-R', '--residues', type=argparse.FileType('w'), metavar="RESIDUES",
                        help="Save pairs of residues in tab separated file %(metavar)s")
    parser.add_argument('-A', '--atoms', type=argparse.FileType('w'), metavar="ATOMS",
                        help="Save pairs of atoms in tab separated file %(metavar)s")

    args = parser.parse_args()
    if not args.name:
        if args.infile == sys.stdin:
            args.name = "Unknown"
        else:
            args.name = os.path.basename(os.path.splitext(args.infile.name)[0])

    parser = PDBParser()
    structure = parser.get_structure(args.name, args.infile)

    atoms_a = potential_interactor_atoms(structure[0]["A"])
    atoms_b = potential_interactor_atoms(structure[0]["B"])

    interactions_search = NeighborSearch(atoms_a + atoms_b)
    interactions = interactions_search.search_all(args.radius, level='R')
    interactions = [residue_pair for residue_pair in interactions if residue_pair[0].get_parent() != residue_pair[1].get_parent()]

    print(f"{len(interactions)}")

    if args.residues:
        args.residues.write("Chain A\tResidue number A\tResidue name A\tChain B\tResidue number B\tResidue name B\n")
        for residue_pair in interactions:
            args.residues.write(f"{residue_pair[0].get_parent().get_id()}")
            args.residues.write("\t")
            args.residues.write(f"{residue_pair[0].get_id()[1]}")
            args.residues.write("\t")
            args.residues.write(f"{residue_pair[0].get_resname()}")
            args.residues.write("\t")
            args.residues.write(f"{residue_pair[1].get_parent().get_id()}")
            args.residues.write("\t")
            args.residues.write(f"{residue_pair[1].get_id()[1]}")
            args.residues.write("\t")
            args.residues.write(f"{residue_pair[1].get_resname()}")
            args.residues.write("\n")

    if args.atoms:
        interactions = interactions_search.search_all(args.radius, level='A')
        interactions = [atom_pair for atom_pair in interactions
                        if atom_pair[0].get_parent().get_parent() != atom_pair[1].get_parent().get_parent()]
        args.atoms.write("Chain A\tResidue number A\tResidue name A\tAtom A\t" +
                         "Chain B\tResidue number B\tResidue name B\tAtom B\n")
        for atom_pair in interactions:
            args.atoms.write(f"{atom_pair[0].get_parent().get_parent().get_id()}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[0].get_parent().get_id()[1]}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[0].get_parent().get_resname()}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[0].get_name()}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[1].get_parent().get_parent().get_id()}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[1].get_parent().get_id()[1]}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[1].get_parent().get_resname()}")
            args.atoms.write("\t")
            args.atoms.write(f"{atom_pair[1].get_name()}")
            args.atoms.write("\n")


if __name__ == '__main__':
    main()
