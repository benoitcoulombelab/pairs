import argparse
import re
import sys
from typing import TextIO, Tuple

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord


class Residue:
    def __init__(self, chain: str, index: int, name: str):
        self.chain = chain
        self.index = index
        self.name = name

    def __str__(self):
        return f"Residue({self.chain},{self.index},{self.name})"

    def __repr__(self):
        return f"Residue({self.chain},{self.index},{self.name})"


class ResiduePair:
    def __init__(self, first: Residue, second: Residue, bond: str = None):
        self.first = first
        self.second = second
        self.bond = bond

    def __str__(self):
        return f"ResiduePair({self.first},{self.second})"

    def __repr__(self):
        return f"ResiduePair({self.first},{self.second})"


class AlignmentRecord:
    def __init__(self, record: SeqRecord):
        self.positions = {}
        position = 0
        for index in range(0, len(record)):
            if record[index] != "-":
                self.positions[position] = index
                position = position + 1

    def position(self, residue_index: int) -> int:
        return self.positions[residue_index - 1] + 1 if self.positions else residue_index


class ConsensusResiduePair:
    def __init__(self, bait_residue_index: int, target_residue_index: int):
        self.bait_residue_index = bait_residue_index
        self.target_residue_index = target_residue_index
        self.residue_pairs = {}

    def append(self, bait: str, target: str, residue_pair: ResiduePair):
        self.residue_pairs[(bait, target)] = residue_pair
        return self

    def __str__(self):
        return f"({self.bait_residue_index},{self.target_residue_index}):" \
               f"[{';'.join([f'({bait},{target})' for bait, target in self.residue_pairs])}]"

    def __repr__(self):
        value = f"ConsensusResiduePair({self.bait_residue_index},{self.target_residue_index})"
        for (bait, target) in self.residue_pairs:
            residue_pair = self.residue_pairs[(bait, target)]
            value = value + f".append({bait}, {target}, {residue_pair})"
        return value


def main(argv: list[str] = None):
    parser = argparse.ArgumentParser(
        description="Extract consensus residue pairs from multiple protein pairs.")
    parser.add_argument('-r', '--residues', type=argparse.FileType('r'), nargs='+',
                        help="Residue-pairs files")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file")
    parser.add_argument('-n', '--name', default=r"([\w-]+)__([\w-]+)",
                        help="Regular expression to obtain protein/gene names based on residue-pairs filename "
                             " (default: %(default)s)")
    parser.add_argument('-c', '--consensus', type=float, default=0.5,
                        help="Minimal ratio to consider a consensus  (default: %(default)s)")
    parser.add_argument('-b', '--baits', type=argparse.FileType('r'),
                        help="Alignment file of baits - if None, baits are assumed to have the same sequence")
    parser.add_argument('-t', '--targets', type=argparse.FileType('r'),
                        help="Alignment file of targets - if None, targets are assumed to have the same sequence")
    parser.add_argument('-f', '--format', default="clustal",
                        help="Alignment file format - see https://biopython.org/wiki/AlignIO  (default: %(default)s)")

    args = parser.parse_args(argv)

    consensus_interface(residue_pair_files=args.residues, output_file=args.output, name=args.name,
                        consensus_ratio=args.consensus,
                        baits_file=args.baits, targets_file=args.targets, alignment_format=args.format)


def consensus_interface(residue_pair_files: TextIO, output_file: TextIO, name: str = r"([\w-]+)__([\w-]+)",
                        consensus_ratio: float = 0.5,
                        baits_file: TextIO = None, targets_file: TextIO = None, alignment_format: str = "clustal"):
    """
    Extract consensus residue pairs of multiple protein pairs.

    :param residue_pair_files: residue-pairs files
    :param output_file: output file
    :param name: regular expression to obtain protein/gene names based on residue-pairs filename
    :param consensus_ratio: minimal ratio to consider a consensus
    :param baits_file: alignment file of baits
    :param targets_file: alignment file of targets
    :param alignment_format: alignment file format
    """
    residue_pairs = {}
    for residue_pairs_file in residue_pair_files:
        re_match = re.search(name, residue_pairs_file.name)
        if not re_match:
            raise AssertionError(f"Expression {name} cannot be found in filename {residue_pairs_file}")
        bait, target = re_match.group(1, 2)
        residue_pairs[(bait, target)] = parse_residue_pairs(residue_pairs_file)

    baits = None
    if baits_file:
        baits = parse_alignment(baits_file, alignment_format)
        for bait, target in residue_pairs:
            if bait not in baits:
                raise AssertionError(f"Bait {bait} not in alignment {baits_file.name}")
    targets = None
    if targets_file:
        targets = parse_alignment(targets_file, alignment_format)
        for bait, target in residue_pairs:
            if target not in targets:
                raise AssertionError(f"Target {target} not in alignment {targets_file.name}")

    consensuses = consensus_residue_pairs(residue_pairs, baits, targets)

    output_file.write("Bait residue index\tTarget residue index\tConsensus count\t"
                      "Baits\tTargets\t"
                      "Chain A\tResidue number A\tResidue name A\t"
                      "Chain B\tResidue number B\tResidue name B\tBond type (guess)\n")
    for consensus in consensuses:
        if len(consensus.residue_pairs) >= len(residue_pairs) * consensus_ratio:
            output_file.write(f"{consensus.bait_residue_index}\t{consensus.target_residue_index}\t")
            output_file.write(f"{len(consensus.residue_pairs)}\t")
            output_file.write(
                f"{','.join([bait for bait, target in consensus.residue_pairs])}\t")
            output_file.write(
                f"{','.join([target for bait, target in consensus.residue_pairs])}\t")
            output_file.write(
                f"{','.join([residue_pair.first.chain for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([str(residue_pair.first.index) for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([residue_pair.first.name for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([residue_pair.second.chain for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([str(residue_pair.second.index) for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([residue_pair.second.name for residue_pair in consensus.residue_pairs.values()])}\t")
            output_file.write(
                f"{','.join([str(residue_pair.bond or '') for residue_pair in consensus.residue_pairs.values()])}\n")


def consensus_residue_pairs(residue_pairs: dict[Tuple[str, str], list[ResiduePair]],
                            baits: dict[str, AlignmentRecord] = None, targets: dict[str, AlignmentRecord] = None) \
        -> list[ConsensusResiduePair]:
    """
    Creates a list of consensus residue pairs based on individual residue-pairs.

    :param residue_pairs: individual residue-pairs
    :param baits: aligned baits
    :param targets: aligned targets
    :return: list of consensus residue pairs based on individual residue-pairs.
    """
    consensuses = {}
    for bait, target in residue_pairs:
        for residue_pair in residue_pairs[(bait, target)]:
            first_index = residue_pair.first.index
            second_index = residue_pair.second.index
            if baits:
                first_index = baits[bait].position(first_index)
            if targets:
                second_index = targets[target].position(second_index)
            if (first_index, second_index) not in consensuses:
                consensuses[(first_index, second_index)] = ConsensusResiduePair(first_index, second_index)
            consensuses[(first_index, second_index)].append(bait, target, residue_pair)
    return list(consensuses.values())


def parse_residue_pairs(residue_pairs_file: TextIO) -> list[ResiduePair]:
    """
    Parses residue-pairs file.

    :param residue_pairs_file: residue-pairs file
    :return: all residue-pairs
    """
    residue_pairs = []
    residue_pairs_file.readline()
    for line in residue_pairs_file:
        columns = line.rstrip('\r\n').split('\t')
        first = Residue(columns[0], int(columns[1]), columns[2])
        second = Residue(columns[3], int(columns[4]), columns[5])
        residue_pairs.append(ResiduePair(first, second, columns[6] if columns[6] else None))
    return residue_pairs


def parse_alignment(alignment_file: TextIO, alignment_format: str) -> dict[str, AlignmentRecord]:
    """
    Parses alignment file.

    :param alignment_file: alignment file
    :param alignment_format: format alignment file, see https://biopython.org/wiki/AlignIO
    :return: list of AlignmentRecord
    """
    alignment = AlignIO.read(alignment_file, alignment_format)
    sequences = {}
    for record in alignment:
        sequences[record.id] = AlignmentRecord(record)
    return sequences


if __name__ == '__main__':
    main()
